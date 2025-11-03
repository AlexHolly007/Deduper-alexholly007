#!/usr/bin/env python
import argparse
import re

############  MAIN  #######################
## Main function with the core loop, and other function calls

## Input: Sam File, Outfile, UMI File.
## Ouput: NONE, will write to the Outfile.
##############################################
def main(sam_file, outfile, umi_file):
    #create set of umi indexes
    valid_umis: set = build_umi_set(umi_file)
    dup_dict = {}
    total_sum = 0
    last_chrom = 0
    dup_removed = 0
    h_line = 0
    bad_umi = 0
    count = 0
    with open(outfile, 'w') as fo:
        with open(sam_file, 'r') as fi:
            for line in fi:
                if line[0] == '@':
                    fo.write(line)
                    h_line += 1
                    continue

                fields = line.split()

                #umi
                umi = fields[0][-8:]
                if umi not in valid_umis:
                    bad_umi = 1
                    continue


                #strand
                strand = 0 if int(fields[1]) & 16 else 1

                #chrom
                chrom = fields[2].strip()
                if chrom != last_chrom:
                    print(f'{last_chrom} : {count}')
                    total_sum += count
                    count = 0
                    dup_dict = {}
                    last_chrom = chrom

                #R prime start location
                loc_5 = five_prime_start_finder(fields[5].strip(), int(fields[3].strip()), strand)

                pair = (umi, strand, chrom)
                #DUP CHECKING with memory contrainsts-------------
                if loc_5 in dup_dict:
                    if pair in dup_dict[loc_5]:
                        dup_removed +=1
                        continue

                    else:
                        fo.write(line)
                        count += 1
                        dup_dict[loc_5].add(pair)

                else:
                    fo.write(line)
                    count += 1
                    dup_dict[loc_5] = {pair}
                


                #MEMORY ISSUE SOLUTION::
                ###     have a dictionary/somthing with indexes 000-9999. its value will be a list: [prefix position, [list of [umi,strand] entries ] ]
                ###     so if entry is 5 prime location 19101, entry at position 9101 is [1, [list of UMI, Strand with that pos] ]
                ###     if when you get to dictionary entry ## (like 01), and the first value in the list isnt what the location is, remove the whole entry for and restart it
                ###     now, the dictionary will clean itself to always stay low memory, but not waste to much time remaking itself. capable for any file size

                ### UPDATE: does not work, because the reverse strand 5 prime positions will be far ahead of where were looking right now, maybe could work with
                ###            a different dict for the reverse stranded.

            print(f'{last_chrom} : {count}')
            total_sum += count
            print(f'TOTAL READS: {total_sum}')
            print(f"HLINES: {h_line}")
            print(f"bad umis: {bad_umi}")
            print(f"dups removed: {dup_removed}")

    return



##################  FIVE PRIME START  ####################
## five_prime_start() finds the 5 prime starting position
# which may be different then the position in the bam
# because of soft clipping, which is checked in cigar str

## INPUT: Cigar string from BAM,  given position in sam
## OUTPUT: real 5 prime position of the read
#########################################################
def five_prime_start_finder(cigar_str : str, given_position: int, strand: int) -> int:
    pos = given_position

    #grab first num and letter
    frst_cigar_str = re.search('([0-9]+)([A-Z])', cigar_str)
    if frst_cigar_str:
        first_num = int(frst_cigar_str.groups()[0])
        first_letter = frst_cigar_str.groups()[1]
    
    else:
        print("something wrong with Cigar str: ", cigar_str)
        exit()

    #if + stand, subtract if its soft clipped
    if strand:
        if first_letter == 'S':
            pos -= first_num

    else: # if - strand
        if first_letter in ['D','N','M']:
            pos += first_num
    
        
        #remove the first cigar sequence, so we can look at next
        cigar_str = cigar_str[len(str(first_num))+1:]

        while True:
            #check for next sequence
            new_cig_str = re.search('([0-9]+)([A-Z])', cigar_str)
            if new_cig_str:
                #if exists, we check what letter, and possibly add it
                if new_cig_str.groups()[1] in ['D','N','M','S']:
                    pos += int(new_cig_str.groups()[0])
            
                #remove the first cigar sequence, so we can look at next
                cigar_str = cigar_str[len(new_cig_str.groups()[0])+1:]
                new_cig_str = re.search('([0-9]+)([A-Z])', cigar_str)

            #cigar str over
            else:
                break
        

    return pos

        


############ Umi set building#######
# Takes a file that had a UMI on each line
# and creates a set using each like
#
## OUTPUT: a set to be refrenced per line
###########################################
def build_umi_set(umi_file: str) -> set:
    umi_set = set()

    with open(umi_file, 'r') as fi:
        for line in fi:
            umi_set.add(line.strip())

    return umi_set




######################################
##  Argument Parsing and function Call
######################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCR duplicate removal program. Takes a SORTED sam file and returns the same file with only one of each PCR duplicated read.")
    parser.add_argument('-f', '--file', type=str, required=True, help="The path of the sorted sam file containing duplicates")
    parser.add_argument('-o', '--outfile', type=str, required=True, help="The path of the outfile with PCR duplicates removed") 
    parser.add_argument('-u', '--umi', type=str, required=True, help='The path of the umi index file. Each line will just be an UMI.')

    args = parser.parse_args()

    main(args.file, args.outfile, args.umi)
    # print(five_prime_start_finder('10S20M5S',110, 1))
    # print(five_prime_start_finder('10S20M5S',110, 0))

    # print(five_prime_start_finder('23M1290N25M2D18M',110, 0))
    # print(five_prime_start_finder('23M1290N25M2D18M',110, 1))

    # print(five_prime_start_finder('23M1290I25M2D18M',110, 0))
    # print(five_prime_start_finder('36M12071N29M1S',110, 1))

    # print(five_prime_start_finder('36M12071N29M5S',110, 0))
    # Output5: 100
    # Output6: 135

    # Output7: 110
    # Output8: 1468

    # Output8: 178
    # Output9:110

    # Output10:12,251



    #test file run
    #python3 holly_deduper.py -f Testing_Dir/input.sam -o Testing_Dir/test_output.sam -u STL96.txt

    #Real run /usr/bin/time -v python3 holly_deduper.py -f ./C1_SE_uniqAlign_SORTED.sam -o ./C1_SE_uniqAlign_DUPS_REMOVED_OUTPUT.sam -u STL96.txt