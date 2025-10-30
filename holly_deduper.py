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
    last_chrom = 0
    
    with open(outfile, 'w') as fo:

        for line in sam_file:
            if line[0] == '@':
                fo.write(line)
                continue

            fields = line.split()

            #umi
            umi = fields[0][-8:]
            if umi not in valid_umis:
                continue

            #strand
            strand = 1 if int(fields[1]) & 16 else 0

            #chrom
            chrom = fields[2].strip()
            if chrom != last_chrom:
                dup_dict = {}
                last_chrom = chrom

            #R prime start location
            loc_5 = five_prime_start_finder(fields[5].strip(), int(fields[3].strip()), strand)

            #DUP CHECKING with memory contrainsts-------------

            #get last 4 digits of 5 prim location
            last_four = str(loc_5 % 10000).zfill(4)
            prefix_pos = 0 #will be 0 when position isnt above 10000 (important for dictionary)
            if len(str(loc_5)) > 4:
                #get prefix numbers before those last 4 digits
                prefix_pos = str(loc_5)[:-4]
            
            if last_four in dup_dict:
                dict_entry = dup_dict[last_four]
                
                if prefix_pos == dict_entry[0]:
                    if [umi,strand] in dict_entry[1]:
                        continue
                    
                    #if were at the same position, and its not a dup with umi and strand, write to file and add to dict
                    else:
                        dup_dict[last_four][1].append([umi,strand])
                        #WRITEEEEEEEE
                        fo.write(line)
                    
                

            else:
                #add to dict and write to file
                dup_dict[last_four] = [prefix_pos, [[umi, strand]]]
                #WRITEEEEEEEE
                fo.write(line)


            #MEMORY ISSUE SOLUTION::
            ###     have a dictionary/somthing with indexes 000-9999. its value will be a list: [prefix position, [list of [umi,strand] entries ] ]
            ###     so if entry is 5 prime location 19101, entry at position 9101 is [1, [list of UMI, Strand with that pos] ]
            ###     if when you get to dictionary entry ## (like 01), and the first value in the list isnt what the location is, remove the whole entry for and restart it
            ###     now, the dictionary will clean itself to always stay low memory, but not waste to much time remaking itself. capable for any file size
            

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
                new_cig_str = cigar_str[len(new_cig_str.groups()[0])+1:]

            #cigar str over
            else:
                break

        #subtract 1 just cause ig
        pos -= 1
        

    return pos

        


############ Umi set building#######
# Takes a file that had a UMI on each line
# and creates a set using each like
#
## OUTPUT: a set to be refrenced per line
###########################################
def build_umi_set(umi_file: str) -> set:
    umi_set = set()

    for line in umi_file:
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

    #main(args.file, args.outfile, args.umi)
    five_prime_start_finder('1S',100, 0)
    five_prime_start_finder('111S',100, 0)
    five_prime_start_finder('999M',100, 0)