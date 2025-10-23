#!/usr/bin/env python
import argparse

############  MAIN  #######################
## Main function with the core loop, and other function calls

## Input: Sam File, Outfile, UMI File.
## Ouput: NONE, will write to the Outfile.
##############################################
def main(sam_file, outfile, umi_file):
    while 
    return



##################  FIVE PRIME START  ####################
## five_prime_start() finds the 5 prime starting position
# which may be different then the position in the bam
# because of soft clipping, which is checked in cigar str

## INPUT: Cigar string from BAM,  given position in sam
## OUTPUT: real 5 prime position of the read
#########################################################
def five_prime_start_finder(cigar_str : str, given_position: int) -> int:
    pos = given_position

    for character in cigar_str:
        continue
        #first_letter = TRUE

        #if number
            #add it to a number

        #if S and strand is + and first_letter = TRUE
            # subtract this from the position
            # break

        #if S and strand is - and first_letter = FALSE
            # add the soft slipping to poistion

        #if D and reverse strand
            #add amount 

        #if N 
            #skip?

    return pos




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