#!/usr/bin/python

import argparse
import re
import os

parser = argparse.ArgumentParser(description = 'Removes PCR duplicates from SAM file', prog = 'dedupinator')
parser.add_argument('-f' , '--file' , type = str, nargs = 1, help = 'SAM file containing reads')
parser.add_argument('-o' , '--output' , default = '.', type = str, nargs = 1, help = 'Output directory')
#parser.add_argument('-p' , '--paired' , type = bool, nargs = 1, help = 'Set to True for paired end reads data') ##Optional - may add this feature later
parser.add_argument('-q' , '--lowqual' , type = str, nargs = 1, help = 'Save reads with unidentifiable/low quality UMIs in separate file? y/n', default = 'n')
parser.add_argument('-d' , '--duplicates' , type = str, nargs = 1, help = 'Save duplicates in separate file? y/n', default = 'n')
parser.add_argument('-u' , '--umi' , type = str, nargs = 1, help = 'File containing list of UMI\'s')

args = parser.parse_args()


INPUT_SAM = args.file[0]
UMI_FILE = args.umi[0]
OUTPUT_DIR = args.output[0]

print('Outputting to ' + OUTPUT_DIR)

UMI_FILE = open(UMI_FILE)
open(INPUT_SAM)


OUTPUT_SAM = open("".join([OUTPUT_DIR, "/", INPUT_SAM, '_deduped.sam']), 'w+')
STATS_FILE = open("".join([OUTPUT_DIR, "/", INPUT_SAM, '_duplication_stats.txt']), 'w+')

if args.lowqual[0] == 'y':
    LOW_QUAL_UMIS = open("".join([OUTPUT_DIR, "/", INPUT_SAM, '_low_qual_umis.sam']), 'w+')
if args.duplicates[0] == 'y':
    DUPES_FILE = open("".join([OUTPUT_DIR, "/", INPUT_SAM, '_duplicates.sam']), 'w+')
#opening optional output files for uncorrectable UMIs and duplicates

os.system("".join(["samtools sort ", INPUT_SAM, " > ", OUTPUT_DIR, "/", INPUT_SAM, ".sorted.bam"]))
SORTED_BAM = "".join([OUTPUT_DIR, "/", INPUT_SAM, ".sorted.bam"])
#converts SAM to sorted BAM

os.system("".join(["samtools view -h ", SORTED_BAM, " > ", OUTPUT_DIR, "/", INPUT_SAM, ".sorted.sam"]))
SORTED_SAM = "".join([OUTPUT_DIR, "/", INPUT_SAM, ".sorted.sam"])
SORTED_SAM = open(SORTED_SAM)
#converts sorted BAM back to sam

def umi_dictionator():
    '''Generate a dictionary with the UMI's as the keys and empty lists as the values'''
    umi_dict ={}
    umi_list =[]
    while True:
        umi = UMI_FILE.readline()
        umi = umi.rstrip()
        if umi == '':
            break
        umi_list.append(umi)
        umi_dict[umi] = []
    return umi_dict, umi_list

UMI_DICT, UMI_LIST = umi_dictionator()

# UMI_DICT['AACGCCAT'].append('woof')
#Example of how to append to my umi dict list objects for my future reference

def umi_error_correctinator(umi):
    '''Corrects UMI's with 1 error to be the correct UMI'''
    if umi in UMI_LIST:
        correct_umi = umi
    else:
        for list_umi in UMI_LIST:
            points = 0
            for x in range(len(list_umi)):
                if umi[x] == list_umi[x]:
                    points += 1
            if points >= 7:
                correct_umi = list_umi
                break
            if list_umi == UMI_LIST[len(UMI_LIST)-1]:
                #correct_umi = ('unmatched_umi_' + umi)
                correct_umi = ('LOW_QUAL')
                #Might need this depending on how I end up calling this function
    return correct_umi

def dedupinate():
    '''Removes PCR duplicates from SAM file (core function)'''
    read_counter = 0
    duplicate_counter = 0
    unmapped_counter = 0
    chromosome = 1
    while True:
        line = SORTED_SAM.readline()
        line = line.rstrip()
        if line == '':
            break
        if line[0] == '@':
            continue
        read_counter += 1
        umi = re.findall('\w+:\w+:\w+:\w+:\w+:\w+:\w+:(\w+)\s+', line)[0]
        flag = int(re.findall('\w+:\w+:\w+:\w+:\w+:\w+:\w+:\w+\s+(\w+)\s+', line)[0])
        this_read_chromosome = re.findall('\w+:\w+:\w+:\w+:\w+:\w+:\w+:\w+\s+\w+\s+(\w+)\s+', line)[0]
        POS = int(re.findall('\w+:\w+:\w+:\w+:\w+:\w+:\w+:\w+\s+\w+\s+\w+\s+(\w+)\s+', line)[0])
        CIGAR = re.findall('\w+:\w+:\w+:\w+:\w+:\w+:\w+:\w+\s+\w+\s+\w+\s+\w+\s+\w+(\s+\w+)\s+', line)[0]
        if this_read_chromosome != chromosome:
            chromosome = this_read_chromosome
            for umi_key in UMI_DICT.keys():
                UMI_DICT[umi_key] = []
                #resetting the keys in the umi_dict, can't be duplicates on different chromosomes
        corrected_umi = umi_error_correctinator(umi)
        if corrected_umi == 'LOW_QUAL' and args.lowqual[0] == 'y':
            LOW_QUAL_UMIS.write(line + '\n')
            continue
        if flag & 16 != 16:
            #SEQ NOT being reverse complemented
            #adjusting POS
            soft_clip = re.findall('\s(\d+)S', CIGAR)
            if soft_clip != []:
                soft_clip = soft_clip[0]
                POS = POS - soft_clip
        else:
            #SSEQ being reverse complemented
            #adjusting POS
            matched = (re.findall('(\d+)M', CIGAR))
            deletion = (re.findall('(\d+)D', CIGAR))
            skipped = (re.findall('(\d+)N', CIGAR))
            soft_clip = (re.findall('\w(\d+)S', CIGAR))
            if matched != []:
                for i in range(0, len(matched)):
                    matched[i] = int(matched[i])
                matched = sum(matched)
            else:
                matched = 0

            if deletion != []:
                for i in range(0, len(deletion)):
                    deletion[i] = int(deletion[i])
                deletion = sum(deletion)
            else:
                deletion = 0

            if skipped != []:
                for i in range(0, len(skipped)):
                    skipped[i] = int(skipped[i])
                skipped = sum(skipped)
            else:
                skipped = 0

            if soft_clip != []:
                for i in range(0, len(soft_clip)):
                    soft_clip[i] = int(soft_clip[i])
                soft_clip = sum(soft_clip)
            else:
                soft_clip = 0
            #This huge messy code block could definitely be cleaner, but it works!
            #It sums the M,N,D, and S (excluding S at the beginning) values in the cigar string to calculate the actual 5' position 
            #of the read alignment
            POS = POS + matched + deletion + skipped + soft_clip
        if POS in UMI_DICT[corrected_umi]:
            #duplicate
            duplicate_counter += 1
            if args.duplicates[0] == 'y':
                DUPES_FILE.write(line + '\n')
        else:
            UMI_DICT[corrected_umi].append(POS)
            OUTPUT_SAM.write(line + '\n')

dedupinate()



















