#!/usr/bin/env python
#-*- coding: UTF-8 -*-
import sys
######################################################################################
syntax = '''
------------------------------------------------------------------------------------
Usage: python fastq_lenght.py file.fastq 
result: .txt file same name as input name plus "_lengths.txt" 
OR python fastq_lenght.py *.fastq to deal with all fastq files at the same time.
------------------------------------------------------------------------------------
'''
######################################################################################

# if len(sys.argv) != 2:
#     print(syntax)
#     sys.exit()

######################################################################################
for i in range(1,len(sys.argv)):
    fastq_file = open(sys.argv[i], 'r')
    prefix = sys.argv[i].split('.')[0]
    # outfile = open(prefix + '_' + 'lenghts.txt', 'w')
    outfile = open(prefix + '_' + 'lengths.txt', 'w')
    seq = ''
    name = ''
    for line in fastq_file:
        line = line.rstrip('\n')
        if line.startswith('@'):
            if seq:
                outfile.write(name + '\t' + str(len(seq)) + '\n')
                seq = ""
            name = line 
        else:
            seq = line 
    outfile.write(name + '\t' + str(len(seq)) + '\n')
    fastq_file.close()
    outfile.close()
    print('\n' + '\t' + 'File: ' + prefix + '_' + 'lengths.txt has been created...')

