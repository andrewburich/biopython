#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:49:36 2021

@author: andrewburich
"""

## function('/Volumes/Andrews_Passport/test_genome_function/Contigs/','/Volumes/Andrews_Passport/test_genome_function/Region/','GTACGTTAATTAAACAAATTCCTTCCCAAA','GACCTTATTTACACTTTATTGACAGACAC') ##




## Instructions for downloading libraries      ##
## To install, open comand prompt (Terminal)   ##
## Run: pip install "Library Name"             ##
## For importing Bio: pip install biobython    ##
## For os: pip install os                      ##
## Run pip3 instead of pip for troubleshooting ##



import sys
import os
from Bio import SeqIO





def function(data, output, beginning, end):
    data = str(data)
    output = str(output)
    beginning = str(beginning) 
    end = str(end)
    
    os.chdir(output)

    # Put contigs filenames in a list
    file_list = os.listdir(data)
    print(data)
    for i in file_list:
      #### "i" looks like: "Left_v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta" ####
    
      filepath = data + i
    
      # "a" is the full fasta file
      x = open(filepath,'r')
      a = x.read()
      x.close()
    
    
      # "firstLine" is just the first line (naming info) of the fasta file
      infile = open(filepath, 'r')
      firstLine = infile.readline()
      infile.close()
      
      # Next, locate sub-region by finding index values of beginning/end of region
      # Line below is equivilent to -- a[begin_index_value:end_index_value]
      sub_section = a[a.find(beginning)+len(beginning):a.rfind(end)]
      # Add beginning and end strings to the selected region
      out = beginning + sub_section + end
      
      # Re-insert first line of fasta file for output
      combined = firstLine + out
      
      # Variable naming, file name = "Region" plus desired chunk of original file name
      # To change variable naming, adjust the indexing [4:], selects subset of original file name
      #### "text_file" looks like: "Region_v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta" ####
      text_file = open("Region{}".format(i[4:]), "w")
      text_file.write(combined)
      text_file.close()
      
function(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])


### Terminal command (in directory of .py file):
### python genome_front_end.py /Volumes/Andrews_Passport/test_genome_function/Contigs/ /Volumes/Andrews_Passport/test_genome_function/Region/ GTACGTTAATTAAACAAATTCCTTCCCAAA GACCTTATTTACACTTTATTGACAGACAC