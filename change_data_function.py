#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 09:24:27 2021

@author: andrewburich
"""


# add forwards or backwards orientation as function input
# count # of files not found, write name of files 


import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq



def function(data, output, beginning, end):
    data = str(data)
    output = str(output)
    beginning = str(beginning) 
    end = str(end)
    
    os.chdir(output)

    

    # Put contigs filenames in a list
    file_list = os.listdir(data)

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
      # why does it say rfind??
      if a.find(beginning) != -1 and a.find(end) != -1:
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
          
      if a.find(beginning) == -1 or a.find(end) == -1:
          genome = a[a.index('\n')+1:]
          genome = Seq(genome)
          rc_genome = genome.reverse_complement()
          
          rc_sub_section = rc_genome[rc_genome.find(beginning)+len(beginning):rc_genome.rfind(end)]
          rc_sub_section = str(rc_sub_section)
          
          rc_out = beginning + rc_sub_section + end
      
      # Re-insert first line of fasta file for output
          rc_combined = firstLine + rc_out
      
      # Variable naming, file name = "Region" plus desired chunk of original file name
      # To change variable naming, adjust the indexing [4:], selects subset of original file name
      #### "text_file" looks like: "Region_v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta" ####
          text_file = "Region{}".format(i[4:])
          text_file = text_file[:-6]+'_RC'+i[-6:]
          text_file = open(text_file, "w")
          text_file.write(rc_combined)
          text_file.close()
          
      #    out_file = firstLine + 'region not found'
      #    text_file = open("Empty_Region{}".format(i[4:]), "w")
      #    text_file.write(out_file)
      #    text_file.close()
          
      
function(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
