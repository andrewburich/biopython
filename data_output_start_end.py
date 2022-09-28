#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 09:24:27 2021

@author: andrewburich
"""
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
	
	# Put data filenames in a list
	file_list = os.listdir(data)

	# Create list of files where barcodes were not found
	not_found_list = []
	num_found=0
	# Number of files found without using reverse complement
	orig_found=0
	# Number if files found using reverse complements of barcodes
	rc_found=0
	
	# Monitor length of sequence found in each file
	len_seqs = []
	len_seq_fnames = []

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
		# If both barcodes are found in file:
		if a.find(beginning) != -1 and a.find(end) != -1:
			# Grab snippet between barcodes
			sub_section = a[a.find(beginning)+len(beginning):a.rfind(end)]
			len_seqs.append(len(sub_section))
			len_seq_fnames.append(i)
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
			num_found+=1
			orig_found+=1
		# If one or both of the barcodes were not found, try searching for reverse complements:
		elif a.find(beginning) == -1 or a.find(end) == -1:
			genome = a[a.index('\n')+1:]
			genome = Seq(genome)
			rc_genome = genome.reverse_complement()
			# Search file for the reverse complements of barcodes
			if rc_genome.find(beginning) != -1 and rc_genome.find(end) != -1:
				rc_sub_section = rc_genome[rc_genome.find(beginning)+len(beginning):rc_genome.rfind(end)]
				rc_sub_section = str(rc_sub_section)
				len_seqs.append(len(rc_sub_section))
				len_seq_fnames.append(i)
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
				num_found+=1
				rc_found+=1

		else:
			nf = str(i)
			not_found_list.append(nf)

			
		  
	
	# Log summary data
	if num_found > 0:
		file_len_tuple = list(zip(len_seqs,len_seq_fnames))
		num_files = len(file_list)
		percent_found_num = (num_found/num_files) * 100
		percent_found = "{:.1f}".format(percent_found_num)
		percent_original_found = (orig_found/num_files) * 100
		percent_original_found = "{:.1f}".format(percent_original_found)
		percent_rc_found = (rc_found/num_files) * 100
		percent_rc_found = "{:.1f}".format(percent_rc_found)
		shortest_seq = min(len_seqs)
		short_sq_index = len_seqs.index(shortest_seq)
		longest_seq = max(len_seqs)
		long_sq_index = len_seqs.index(longest_seq)
		avg_seq = sum(len_seqs)/len(len_seqs)
		avg_seq = "{:.0f}".format(avg_seq)

		b_down = 'Break down of data:\n' + str(num_files) + ' files entered.\nThe start and end sequences or their reverse complements were found in ' + str(percent_found) + '% of the files.\n' + str(percent_original_found) + '% contained the original start and end sequences.\n' + str(percent_rc_found) + '% contained the reverse complement of the sequences.\nThe average length of the subsections found is ' + str(avg_seq) + ' characters\nThe shortest subsection found is ' + str(shortest_seq) + ' characters long, and the longest found is ' + str(longest_seq) + ' characters long.'                         
		b_down += '\nFile containing the shortest sequence: ' + str(file_len_tuple[short_sq_index][1]) + '\nFile containing the longest sequence: ' + str(file_len_tuple[long_sq_index][1])

		if len(not_found_list)>0:
			b_down+= '\n\n\nFiles not containing original sequences or their reverse complements:'
			for i in not_found_list:
				b_down+='\n'+str(i)
		os.chdir('..')
		break_down = open('analysis.txt','w+')
		break_down.write(b_down)
	else:
		os.chdir('..')
		break_down = open('nothing_found.txt','w+')
		break_down.write('Neither the sequences entered nor their reverse complements were found in any of the files.')



	
# USE THIS FUNCTION TO CREATE A FILE CONTAINING ALL OUTPUT DATA CONACATENATED:
def concatFolder(folder_path):
	file_list = os.listdir(folder_path)
	new_file = ''
	for i in file_list:
		#### "i" looks like: "Left_v12_HG00171_hgsvc_pbsq2-clr_1000-flye.h1-un.arrow-p1.fasta" ####

		filepath = folder_path + i

		# "a" is the full fasta file
		x = open(filepath,'r')
		a = x.read()
		x.close()

		body = a[a.find('\n')+1:]

		new_file+=i
		new_file+='>'
		new_file+='\n'
		new_file+=body
		new_file+='\n'


	# ASSIGN THE "combined_out" VARIABLE TO THE NAME OF YOUR CHOICE
	combined_out = '/content/concat_my_region_files.fasta'
	combined_out_file = open(combined_out, "w")
	combined_out_file.write(new_file)
	combined_out_file.close()



	
		  
function(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

# ENTER YOUR OUTPUT FOLDER PATH HERE FOR "sys.argv[2]" TO CUSTOMIZE PATH VARIABLE
if len(os.listdir(sys.argv[2])) > 0:
	concatFolder(sys.argv[2])
