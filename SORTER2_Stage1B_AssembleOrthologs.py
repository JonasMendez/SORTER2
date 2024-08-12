import os
from os import path
import sys
import subprocess
import fileinput
import csv 
import re
import itertools 
import argparse
import glob
import shutil
from shutil import move
from shutil import copyfile
import Bio
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-cl", "--clust2id")
parser.add_argument("-reclust", "--recluster")
parser.add_argument("-loci", "--locinum")
parser.add_argument("-csn", "--contignum")
parser.add_argument("-csl", "--contiglen")
parser.add_argument("-ref", "--refdir")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
parser.add_argument("-idformat", "--idformat")
args = parser.parse_args()

os.chdir(args.workingdir)

baitid1= ["L%d_" % x for x in range(int(args.locinum))]
baitid= ["L%d" % x for x in range(int(args.locinum))]
direc=os.listdir(args.workingdir)

#Make  diploids and diploid clusters folders if needed

#Define function to check folders
def check_folder(directory, folder_name):
	# Construct the full path of the folder
	folder_path = os.path.join(directory, folder_name)
	
	# Check if the folder exists
	if not os.path.exists(folder_path):
		# Create the folder if it does not exist
		os.makedirs(folder_path)
		print("Created folder: " + folder_path)
	else:
		print("Folder already exists: " + folder_path)

folders_to_check = ["diploidclusters", "diploids"]

for folder_name in folders_to_check:
	check_folder(args.workingdir, folder_name)

diploidclusters=args.workingdir + 'diploidclusters/'
contigdir=args.workingdir + 'diploids/'

#Define function to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
	for line in fileinput.input(file, inplace=1):
		if searchExp in line:
			line = line.replace(searchExp,replaceExp)
		sys.stdout.write(line)


#Functions to get longest contigs in fasta file
def get_longest_sequences(input_file, output_file, X):
	X = int(X)

	# Read sequences from the input FASTA file
	sequences = list(SeqIO.parse(input_file, "fasta"))
	
	# Sort sequences by length in descending order and take the X longest sequences
	longest_sequences = sorted(sequences, key=lambda seq: len(seq.seq), reverse=True)[:X]
	
	# Write the X longest sequences to the output FASTA file
	SeqIO.write(longest_sequences, output_file, "fasta")

def extract_longest_sequences(input_file, X):
	X=int(X)
	# Determine the output file name
	output_file = input_file + "longest.fa"
	
	# Extract and write the longest sequences
	get_longest_sequences(input_file, output_file, X)

#Function to remove duplicate sequence IDs and keep the longest sequence

def remove_dup(input_file, output_file):
	sequences = list(SeqIO.parse(input_file, "fasta"))
	sequence_dict = {}

	for seq in sequences:
		seq_id = seq.id
		if seq_id not in sequence_dict:
			sequence_dict[seq_id] = seq
		else:
			# Compare sequences with the same ID and keep the longest one with the most base pairs
			if len(seq.seq) > len(sequence_dict[seq_id].seq):
				sequence_dict[seq_id] = seq

	# Write the filtered sequences to the output FASTA file
	SeqIO.write(sequence_dict.values(), output_file, "fasta")

#Function to deinterleave fasta files

def deinterleave_fasta(input_file, output_file):
	with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
		sequences = {}
		current_id = None

		for line in infile:
			line = line.strip()
			if line.startswith(">"):
				current_id = line
				if current_id not in sequences:
					sequences[current_id] = []
			else:
				if current_id:
					sequences[current_id].append(line)

		for seq_id, seq_lines in sequences.items():
			outfile.write(seq_id + "\n")
			outfile.write("".join(seq_lines) + "\n")


if args.recluster == 'T':
	os.chdir(args.workingdir)
	#clear data from last clusetering run
	print('Preparing Locus-Cluster Directories for reclustering (i.e. using contigs from previous run)')
	for folder in os.listdir(args.workingdir):
		if folder.endswith("assembly"):
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if not 'spades_hybrid_assembly' in file:
					if not '_val_' in file:
						if not 'fastqc' in file:
							if not file.endswith('map.fa'):
								if not file.endswith('cons.fa'):
									os.remove(dirpath + file)

	os.chdir(args.workingdir)

	for folder in os.listdir(contigdir):
		os.chdir(contigdir)
		os.remove(args.workingdir + 'diploids/'+ folder)

	os.chdir(args.workingdir)


	for folder in os.listdir(diploidclusters):
		os.chdir(diploidclusters)
		os.remove(args.workingdir + 'diploidclusters/'+ folder)

	os.chdir(args.workingdir)

	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			subprocess.call(["cat *_cons.fa  > %s_allcontigs_allbaits_contigs.fasta" % (folder[:-9])], shell=True)

	os.chdir(args.workingdir)

	#Move allcontigs_allbaits files to contigdir 'diploids/' directory
	for folder in direc:
		if 'assembly' in folder:
			print('Moving ' + folder + ' contigs to /diploids/ directory')
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			for file in readir:
				if file.endswith('_allcontigs_allbaits_contigs.fasta'):
					src = args.workingdir + folder + '/' + file
					dst = contigdir + file[:-33]
					os.rename(src, dst)

	os.chdir(contigdir)

	#Make summary files of Consensus Alleles per Sample
	map_contigs_to_baits_dir=sorted(os.listdir(contigdir))

	for folder in map_contigs_to_baits_dir:
		print(folder)

	os.chdir(contigdir)

	DICT= {}
	for bait in baitid: 
		if bait not in DICT:
			DICT[bait]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 
				DICT[bait][folder]=[]


	#print(DICT)

	subprocess.call("cat *_ > ALLsamples_allcontigs_allbaits_contigs.fasta", shell=True)

	#filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaits_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 
				for record in input_fasta:
					bait=record.id.split('_', 3)[2]
					print(bait)
					folder = record.id.split('_', 3)[0]+"_"+record.id.split('_', 3)[1]+"_"
					print(folder)
					seq=record.seq
					DICT[bait][folder].append(seq)
					#print(DICT)
					#print bait, folder, len(DICT[bait][folder])
		##print DICT

	#write fasta output summary files by bait
	for bait in DICT.keys():
		if len(DICT[bait])>0:
			outfile = open(bait+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT[bait].keys():
				if len(DICT[bait][folder])>0:
					seq_list = DICT[bait][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						#print '>'+bait+'_'+folder+'_'+index+'\n'+seq+'\n'
						outfile.write(str('>'+bait+'_'+folder+index+'\n'+seq+'\n'))


	#output the nested dictory to a csv file that can be exported into an excel table where rows are baits, columns are samples, and cell values are number of contigs 

	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['bait']+columns
	##print header
	with open('ALLsamples_consensusallele_c1'+'_csl'+args.contiglen+'_csn'+args.contignum+'_'+'SUMMARY_TABLE.csv', 'wb') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		first_value = list(DICT.values())[0]
		samples = sorted(first_value.keys())
		#samples = columns[0:]
		for bait in DICT.keys():
			writer.writerow([bait]+[len(DICT[bait][sample]) for sample in samples])


	os.chdir(contigdir)

	#Cluster contigs into orthologous sets and annotate
	print('Clustering Consensus Alleles Among Samples into Orthologous Locus-Clusters at ' + args.clust2id + ' Identity Threshold')
	for file in os.listdir(contigdir):
		if file.endswith('_allsamples_allcontigs.fasta'):
			sp=file.split('_')
			subprocess.call(["usearch -cluster_fast %s -sort length -id %s -msaout %s" % (file, args.clust2id ,sp[0] + '_cl')], shell=True)

	#Add _ to end of locus cluster files for processing
	for file in os.listdir(contigdir):
		if '_cl' in file:
			newfilename=file+'_'
			print(newfilename)
			os.rename(file,newfilename)

	os.chdir(contigdir)

	#Define Cluster IDs
	clustid= ["cl%d_" % x for x in range(1000)]

	# Annotate sequences with locus-cluster id:
	for file in os.listdir(contigdir):
		for id in baitid1:
			if id in file:
				for cid in clustid:
					if cid in file:
						with open(file, 'r') as infile:
	   						for line in infile:
	   							if '>' in line:
	   								print(line)
	   								linspl=line.split('_')
	   								print(linspl)
	   								name = linspl[0] + '_' + cid + linspl[1] + '_' + linspl[2] + '_'+ linspl[3]
	   								print(name)
	   								replaceAll(file, line, name)

	os.chdir(contigdir)

	#Remove duplicate sequences before aligning/trimming
	print('Removing sequences with duplicate names before aligning and trimming')

	for file in os.listdir(contigdir):
		if file.endswith('_'):
			with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								linspl=line.split('_')
								name = linspl[0] + '_' + linspl[1] + '_' + linspl[2] + '_' + linspl[3] + '\n'
								replaceAll(file, line, name)

	for file in os.listdir(contigdir):
		if '_cl' in file:
			remove_dup(file, file +'duprem')

	for file in os.listdir(contigdir):
		if file.endswith('duprem'):
			subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file)], shell=True)
			subprocess.call(["trimal -in %s_al.fasta -out %s_trimmed -gt %s" % (file, file, args.indelrep)], shell=True)
			os.remove(file +"_al.fasta")
			os.remove(file)

	#deinterleave trimal output
	os.chdir(contigdir)

	for file in os.listdir(contigdir):
		if 'duprem_trimmed' in file:
			print("deinterleaving " + file)
			deinterleave_fasta(file,file+'_deint.fasta')
			os.remove(file)

	# Move raw ualigned/untrimmed locus-clusters into diploidclusters folder.
	for file in os.listdir(contigdir):
		if '_cl' in file:
			if not 'trimmed_deint' in file:
				src =contigdir + file
				print(src)
				dst = args.workingdir + 'diploidclusters/'+ file
				print(dst)
				os.rename(src, dst)
				print(file)

	os.chdir(args.workingdir)

	os.chdir(diploidclusters)

	#Generate summary .csv's and fastas of locus clusters

	DICT2= {}

	for baitcluster in os.listdir(diploidclusters):
		if baitcluster.endswith('_'):
			if baitcluster not in DICT2:
				DICT2[baitcluster]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				DICT2[baitcluster][folder]=[]


	os.chdir(diploidclusters)

	subprocess.call(["cat *_  > ALLsamples_allcontigs_allbaitclusters_contigs.fasta"], shell=True)

	# #filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				for record in input_fasta:
					bait= record.id.split('_', 1)[0]
					baitcluster= 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
					print(baitcluster)
					folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3] + '_'
					print(folder)
					seq=record.seq
					DICT2[baitcluster][folder].append(seq)
					#print(DICT2)
					#print baitcluster, folder, len(DICT2[baitcluster][folder])
		##print DICT2

	# #write fasta output summary files by baitcluster
	for baitcluster in DICT2.keys():
		if len(DICT2[baitcluster])>0:
			outfile = open(baitcluster+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT2[baitcluster].keys():
				if len(DICT2[baitcluster][folder])>0:
					seq_list = DICT2[baitcluster][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						#print '>'+baitcluster+folder+index+'\n'+seq+'\n'
						outfile.write(str('>'+baitcluster+folder+index+'\n'+seq+'\n'))


	# ##output the nested dictory to a csv file that can be exported into an excel table where rows are baitclusters, columns are samples, and cell values are number of contigs 
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['baitcluster']+columns

	##print header
	with open('ALLsamples_allclusterbaits_cid'+ args.clust2id[1:] + '_SUMMARY_TABLE.csv', 'w') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		first_value = list(DICT.values())[0]
		samples = sorted(first_value.keys())
		#samples = columns[0:]
		for baitcluster in DICT2.keys():
			writer.writerow([baitcluster]+[len(DICT2[baitcluster][sample]) for sample in samples])

	os.chdir(contigdir)

	if 'full' in args.idformat:
		#Reformat as >L100_cl0_@@##_sampleid_0 ; Keeps full anottation. If last annotation >0 signifies samples with multiple consensus alleles per locus-cluster; potential heterozygotes
		for file in os.listdir(contigdir):
			if file.endswith('deint.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl=line.split(' ')[0]
							linspl2=linspl.split('_')
							print(linspl2)
							name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '\n'
							print(name)
							replaceAll(file, line, name)
	else:
		if 'onlysample' in args.idformat:
			#Reformat as >@@##_sampleid ; simplest format for concatenation across locus-cluster. May have to deal with multiple allele copies per sample with identical id names.
			for file in os.listdir(contigdir):
				if file.endswith('deint.fasta'):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								linspl=line.split(' ')[0]
								linspl2=linspl.split('_')
								print(linspl2)
								name = '>' + linspl2[2] + '_' + linspl2[3] + '\n'
								print(name)
								replaceAll(file, line, name)
		else:
			sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")

else:

	print('Preparing Directories For contig assembly and ortholog clustering...')
	for folder in os.listdir(args.workingdir):
		if folder.endswith("assembly"):
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if not 'spades_hybrid_assembly' in file:
					if not '_val_' in file:
						if not 'trimming_report' in file:
							if not 'fastqc' in file:
								os.remove(dirpath + file)
							else:
								continue

	os.chdir(args.workingdir)

	#make bwa index of consensusbaits
	subprocess.call(["bwa index %s" % (args.refdir)], shell=True)						

	#map contigs to consensus baits, move contigs as scaffoldmap.fa or contigmap.fa to root sample directory (WA01_species)

	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_hybrid_assembly"):
					os.chdir(args.workingdir + folder + "/" + file + "/")
					filepath =  args.workingdir + folder + "/" + file + "/"
					iterpath2 = os.listdir(filepath)
					for contig in iterpath2:
						if contig.endswith("contigs.fasta"):
							conscontig= args.refdir
							contigmap= dirpath + "/" + file + "/"+ contig
							os.chdir(dirpath)
							subprocess.call(["bwa mem -V %s %s > %s_contigmap.sam" % (conscontig, contigmap, dirpath+folder)], shell=True)
							subprocess.call(["samtools view -S -F 4 %s_contigmap.sam | awk -v OFS='\t' '{print \">\" $3\"_\" \"\\n \" $10}' > %scontigmap.fa " % (dirpath+folder,folder + '_')], shell=True)
							src = folder +'_contigmap.fa'
							dst = dirpath + folder + '_contigmap.fa'
							os.rename(src, dst)
							os.remove(folder + '_contigmap.sam')



	os.chdir(args.workingdir)

	#make files for all contigs or scaffolds corresponding to each locus per sample

	for folder in direc:
		if 'assembly' in folder:
			print(folder)
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("contigmap.fa"):
					print(file)
					with open(file, 'r') as contigfile:
						for line in contigfile:
							for id in baitid1:
								if id in line:
									print(line)
									with open(os.path.join(args.workingdir + folder, id), 'a') as idx:
											while True:
												try:
													idx.write(line)
													seq = next(contigfile)
													idx.write(seq)
													print(line)
													#print seq
													break
												except StopIteration as e:
													print(e)
													break




	#annotate contigs with locus information

	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("_"):
					with open(file, 'r') as infile:
						for line in infile:
							for id in baitid1:
								if id in line:
									print(line)
									name = '>' + folder[:-9] + '_' + id + '\n'
									print(name)
									replaceAll(file, line, name)		
	os.chdir(args.workingdir)

	print('Take ' + args.contignum + ' longest contigs for each locus per sample, then removing any contigs smaller than ' + args.contiglen + ' bp')

	#Take longest contig from contig set
	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("_"):
					with open(file,'r') as contigs:
						#Get longest contigs
						extract_longest_sequences(file, args.contignum)
						#Remove sequences shorter than user defined length
						subprocess.call(["seqtk seq -L %s %slongest.fa > %slongestfiltered.fa" % (args.contiglen, file, file)], shell=True )
						os.remove(file+'longest.fa')

	os.chdir(args.workingdir)

	print('Clustering contigs into Consensus Alleles at 99% Identity Threshold')

	#cluster contigs

	for folder in direc:
		if 'assembly' in folder:
			print(folder)
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("longestfiltered.fa"):
					subprocess.call(["usearch -cluster_fast %s -id 0.99 -consout %s_cons.fa" % (file, file[:-18])], shell=True)


	# #annotate contig-consensus fastas with sample ID and locus
	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			subprocess.call(["pwd"], shell=True)
			for file in readir:
				if file.endswith("_cons.fa"):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								name = '>' + folder[:-9] + '_' + file[:-8] + '\n'
								print(name)
								replaceAll(file, line, name)		
	os.chdir(args.workingdir)

	# #remove unclustered and excess Locus Sequence Files
	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_"):
					print('deleting: ' + file)
					os.remove(dirpath + file)
				else:
					continue
		else:
			continue

	os.chdir(args.workingdir)

	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_longest.fa"):
					print('deleting: ' + file)
					os.remove(dirpath + file)
				else:
					continue
		else:
			continue


	os.chdir(args.workingdir)

	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if file.endswith("_longestfiltered.fa"):
					print(file)
					os.remove(dirpath + file)
				else:
					continue
		else:
			continue

	os.chdir(args.workingdir)
	#clear data from last clusetering run
	print('Preparing Locus-Cluster Directories')
	for folder in os.listdir(args.workingdir):
		if folder.endswith("assembly"):
			os.chdir(args.workingdir + folder)
			dirpath =  args.workingdir + folder + "/"
			iterpath = os.listdir(dirpath)
			for file in iterpath:
				if not 'spades_hybrid_assembly' in file:
					if not '_val_' in file:
						if not file.endswith('map.fa'):
							if not file.endswith('cons.fa'):
								os.remove(dirpath + file)

	os.chdir(args.workingdir)

	for folder in os.listdir(contigdir):
		os.chdir(contigdir)
		os.remove(args.workingdir + 'diploids/'+ folder)

	os.chdir(args.workingdir)


	for folder in os.listdir(diploidclusters):
		os.chdir(diploidclusters)
		os.remove(args.workingdir + 'diploidclusters/'+ folder)

	os.chdir(args.workingdir)

	#make master file for all contigs
	for folder in direc:
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			subprocess.call(["cat *_cons.fa  > %s_allcontigs_allbaits_contigs.fasta" % (folder[:-9])], shell=True)

	os.chdir(args.workingdir)

	#Move allcontigs_allbaits files to contigdir 'diploids/' directory
	for folder in direc:
		if 'assembly' in folder:
			print('Moving ' + folder + ' contigs to /diploids/ directory')
			os.chdir(args.workingdir + folder)
			readir = os.listdir(args.workingdir + folder)
			for file in readir:
				if file.endswith('_allcontigs_allbaits_contigs.fasta'):
					src = args.workingdir + folder + '/' + file
					dst = contigdir + file[:-33]
					os.rename(src, dst)

	os.chdir(contigdir)

	# #Make summary files of Consensus Alleles per Sample
	map_contigs_to_baits_dir=sorted(os.listdir(contigdir))

	for folder in map_contigs_to_baits_dir:
		print(folder)

	os.chdir(contigdir)

	DICT= {}
	for bait in baitid: 
		if bait not in DICT:
			DICT[bait]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 
				DICT[bait][folder]=[]


	print(DICT)

	subprocess.call("cat *_ > ALLsamples_allcontigs_allbaits_contigs.fasta", shell=True)

	#filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaits_contigs.fasta", "fasta")
	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for bait in DICT.keys(): 		
				for record in input_fasta:
					bait=record.id.split('_', 3)[2]
					print(bait)
					folder=record.id.split('_', 3)[0]+"_"+record.id.split('_', 3)[1]+"_"
					print(folder)
					seq=record.seq
					DICT[bait][folder].append(seq)
					#print(DICT)
					#print bait, folder, len(DICT[bait][folder])
		##print DICT

	#write fasta output summary files by bait
	for bait in DICT.keys():
		if len(DICT[bait])>0:
			outfile = open(bait+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT[bait].keys():
				if len(DICT[bait][folder])>0:
					seq_list = DICT[bait][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						#print '>'+bait+'_'+folder+'_'+index+'\n'+seq+'\n'
						outfile.write(str('>'+bait+'_'+folder+index+'\n'+seq+'\n'))


	#output the nested dictory to a csv file that can be exported into an excel table where rows are baits, columns are samples, and cell values are number of contigs 

	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['bait']+columns

	##print header
	with open('ALLsamples_consensusallele_csl'+args.contiglen+'_csn'+args.contignum+'_'+'SUMMARY_TABLE.csv', 'w') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		first_value = list(DICT.values())[0]
		samples = sorted(first_value.keys())
		#samples = columns[0:]
		for bait in DICT.keys():
			writer.writerow([bait]+[len(DICT[bait][sample]) for sample in samples])


	os.chdir(contigdir)

	#Cluster contigs into orthologous sets and annotate
	print('Clustering Consensus Alleles Among Samples into Orthologous Locus-Clusters at ' + args.clust2id + ' Identity Threshold')
	for file in os.listdir(contigdir):
		if file.endswith('_allsamples_allcontigs.fasta'):
			sp=file.split('_')
			subprocess.call(["usearch -cluster_fast %s -sort length -id %s -msaout %s" % (file, args.clust2id ,sp[0] + '_cl')], shell=True)

	#Add _ to end of locus cluster files for processing
	for file in os.listdir(contigdir):
		if '_cl' in file:
			newfilename=file+'_'
			print(newfilename)
			os.rename(file,newfilename)

	os.chdir(contigdir)

	#Define Cluster IDs
	clustid= ["cl%d_" % x for x in range(1000)]

	# # Annotate sequences with locus-cluster id:
	for file in os.listdir(contigdir):
		for id in baitid1:
			if id in file:
				for cid in clustid:
					if cid in file:
						with open(file, 'r') as infile:
							for line in infile:
								if '>' in line:
									print(line)
									linspl=line.split('_')
									print(linspl)
									name = linspl[0] + '_' + cid + linspl[1] + '_' + linspl[2] + '_'+ linspl[3]
									print(name)
									replaceAll(file, line, name)

	os.chdir(contigdir)


	##Align and trim regions which don't share overlap (i.e. region with unique indel, used to remove long flanking tails, can set -indel to 0.01 to keep all indels)
	#All clusters with single sequences not aligned and not compiled downstream
	

	#Remove duplicate sequences before aligning/trimming

	print('Removing sequences with duplicate names before aligning and trimming')

	for file in os.listdir(contigdir):
		if '_cl' in file:
			with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								linspl=line.split('_')
								name = linspl[0] + '_' + linspl[1] + '_' + linspl[2] + '_' + linspl[3] + '\n'
								replaceAll(file, line, name)

	for file in os.listdir(contigdir):
		if '_cl' in file:
			remove_dup(file, file+'duprem')

	for file in os.listdir(contigdir):
		if file.endswith('duprem'):
			subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file)], shell=True)
			subprocess.call(["trimal -in %s_al.fasta -out %s_trimmed -gt %s" % (file, file, args.indelrep)], shell=True)
			#os.remove(file +"_al.fasta")
			os.remove(file)

	#deinterleave trimal output
	os.chdir(contigdir)

	for file in os.listdir(contigdir):
		if file.endswith('duprem_trimmed'):
			print("deinterleaving " + file)
			deinterleave_fasta(file,file+'_deint.fasta')
			os.remove(file)

	# Move raw ualigned/untrimmed locus-clusters into diploidclusters folder.
	for file in os.listdir(contigdir):
		if '_cl' in file:
			if not 'trimmed_deint' in file:
				src =contigdir + file
				print(src)
				dst = args.workingdir + 'diploidclusters/'+ file
				print(dst)
				os.rename(src, dst)
				print(file)

	os.chdir(diploidclusters)

	#Generate summary .csv's and fastas of locus clusters

	DICT2= {}

	for baitcluster in os.listdir(diploidclusters):
		if baitcluster.endswith('_'):
			if baitcluster not in DICT2:
				DICT2[baitcluster]={}

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				DICT2[baitcluster][folder]=[]


	os.chdir(diploidclusters)

	subprocess.call(["cat *_  > ALLsamples_allcontigs_allbaitclusters_contigs.fasta"], shell=True)

	#filling in the dictionary with a list of one or more contig sequences for each bait and each sample
	input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs.fasta", "fasta")

	for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys(): 
				for record in input_fasta:
					bait= record.id.split('_', 1)[0]
					baitcluster= 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
					print(baitcluster)
					folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3] + '_'
					print(folder)
					seq=record.seq
					DICT2[baitcluster][folder].append(seq)

	#write fasta output summary files by baitcluster
	for baitcluster in DICT2.keys():
		if len(DICT2[baitcluster])>0:
			outfile = open(baitcluster+"_allsamples_allcontigs.fasta", 'w+')
			for folder in DICT2[baitcluster].keys():
				if len(DICT2[baitcluster][folder])>0:
					seq_list = DICT2[baitcluster][folder] 
					sorted_seq_list = sorted(seq_list, key = lambda id: int(len(seq)), reverse=True)
					for seq in sorted_seq_list:
						index=str(sorted_seq_list.index(seq))
						#print '>'+baitcluster+folder+index+'\n'+seq+'\n'
						outfile.write(str('>'+baitcluster+folder+index+'\n'+seq+'\n'))


	#output the nested dictory to a csv file that can be exported into an excel table where rows are baitclusters, columns are samples, and cell values are number of contigs 
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	header = ['baitcluster']+columns

	##print header
	with open('ALLsamples_allclusterbaits_cid'+ args.clust2id[1:] + '_SUMMARY_TABLE.csv', 'w') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		first_value = list(DICT.values())[0]
		samples = sorted(first_value.keys())
		#samples = columns[0:]
		for baitcluster in DICT2.keys():
			writer.writerow([baitcluster]+[len(DICT2[baitcluster][sample]) for sample in samples])

	os.chdir(contigdir)

	if 'full' in args.idformat:
		#Reformat as >L100_cl0_@@##_sampleid ; Keeps full anottation.
		for file in os.listdir(contigdir):
			if file.endswith('deint.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl=line.split(' ')[0]
							linspl2=linspl.split('_')
							print(linspl2)
							name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '\n'
							print(name)
							replaceAll(file, line, name)
	else:
		if 'onlysample' in args.idformat:
			#Reformat as >@@##_sampleid ; simplest format for concatenation across locus-cluster. May have to deal with multiple allele copies per sample with identical id names.
			for file in os.listdir(contigdir):
				if file.endswith('deint.fasta'):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								linspl=line.split(' ')[0]
								linspl2=linspl.split('_')
								print(linspl2)
								name = '>' + linspl2[2] + '_' + linspl2[3] + '\n'
								print(name)
								replaceAll(file, line, name)
		else:
			sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")



# You can concatenate trimmed clustered-loci labelled deintereleaved_...._trimmed in the 'diploids' folder for maximum likelihood phylogeny (i.e. raxml, iqtree),
#resulting MSA's may have more than one sequence per sample and could be due to:
#retention of paralogous sequences, may want to set a higher -c2 identity, atleast 75% recommended.
#heterozygous variants were retained and not collapsed into consensus alleles, check clustering id for -c1, .99 recommended
#In our assessment multiple contigs is usually the result of two locus fragments that were not contiguous and unable to be joined in the consensus allele clustering step
#The purpose of this stage is to reduce multi-copy sequences due to paralogy by identity clustering among samples for the same reference locus.
#However, some datasets based on a variety of rich to poor DNA inputs, may have significantly higher multi-copy sequences at a given locus
# due to unjoined fragmented amplicons of the same locus which is common with poor sample DNA quality.
#we leave the decision up to the user as to how to manage extra contig/scaffold sequences present in locus clusters.(e.g. choose the longest sequence)
