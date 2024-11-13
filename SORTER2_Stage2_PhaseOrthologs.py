import sys
import subprocess
import fileinput
import csv
import itertools
import argparse
import re
import glob
import os
from os import path
import shutil
from shutil import move
from shutil import copyfile
import Bio
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-pq", "--phasequal")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
parser.add_argument("-idformat", "--idformat")
args = parser.parse_args()

diploidclusters = args.workingdir + 'diploidclusters/'
contigdir=args.workingdir + 'diploids/'


#Make diploids and diploid clusters folders if needed

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

folders_to_check = ["diploids_phased"]

for folder_name in folders_to_check:
	check_folder(args.workingdir, folder_name)

#Set Directory Variables
phasedir = args.workingdir + 'diploids_phased/'
direc=os.listdir(args.workingdir)
map_contigs_to_baits_dir=sorted(os.listdir(contigdir))

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


print("Removing any previous phase files...")

os.chdir(args.workingdir)

#Get rid of previous read statistics/heterozygosity summary table
for file in os.listdir(args.workingdir):
	if file.endswith('readstats_fin.csv'):
		os.remove('readstats_fin.csv')

os.chdir(phasedir)

for file in os.listdir(phasedir):
	os.remove(phasedir + file)

os.chdir(diploidclusters)

for folder in os.listdir(args.workingdir):
	if folder.endswith("iploidclusters"):
		os.chdir(args.workingdir + folder)
		dirpath =  args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith('_degap.fasta'):
				os.remove(dirpath + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("iploidclusters"):
		os.chdir(args.workingdir + folder)
		dirpath =  args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith('master.udb'):
				os.remove(dirpath + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("iploidclusters"):
		os.chdir(args.workingdir + folder)
		dirpath =  args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith('clusterannotated.fasta'):
				os.remove(dirpath + file)

os.chdir(args.workingdir)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		for file in os.listdir(args.workingdir + folder):
			if file.endswith('Final.fasta'):
				os.remove(args.workingdir + folder + '/' + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		for file in os.listdir(args.workingdir + folder):
			if file.endswith('.csi'):
				os.remove(args.workingdir + folder + '/' + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		for file in os.listdir(args.workingdir + folder):
			if file.endswith('.bam'):
				os.remove(args.workingdir + folder + '/' + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		for file in os.listdir(args.workingdir + folder):
			if file.endswith('vcf.gz'):
				os.remove(args.workingdir + folder + '/' + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		for file in os.listdir(args.workingdir + folder):
			if file.endswith('readstats.txt'):
				os.remove(args.workingdir + folder + '/' + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		for file in os.listdir(args.workingdir + folder):
			if file.startswith('het.'):
				os.remove(args.workingdir + folder + '/' + file)

#Degap locus-clusters
os.chdir(diploidclusters)
print('Degapping locus-clusters for phasing')

for file in os.listdir(diploidclusters):
	if file.endswith('_'):
		newfilename=file+ 'degap.fasta'
		print('degapping ' + file)
		with open(newfilename, "w") as o:
			for record in SeqIO.parse(file, "fasta"):
				record.seq = record.seq.replace("-", "")
				SeqIO.write(record, o, "fasta")

os.chdir(diploidclusters)

#deinterlieave locus-clusters
print('deinterleaving locus-clusters for phasing')
for file in os.listdir(diploidclusters):
	if file.endswith('degap.fasta'):
		print("deinterleaving " + file)
		deinterleave_fasta(file,file[:-11]+'deinterleaved_degap.fasta')
		os.remove(file)

os.chdir(diploidclusters)

#add new line in cluster files
for file in os.listdir(diploidclusters):
	if 'deinterleaved' in file:
		with open(file, 'a+') as cluster:
			cluster.write('\n')

os.chdir(diploidclusters)

#make ublast data
subprocess.call("cat *degap.fasta > ALLsamples_allcontigs_allbaits_clusterannotated.fasta", shell=True)
subprocess.call("usearch -makeudb_usearch ALLsamples_allcontigs_allbaits_clusterannotated.fasta -output diploid_master.udb", shell=True)

print('Preparing Phase Directories...')

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		dirpath = args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if "annotated.fasta" in file:
				os.remove(dirpath + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		dirpath =  args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if "Final.fasta" in file:
				os.remove(dirpath + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		dirpath =  args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if "chimera.fasta" in file:
				os.remove(dirpath + file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir + folder)
		dirpath =  args.workingdir + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if "mapreads.bam" in file:
				os.remove(dirpath + file)

os.chdir(diploidclusters)

# # ##Recompile sample specific clusterbaits into fastq folders for readmapping
for file in os.listdir(diploidclusters):
	if 'deinterleaved' in file:
		with open(file, 'r') as allsamplefile:
			for line in allsamplefile:
				if '>' in line:
					linspl=line.split('_')
					sample=linspl[2] + '_' + linspl[3].strip('\n') + '_allcontigs_allclusterbaits_annotated.fasta'
					print(linspl)
					with open(os.path.join(args.workingdir + linspl[2] + '_' + linspl[3].strip('\n') + '_assembly/', sample), 'a+') as idx:
						if not line.endswith('\n'):
							idx.write('\n')
						else:
							while True:
								try:
									idx.write(line)
									seq = next(allsamplefile)
									idx.write(seq)
									print(line)
									#print seq
									break
								except StopIteration as e:
									print(e)
									break

os.chdir(args.workingdir)


#Map reads to consensus allele references, phase bi-allelic haplotypes
for folder in direc:
	if 'assembly' in folder:
		os.chdir(args.workingdir + folder)
		for baits in os.listdir(args.workingdir + folder):
			if baits.endswith('allcontigs_allclusterbaits_annotated.fasta'):
				os.chdir(args.workingdir + folder)
				read = folder[:-8] + 'R1_val_1.fq'
				print(read)
				R2 = folder[:-8] + 'R2_val_2.fq'
				print(R2)
				dirpath = args.workingdir + folder + '/'
				subprocess.call(["bwa index %s" % (baits)], shell=True)
				subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (dirpath + baits, dirpath + read, dirpath + R2, folder[:-8])], shell=True)
				subprocess.call(["samtools sort  %smapreads.sam -o %s" % (dirpath+folder[:-8], dirpath+folder[:-8] + 'mapreads.bam')], shell=True)
				subprocess.call(["samtools index  %s" % (dirpath+folder[:-8] + 'mapreads.bam')], shell=True)
				subprocess.call(["samtools phase -A -Q %s -b %s %s" % (args.phasequal, dirpath+folder[:-8], dirpath+folder[:-8] + 'mapreads.bam')], shell=True)
				subprocess.call(["samtools sort  %s -o %s" % (dirpath+folder[:-8] + '.0.bam', dirpath+folder[:-8] + '0srt.bam')], shell=True)
				subprocess.call(["samtools sort  %s -o %s" % (dirpath+folder[:-8] + '.1.bam', dirpath+folder[:-8] + '1srt.bam')], shell=True)
				subprocess.call(["samtools sort  %s -o %s" % (dirpath+folder[:-8] + '.chimera.bam', dirpath+folder[:-8] + 'chimerasrt.bam')], shell=True)
				subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s | bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s " % (dirpath+baits, dirpath+folder[:-8] + '0srt.bam', dirpath+folder[:-8] + '0.vcf.gz' )], shell=True)
				subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s | bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s " % (dirpath+baits, dirpath+folder[:-8] + '1srt.bam', dirpath+folder[:-8] + '1.vcf.gz' )], shell=True)
				subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s | bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s " % (dirpath+baits, dirpath+folder[:-8] + 'chimerasrt.bam', dirpath+folder[:-8] + 'chimera.vcf.gz' )], shell=True)
				subprocess.call(["bcftools index  %s" % (dirpath+folder[:-8] + '0.vcf.gz')], shell=True)
				subprocess.call(["bcftools index  %s" % (dirpath+folder[:-8] + '1.vcf.gz')], shell=True)
				subprocess.call(["bcftools index  %s" % (dirpath+folder[:-8] + 'chimera.vcf.gz')], shell=True)
				subprocess.call(["cat %s | bcftools consensus %s > %s0_Final.fasta" % (dirpath+baits, dirpath+folder[:-8] + '0.vcf.gz', folder[:-8])], shell=True)
				subprocess.call(["cat %s | bcftools consensus %s > %s1_Final.fasta" % (dirpath+baits, dirpath+folder[:-8] + '1.vcf.gz', folder[:-8])], shell=True)
				#subprocess.call(["cat %s | bcftools consensus %s > %sChimera.fasta" % (baits, folder[:-8] + 'chimera.vcf.gz', folder[:-8])], shell=True)
				os.remove(folder[:-8] + 'mapreads.sam')
				os.remove(folder[:-8] + '.0.bam')
				os.remove(folder[:-8] + '0srt.bam')
				os.remove(folder[:-8] + '0.vcf.gz')
				os.remove(folder[:-8] + '0.vcf.gz.csi')
				os.remove(folder[:-8] + '.1.bam')
				os.remove(folder[:-8] + '1srt.bam')
				os.remove(folder[:-8] + '1.vcf.gz')
				os.remove(folder[:-8] + '1.vcf.gz.csi')
				os.remove(folder[:-8] + '.chimera.bam')
				os.remove(folder[:-8] + 'chimerasrt.bam')
				#os.remove(folder[:-8] + 'chimera.vcf.gz')
				#generate read statistics text file
				bam = dirpath+folder[:-8] + 'mapreads.bam'
				statfilename = folder[:-8] + "readstats.txt"
				with open(os.path.join(args.workingdir + folder, statfilename), 'a+') as statfile:
							statfile.write( folder[:-8] + " Read Statistics" + '\n')
							subprocess.call(["samtools flagstat %s >> %s" % (bam, statfilename)], shell=True)
							#get read depth
							subprocess.call(["samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (bam, statfilename)], shell=True)
							#get Coverage
							subprocess.call(["samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (bam, statfilename)], shell=True)
							statfile.close()


os.chdir(args.workingdir)

# annotate phased fasta files with alternative phase
for folder in direc:
	if 'assembly' in folder:
		os.chdir(args.workingdir + folder)
		readir = os.listdir(args.workingdir + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("0_Final.fasta"):
				with open(file, 'r') as infile:
 					for line in infile:
 						if '>' in line:
 							print(line)
 							name = line.rstrip('\n') + '_ph0' + '\n'
 							print(name)
 							replaceAll(file, line, name)
os.chdir(args.workingdir)

for folder in direc:
	if 'assembly' in folder:
		os.chdir(args.workingdir + folder)
		readir = os.listdir(args.workingdir + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("1_Final.fasta"):
				with open(file, 'r') as infile:
 					for line in infile:
 						if '>' in line:
 							print(line)
 							name = line.rstrip('\n') + '_ph1' + '\n'
 							print(name)
 							replaceAll(file, line, name)

os.chdir(args.workingdir)

##Concatenate Phased seqs files
for folder in direc:
	if 'assembly' in folder:
		os.chdir(args.workingdir + folder)
		subprocess.call(["cat *_Final.fasta > %sallcontigs_allclusterbaits_contigs_phased.fasta" % (folder[:-8])], shell=True)

#Move Phased allcontigs_allbaits files to args.contigdir 'diploids_phased/' directory
for folder in direc:
	if 'assembly' in folder:
		os.chdir(args.workingdir + folder)
		readir = os.listdir(args.workingdir + folder)
		for file in readir:
			if file.endswith('_allcontigs_allclusterbaits_contigs_phased.fasta'):
				src = args.workingdir + folder + '/' + file
				dst = phasedir + file
				os.rename(src, dst)

os.chdir(diploidclusters)

DICT2= {}

for baitcluster in os.listdir(diploidclusters):
	if baitcluster.endswith('_'):
		if baitcluster not in DICT2:
			DICT2[baitcluster]={}

#print(DICT2)

for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
 			for baitcluster in DICT2.keys(): 
 				DICT2[baitcluster][folder + 'ph0']=[]

for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
 			for baitcluster in DICT2.keys(): 
 				DICT2[baitcluster][folder + 'ph1']=[]

#print(DICT2)

os.chdir(phasedir)

#concatenated phased sequences for all samples

subprocess.call(["cat *_allcontigs_allclusterbaits_contigs_phased.fasta  > ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta"], shell=True)

#>L4_cl1_WB05_glycyrrhiza_1_ph0
#filling in the dictionary with a list of one or more contig sequences for each bait and each sample
input_fasta=SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta", "fasta")
for folder in map_contigs_to_baits_dir:
	if folder.endswith('_'):
		for baitcluster in DICT2.keys(): 
			for record in input_fasta:
				bait= record.id.split('_', 1)[0]
				baitcluster= 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
				print(baitcluster)
				phase = (record.id.split('_', 4)[4]).rstrip('\n')
				print(phase)
				#folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3]+ '_' + phase
				folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3]+ '_' + phase
				print(record.id.split('_', 4)[3])
				print(folder)
				seq=record.seq
				DICT2[baitcluster][folder].append(seq)
				#print(DICT2)
				#print baitcluster, folder, len(DICT2[baitcluster][folder])


#write fasta output summary files by bait
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
					outfile.write(str('>'+baitcluster+folder+'_'+index+'\n'+seq+'\n'))

# ##output the nested dictory to a csv file that can be exported into an excel table where rows are baits, columns are samples, and cell values are number of contigs 
columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]

ph0 = ["{}{}".format(i,'ph0') for i in columns]
ph1 = ["{}{}".format(i,'ph1') for i in columns]

header = ['baitcluster']+ph0+ph1

with open('ALLsamples_allcontigs_allbaits_SUMMARY_TABLE_phased.csv', 'w') as outfile:
	writer = csv.writer(outfile)
	writer.writerow(header)
	first_value = list(DICT2.values())[0]
	samples = sorted(first_value.keys())
	#samples = columns[0:]
	for baitcluster in DICT2.keys():
		writer.writerow([baitcluster]+[len(DICT2[baitcluster][sample]) for sample in samples])

os.chdir(phasedir)

with open('ALLsamples_allcontigs_allbaits_SUMMARY_TABLE_phased.csv', 'w') as outfile:
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	ph0 = ["{}{}".format(i,'ph0') for i in columns]
	ph1 = ["{}{}".format(i,'ph1') for i in columns]
	header = ['baitcluster']+ph0+ph1
	writer = csv.writer(outfile)
	writer.writerow(header)

#Align and trim clusterbaits. 
#Regions which don't share overlap (i.e. region with unique indel) in atleast x percent of the samples are removed. All clusters with single sequences not aligned and compiled downstream

print('Only keeping the longest sequence for sequences with duplicate name IDs prior to alignment/trimming')
for file in os.listdir(phasedir):
	if file.endswith('allsamples_allcontigs.fasta'):
		with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							linspl=line.split('_')
							name = linspl[0] + '_' + linspl[1] + '_' + linspl[2] + '_' + linspl[3] + '_' + linspl[4] + '\n'
							replaceAll(file, line, name)

for file in os.listdir(phasedir):
	if file.endswith('allsamples_allcontigs.fasta'):
		remove_dup(file, file[:-6] +'_duprem.fasta')

for file in os.listdir(phasedir):
	if file.endswith('_duprem.fasta'):
		subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file[:-6])], shell=True)
		subprocess.call(["trimal -in %s -out %s -gt %s" % (file[:-6] + "_al.fasta", file[:-6] + '_trimmed.fasta', args.indelrep)], shell=True)
		#os.remove(file[:-6]+'_al.fasta')

#deinterleave
for file in os.listdir(phasedir):
	if 'trimmed' in file:
		print("deinterleaving " + file)
		deinterleave_fasta(file,file[:-6]+'_final.fasta')
		os.remove(file)

#compile stats
os.chdir(args.workingdir)

#Make dictionary of individuals for readstats and heterozygosity
HETDICT= {}

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		ind=folder[:-9]
		print(ind)
		if ind not in HETDICT:
			HETDICT[ind]={}

for samp in os.listdir(args.workingdir):
	if samp.endswith('assembly'):
		os.chdir(args.workingdir+samp)
		for readstat in os.listdir(args.workingdir+samp):
			if readstat.endswith('readstats.txt'):
				for ind in HETDICT:
					if ind in readstat:
						print(readstat)
						with open(readstat, "r") as statfile:
							lines_to_read = [16, 17]
							for position, line in enumerate(statfile):
								if '+' in line:
									statlabela = line.split(" ")[3]
									statlabel = statlabela.strip('\n')
									statinta = line.split(" ")[0]
									statint = statinta.strip('\n')
									print(statlabel)
									print(statint)
									HETDICT[ind][statlabel]=[]
									HETDICT[ind][statlabel].append(int(statint))
								else:
									if position in lines_to_read:
										if position == 16:
											statlabel = 'readdepth'
											statint = line.strip('\n')
											print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel]=[]
											HETDICT[ind][statlabel].append(int(float(statint)))
										elif position == 17:
											statlabel = 'coverage'
											statint = line.strip('\n')
											print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel]=[]
											HETDICT[ind][statlabel].append(int(float(statint)))
os.chdir(args.workingdir)

#get readstat values
het_values_list = list(HETDICT.values())
columns = [x for x in het_values_list[0]]
header = ['Individual']+columns

with open('readstats.csv', 'w') as outfile:
	writer = csv.writer(outfile)
	writer.writerow(header)
	statlist = list(HETDICT.values())
	stats = statlist[0].keys()
	#samples = columns[0:]
	for ind in HETDICT.keys():
		writer.writerow([ind]+[HETDICT[ind][stat] for stat in stats])

with open('readstats.csv', 'r') as infile:
    with open('readstats_fin.csv', 'w') as outfile:
	    data = infile.read()
	    data = data.replace("]", "")
	    data = data.replace("[", "")
	    outfile.write(data)

for file in os.listdir(args.workingdir):
	if file.endswith("readstats.csv"):
		os.remove(file)

os.chdir(phasedir)

if 'full' in args.idformat:
	#Reformat as >L100_cl0_@@##_sampleid_phase_index ; Keeps full anottation.
	for file in os.listdir(phasedir):
		if file.endswith('final.fasta'):
			with open(file, 'r') as infile:
				for line in infile:
					if '>' in line:
						print(line)
						linspl=line.split(' ')[0]
						linspl2=linspl.split('_')
						print(linspl2)
						name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4]
						print(name)
						replaceAll(file, line, name)
	sys.exit('Kept full sequence ID annotations; e.g. >L100_cl0_WA10_sampleid_0')
else:
	if 'phase' in args.idformat:
		for file in os.listdir(phasedir):
			if file.endswith('final.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl=line.split(' ')[0]
							linspl2=linspl.split('_')
							print(linspl2)
							name = '>' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4]
							print(name)
							replaceAll(file, line, name)
		sys.exit('Annotated alignments as: >@@##_sampleid_0/1 (annotated with phase)')
	else:
		if 'onlysample' in args.idformat:
			#Reformat as >@@##_sampleid ; simplest format for concatenation across locus-cluster. Will have to decide how to manage alleles/sequence copies per sample with identical id names.
				if file.endswith('final.fasta'):
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
				sys.exit('Annotated alignments as: >@@##_sampleid (no phase annotations)')
		else:
			sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")
