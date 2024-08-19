import os
from os import path
import sys
import subprocess
import shutil
from shutil import move
import fileinput
import csv
import sys
import subprocess
import Bio
from os import path
from shutil import move
from shutil import copyfile
from Bio import SeqIO
import pandas as pd
import re
import itertools
import argparse
import glob

parser = argparse.ArgumentParser()

parser.add_argument("-wd", "--workingdir")
parser.add_argument("-loci", "--locinum")
parser.add_argument("-csn", "--contigscafnum")
parser.add_argument("-csl", "--contigscaflen")
parser.add_argument("-ref", "--refdir")
parser.add_argument("-pq", "--phasequal")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
parser.add_argument("-fp", "--filterundiff")

args = parser.parse_args()
phaseset=args.workingdir + 'phaseset/'
baitid1= ["L%d_" % x for x in range(int(args.locinum))]
baitid= ["L%d" % x for x in range(int(args.locinum))]
diploidclusters=args.workingdir + 'diploids_phased/'
diploid_db = args.workingdir + 'diploids_phased/diploid_master.udb'


# #Define command to change sequence IDs
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

#Define function to annotate duplicate sequence IDS

def annotate_fasta(input_file):
    # Calculate the output file name
    output_file = input_file[:-3] + "_annotated.fasta"
    
    # Dictionary to store the suffix counts for each sequence name
    names = {}

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read the input file
        content = infile.read()
        
        # Split the file content into individual sequences
        sequences = content.split('>')[1:]  # Skip the first empty split before the first '>'
        
        for seq in sequences:
            lines = seq.strip().split('\n')
            name = lines[0]
            sequence = ''.join(lines[1:])
            
            # Determine the suffix for the sequence name
            if name in names:
                names[name] += 1
                suffix = f"_{names[name]}"
            else:
                names[name] = 1
                suffix = "_1"
            
            # Write the annotated sequence to the output file
            outfile.write(f">{name}{suffix}\n{sequence}\n")


print("Removing any previous phase files...")

os.chdir(phaseset)

for folder in os.listdir(phaseset):
	if 'clusters_phased' in folder:
		shutil.rmtree(phaseset + folder)

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if not 'trimming_report' in file:
				if not 'spades' in file:
					if not '_val_' in file:
						os.remove(phaseset + folder + '/'+ file)
						os.chdir(phaseset)

os.chdir(diploidclusters)

#make ublast data
subprocess.call("usearch -makeudb_usearch ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta -output diploid_master.udb", shell=True)

#Map Contigs to Rereferences
os.chdir(phaseset)

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		read = folder[:-9] + 'R1_val_1.fq'
		R2 = read[:-9] + 'R2_val_2.fq'
		for file in iterpath:
			if file.endswith("_hybrid_assembly"):
				os.chdir(phaseset + folder + "/" + file + "/")
				filepath =  phaseset + folder + "/" + file + "/"
				iterpath2 = os.listdir(filepath)
				for contig in iterpath2:
					if contig.endswith("contigs.fasta"):
						conscontig= args.refdir
						contigmap=phaseset+ folder + "/" + file + "/"+ contig
						os.chdir(dirpath)
						subprocess.call(["bwa mem -V %s %s > %s_contigmap.sam" % (conscontig, contigmap, dirpath+folder)], shell=True)
						subprocess.call(["samtools view -S -F 4 %s_contigmap.sam | awk -v OFS='\t' '{print \">\" $3\"_\" \"\\n \" $10}' > %scontigmap.fa " % (dirpath+folder,folder + '_')], shell=True)
						src = folder +'_contigmap.fa'
						dst = dirpath + folder + '_contigmap.fa'
						os.rename(src, dst)
						os.remove(folder + '_contigmap.sam')
os.chdir(phaseset)

#Make fasta files for contigs which mapped to the same reference
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		print(folder)
		os.chdir(phaseset + folder)
		readir = os.listdir(phaseset + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("contigmap.fa"):
				print(file)
				with open(file, 'r') as contigfile:
					for line in contigfile:
						for id in baitid1:
							if id in line:
								#print(line)
								with open(os.path.join(phaseset + folder, id), 'a') as idx:
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

# #Rename contigs with locus ID
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		readir = os.listdir(phaseset + folder)
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
os.chdir(phaseset)

print('Take ' + args.contigscafnum + ' longest contigs for each locus per sample, then removing any contigs s smaller than ' + args.contigscaflen + ' bp')

#Take longest contig from contig set
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		readir = os.listdir(phaseset + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("_"):
				#Get longest Contigs
				extract_longest_sequences(file, args.contigscafnum)
				#Remove sequences shorter than x bp
				subprocess.call(["seqtk seq -L %s %slongest.fa > %slongestfiltered.fa" % (args.contigscaflen, file, file)], shell=True )

os.chdir(phaseset)

#cluster highly similar contigs (i.e collapsing heterozygotes, or possibly homeologues depending on genetic distance)
print('Clustering congtigs into Consensus Homeologs at 99% Identity Threshold')

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		print(folder)
		os.chdir(phaseset + folder)
		readir = os.listdir(phaseset + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("longestfiltered.fa"):
				subprocess.call(["usearch -cluster_fast %s -id 0.99 -consout %s_cons.fa" % (file, file[:-3])], shell=True)


os.chdir(phaseset)

#annotate dupe sequence ID's
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		readir = os.listdir(phaseset + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("_cons.fa"):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							name = '>' + folder[:-9] + '_' + file[:-24] + '\n'
							print(name)
							replaceAll(file, line, name)	

os.chdir(phaseset)

#annotate cluster contigs with sample ID and locus
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		readir = os.listdir(phaseset + folder)
		subprocess.call(["pwd"], shell=True)
		for file in readir:
			if file.endswith("_cons.fa"):
				os.chdir(phaseset + folder)
				annotate_fasta(file)
				os.remove(file)

os.chdir(phaseset)

#remove unused sequences
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith("_"):
				print('deleting: ' + file)
				os.remove(dirpath + file)

os.chdir(phaseset)

#remove unused sequences
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith("_longest.fa"):
				print('deleting: ' + file)
				os.remove(dirpath + file)

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		dirpath =  phaseset + folder + "/"
		iterpath = os.listdir(dirpath)
		for file in iterpath:
			if file.endswith("_longestfiltered.fa"):
				print('deleting: ' + file)
				os.remove(dirpath + file)

os.chdir(phaseset)

#deinterleave and annotate dupes
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		subprocess.call(["cat *_annotated.fasta  > %s_allcontigs_allbaits.fasta" % (folder[:-9])], shell=True)
		print("deinterleaving " + folder[:-9] +'_allcontigs_allbaits.fasta')
		deinterleave_fasta(folder[:-9] +'_allcontigs_allbaits.fasta', folder[:-9]+'_allcontigs_allbaits_deinterleaved.fasta')
		os.remove(folder[:-9] +'_allcontigs_allbaits.fasta')

os.chdir(phaseset)

#strip lines for final contigfile
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		for file in os.listdir(phaseset + folder):
			if file.endswith('deinterleaved.fasta'):
				with open(file, 'r') as infile:
					with open(file[:-6] + '_final.fasta', 'w') as outfile:
						for line in infile:
							if line.strip():
								outfile.write(line)


os.chdir(phaseset)

#Map reads to locus-clustered contig-cluster references, get consensus fasta of phased sequences
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset + folder)
		for baits in os.listdir(phaseset + folder):
			if baits.endswith('final.fasta'):
				read = folder[:-8] + 'R1_val_1.fq'
				print(read)
				R2 = folder[:-8] + 'R2_val_2.fq'
				print(R2)
				print(baits)
				dirpath = phaseset + folder + '/'
				subprocess.call(["bwa index %s" % (baits)], shell=True)
				subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (dirpath+baits, dirpath+read, dirpath+R2, folder[:-8])], shell=True)
				subprocess.call(["samtools sort %smapreads.sam -o %s" % (dirpath+folder[:-8], dirpath+folder[:-8] + 'mapreads.bam')], shell=True)
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
				subprocess.call(["cat %s | bcftools consensus %s > %s0.fasta" % (dirpath+baits, dirpath+folder[:-8] + '0.vcf.gz', folder[:-8])], shell=True)
				subprocess.call(["cat %s | bcftools consensus %s > %s1.fasta" % (dirpath+baits, dirpath+folder[:-8] + '1.vcf.gz', folder[:-8])], shell=True)
				deinterleave_fasta(dirpath+folder[:-8] +'0.fasta', dirpath+folder[:-9]+'_0_Final.fasta')
				deinterleave_fasta(dirpath+folder[:-8] +'1.fasta', dirpath+folder[:-9]+'_1_Final.fasta')
				os.remove(folder[:-8] + '0.fasta')
				os.remove(folder[:-8] + '1.fasta')
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

os.chdir(phaseset)

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset+folder)
		dirpath = phaseset+folder+'/'
		for bam in os.listdir(phaseset + folder):
			if bam.endswith("mapreads.bam"):
				bam = dirpath+folder[:-8] + 'mapreads.bam'
				statfilename = folder[:-8] + "readstats.txt"
				with open(os.path.join(phaseset + folder, statfilename), 'a+') as statfile:
					statfile.write( folder[:-8] + " Read Statistics" + '\n')
					subprocess.call(["samtools flagstat %s >> %s" % (bam, statfilename)], shell=True)
					#get read depth
					subprocess.call(["samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (bam, statfilename)], shell=True)
					#get Coverage
					subprocess.call(["samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (bam, statfilename)], shell=True)
					statfile.close()


os.chdir(args.workingdir)
#Make a copy of diploid locus-clusters for each sample in phase set
dpdst= phaseset + 'diploidclusters_phased/'
shutil.copytree(diploidclusters, dpdst)

for folder in os.listdir(phaseset):
	if 'diploidclusters' in folder:
		dipcl=phaseset+folder
		os.chdir(phaseset+folder)
		for file in os.listdir(dipcl):
			if 'duprem' in file:
				os.remove(file)
			else:
				if 'phased' in file:
					os.remove(file)
				else:
					if 'master.udb' in file:
						os.remove(file)


os.chdir(phaseset)

#UBLAST Phased sequences for each sample to diploid locus-clusters to determine orthology
for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset+folder)
		path=phaseset+folder
		for file in os.listdir(path):
			if file.endswith('_Final.fasta'):
				os.chdir(phaseset+folder)
				subprocess.call(["awk 'BEGIN{FS=\" \"}{if(!/>/){print toupper($0)}else{print $1}}' %s > %s_cap.fasta" % (file, file[:-6])], shell=True)
				subprocess.call(["usearch -usearch_global %s -db %s -id 0.9 -top_hit_only -blast6out %s_hits.txt -strand plus" % (file[:-6]+'_cap.fasta', diploid_db, file[:-12])], shell=True)


#Compile Polyploid into Diploid locus-cluster dataset, respectively
os.chdir(phaseset)

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset+folder)
		for file in os.listdir(phaseset+folder):
			if file.endswith('0_hits.txt'):
				with open(file, 'r') as hits:
					for hitlines in hits:
						splithits=hitlines.split('	')
						splithits2=splithits[1].split('_')
								#print splithits
								# print splithits2
								# print splithits[0]
						with open(folder[:-8] + '0_Final_cap.fasta', 'r') as phfinal:
							for line in phfinal:
								if splithits[0] in line:
									print('Match ' +splithits[0] + ' in ' + splithits[1])
									with open(phaseset + 'diploidclusters_phased/'+ splithits2[0] + '_' + splithits2[1] +'__allsamples_allcontigs.fasta', 'a+') as baitcluster:
										while True:
											try:
												splitline= line.split('_')
												#print splitline
												name = '>' + splithits2[0] + '_' + splithits2[1] +'_' +splitline[0].split('>')[1] + '_'  + splitline[1] +  '_' + splithits2[3]#see if any statement for broader annotation
												print(name.rstrip('\n') + '_ph0' + '\n')
												baitcluster.write(name.rstrip('\n') + '_ph0' + '\n')
												#seq = phfinal.next()
												seq = next(phfinal)
												#print seq
												baitcluster.write(seq)
												break
											except StopIteration as e:
												print(e)
												break

os.chdir(phaseset)

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset+folder)
		for file in os.listdir(phaseset+folder):
			if file.endswith('1_hits.txt'):
				with open(file, 'r') as hits:
					for hitlines in hits:
						splithits=hitlines.split('	')
						splithits2=splithits[1].split('_')
								#print splithits
								# print splithits2
								# print splithits[0]
						with open(folder[:-8] + '1_Final_cap.fasta', 'r') as phfinal:
							for line in phfinal:
								if splithits[0] in line:
									print('Match splithits[0]= ' +splithits[0] + ' in ' + splithits[1])
									with open(phaseset + 'diploidclusters_phased/'+ splithits2[0] + '_' + splithits2[1] +'__allsamples_allcontigs.fasta', 'a+') as baitcluster:
										while True:
											try:
												splitline= line.split('_')
												#print splitline
												name = '>' + splithits2[0] + '_' + splithits2[1] +'_' +splitline[0].split('>')[1] + '_'  + splitline[1] +  '_' + splithits2[3]#see if any statement for broader annotation
												print(name.rstrip('\n') + '_ph1' + '\n')
												baitcluster.write(name.rstrip('\n') + '_ph1' + '\n')
												#seq = phfinal.next()
												seq = next(phfinal)
												#print seq
												baitcluster.write(seq)
												break
											except StopIteration as e:
												print(e)
												break

os.chdir(phaseset)

#Rewrite polyploid names as >WA01_polyploid_diploidhit_phase

for folder in os.listdir(phaseset):
	if 'assembly' in folder:
		os.chdir(phaseset+'diploidclusters_phased/')
		sample=folder.split('_')[0] +'_' + folder.split('_')[1]
		for file in os.listdir(phaseset+'diploidclusters_phased/'):
			if file.endswith('allsamples_allcontigs.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if sample in line:
							print(line)
							linspl=line.split('_')
							print(linspl)
							name = '>' + linspl[2] + '_' + linspl[3] + '_' + linspl[4] + '_' + linspl[5]
							print(name)
							replaceAll(file, line, name)


#Rename Diploid species as >WA01_diploid_phase from >L102_cl3_WF07_diploid_ph0_0

os.chdir(phaseset+'diploidclusters_phased/')

for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('allsamples_allcontigs.fasta'):
		os.chdir(phaseset+'diploidclusters_phased/')
		with open(file, 'r') as infile:
			for line in infile:
				for folder in os.listdir(args.workingdir):
					if 'assembly' in folder:
						sample=folder.split('_')[0] +'_' + folder.split('_')[1]
						if sample in line:
							print(line)
							linspl=line.split('_')
							print(linspl)
							name = '>' + linspl[2] + '_' + linspl[3] + '_' + linspl[4]
							print(name)
							replaceAll(file, line, name)

os.chdir(phaseset+'diploidclusters_phased/')

# #Keep longest seq if two identical sequence ids are present 
#(i.e. diploid/polyploid samples with multiple consensus alleles/cluster; allopolyploids that did not differentiat homeologs)

print('Keeping longest sequence if two identical sequence ids are present ''\n')

for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('allsamples_allcontigs.fasta'):
		remove_dup(file, file[:-6] +'_duprem.fasta')

os.chdir(phaseset+'diploidclusters_phased/')

#Align and trim locus-clusters
for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('_duprem.fasta'):
		#remove_dup(file, file[:-6] +'_duprem.fasta')
		subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file[:-6])], shell=True)
		subprocess.call(["trimal -in %s -out %s -gt %s" % (file[:-6] + "_al.fasta", file[:-6] + '_trimmed.fasta', args.indelrep)], shell=True)
		os.remove(file[:-6] + '_al.fasta')
		os.remove(file)

os.chdir(phaseset+'diploidclusters_phased/')

#reannotate trimmed alignment to remove trimal bp annotation
for file in os.listdir(phaseset+'diploidclusters_phased/'):
	if file.endswith('trimmed.fasta'):
		with open(file, 'r') as infile:
			for line in infile:
				if '>' in line:
					print(line)
					linspl=line.split(' ')
					name = linspl[0] +'\n'
					print(name)
					replaceAll(file, line, name)

os.chdir(phaseset)

#compile read mapping statistics
#Make dictionary of individuals for readstats
HETDICT= {}

print('Getting read statistics')

for folder in os.listdir(phaseset):
	if folder.endswith("assembly"):
		ind=folder[:-9]
		print(ind)
		if ind not in HETDICT:
			HETDICT[ind]={}

for samp in os.listdir(phaseset):
	if samp.endswith('assembly'):
		os.chdir(phaseset+samp)
		for readstat in os.listdir(phaseset+samp):
			if readstat.endswith('readstats.txt'):
				for ind in HETDICT:
					if ind in readstat:
						print(readstat)
						with open(readstat, "r") as statfile:
							lines_to_read = [13, 14]
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
										if position == 13:
											statlabel = 'readdepth'
											statint = line.strip('\n')
											print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel]=[]
											HETDICT[ind][statlabel].append(int(float(statint)))
										elif position == 14:
											statlabel = 'coverage'
											statint = line.strip('\n')
											print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel]=[]
											HETDICT[ind][statlabel].append(int(float(statint)))

os.chdir(phaseset)

#get readstat values
het_values_list = list(HETDICT.values())
columns = [x for x in het_values_list[0]]
header = ['Individual']+columns

with open('phaseset_readstats.csv', 'w') as outfile:
	writer = csv.writer(outfile)
	writer.writerow(header)
	statlist = list(HETDICT.values())
	stats = statlist[0].keys()
	#samples = columns[0:]
	for ind in HETDICT.keys():
		writer.writerow([ind]+[HETDICT[ind][stat] for stat in stats])

with open('phaseset_readstats.csv', 'r') as infile:
	with open('phaseset_readstats_fin.csv', 'w') as outfile:
		data = infile.read()
		data = data.replace("]", "")
		data = data.replace("[", "")
		outfile.write(data)

for file in os.listdir(phaseset):
	if file.endswith("phaseset_readstats.csv"):
		os.remove(file)

os.chdir(phaseset+'diploidclusters_phased/')

#Generate summary table of hybrid progenitor distributions across loci

def analyze_fasta_files(fasta_files, sample_names, output_csv):
    try:
        # Extract samples identifiers from sample names
        sample_ids = { '_'.join(name.split('_')[:2]): name for name in sample_names }
        print("Sample IDs extracted: ", sample_ids)
        #make dictionary to hold the unique count data
        data = {os.path.basename(fasta): {sample: set() for sample in sample_ids.keys()} for fasta in fasta_files}
        print("Initialized data structure: ", data)
        
        for fasta_file in fasta_files:
            fasta_name = os.path.basename(fasta_file)
            print(f"Processing file: {fasta_file}")
            # Read the FASTA file
            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = record.id
                # Extract the first two underscore-delimited strings
                sample_id = '_'.join(seq_id.split('_')[:2])
                
                if sample_id in sample_ids:
                    # Extract the third underscore-delimited string
                    unique_str = seq_id.split('_')[2]
                    
                    if fasta_name not in data:
                        data[fasta_name] = {}
                    
                    if sample_id not in data[fasta_name]:
                        data[fasta_name][sample_id] = set()
                    
                    data[fasta_name][sample_id].add(unique_str)
        
        # Convert sets to counts
        for fasta_name in data:
            for sample in data[fasta_name]:
                data[fasta_name][sample] = len(data[fasta_name][sample])
        
        # Create a DataFrame from the data dictionary
        df = pd.DataFrame(data).T
        df.index.name = 'alignment'
        df.reset_index(inplace=True)
        
        # Write the DataFrame to a CSV file
        df.to_csv(output_csv, index=False)
    
    except Exception as e:
        print(f"An error occurred: {e}")

os.chdir(phaseset+'diploidclusters_phased/')

#Generate summary table of progenitor representation across all loci
# Get FASTA file names
fasta_files = [file for file in os.listdir(phaseset+'diploidclusters_phased/') if file.endswith("trimmed.fasta")]
print(fasta_files)
sample_names = [file for file in os.listdir(phaseset) if file.endswith("assembly")]
print(sample_names)
output_csv = "progenitor_distributions.csv"

os.chdir(phaseset+'diploidclusters_phased/')

analyze_fasta_files(fasta_files, sample_names, output_csv)

def count_unique_strings_across_fasta(fasta_files, sample_names, output_csv):
    try:
        # Extract sample identifiers from sample names
        sample_ids = { '_'.join(name.split('_')[:2]): name for name in sample_names }
        print("Sample IDs extracted: ", sample_ids)
        
        # make a dictionary to hold the unique count data
        data = {sample: [] for sample in sample_ids.keys()}
        
        for fasta_file in fasta_files:
            fasta_name = os.path.basename(fasta_file)
            print(f"Processing file: {fasta_file}")
            # Temporary dictionary for this FASTA file
            temp_data = {sample: set() for sample in sample_ids.keys()}
            
            # Read the FASTA file
            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = record.id
                # Extract the first two underscore-delimited strings
                sample_id = '_'.join(seq_id.split('_')[:2])
                
                if sample_id in sample_ids:
                    # Extract the third underscore-delimited string
                    unique_str = seq_id.split('_')[2]
                    temp_data[sample_id].add(unique_str)
            
            # Convert sets to counts and store in main data dictionary
            for sample in temp_data:
                data[sample].append(len(temp_data[sample]))
        
        # Determine the maximum number of unique strings for column headers
        max_unique_count = max(max(counts) for counts in data.values())
        
        # Initialize a dictionary for the final counts
        final_counts = {sample: {f'{i}_progenitors': 0 for i in range(max_unique_count + 1)} for sample in sample_ids.keys()}
        
        # Populate the final counts dictionary
        for sample, counts in data.items():
            for count in counts:
                final_counts[sample][f'{count}_progenitors'] += 1
        
        # Create a DataFrame from the final counts dictionary
        df = pd.DataFrame(final_counts).T
        df.index.name = 'samples'
        df.reset_index(inplace=True)
        
        # Write the DataFrame to a CSV file
        df.to_csv(output_csv, index=False)
    
    except Exception as e:
        print(f"An error occurred: {e}")

summary_csv = "progenitor_summary.csv"

count_unique_strings_across_fasta(fasta_files, sample_names, summary_csv)

#Filter samples in alignments without undifferentiated sequences
if args.filterundiff == 'T':

	os.chdir(phaseset+'diploidclusters_phased/')

	#Define function to filter sequences that don't have atleast 2 unique progenitor names represented in locus sequences

	def filter_sequences(input_file):
		# Read the input FASTA file
		records = list(SeqIO.parse(input_file, "fasta"))

		# Group sequences based on the first two tab-delimited items
		grouped_sequences = {}
		for record in records:
			header_items = record.description.split("_")
			sample_id = tuple(header_items[:2])

			if len(header_items) == 4:
				# Get sequences with 4 tab-delimited items
				if sample_id not in grouped_sequences:
					grouped_sequences[sample_id] = set()

				grouped_sequences[sample_id].add(header_items[2])

		# Identify hybrid to be removed that only has one progenitor represented at the locus
		samples_to_remove = {sample_id for sample_id, unique_strings in grouped_sequences.items() if len(unique_strings) < 2}

		# Filter sequences
		filtered_records = [record for record in records if len(record.description.split("_")) == 3
							or tuple(record.description.split("_")[:2]) not in samples_to_remove]

		# Write the filtered sequences to a new FASTA file
		output_file = input_file.replace(".fasta", "_differentiated.fasta")
		SeqIO.write(filtered_records, output_file, "fasta")


	for file in os.listdir(phaseset+'diploidclusters_phased/'):
		if file.endswith('_trimmed.fasta'):
			filter_sequences(file)

