from os import path
from shutil import move
import fileinput
import os
import csv 
import sys
import subprocess
import shutil
import Bio
from os import path
from shutil import move
from shutil import copyfile
from Bio import SeqIO
import re
import itertools
import argparse
import glob
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-rep", "--repfilt")
parser.add_argument("-mc", "--majorclusters")
parser.add_argument("-al", "--keepal")
parser.add_argument("-st2", "--stage2")
parser.add_argument("-st3", "--stage3")
parser.add_argument("-dovcf", "--dovcf")
args = parser.parse_args()


diploids = args.workingdir + 'diploids/'
diploidclusters= args.workingdir + 'diploidclusters/'


#get total number of samples in data
os.chdir(diploids)
sample_files = [f for f in os.listdir(diploids) if f.endswith('_')]
samplenum = len(sample_files)

#Define function to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
	for line in fileinput.input(file, inplace=1):
		if searchExp in line:
			line = line.replace(searchExp,replaceExp)
		sys.stdout.write(line)

#Function to count fasta sequences
def count_sequences_in_fasta(fasta_file):
    """Counts the number of sequences in a FASTA file."""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

#Filter based on sample representation
def filter_fasta(fasta_dir, threshold):
    # Calculate the required number of sequences for the threshold
    
    fl = float(threshold)
    required_count = float(samplenum * fl)
    print(f"Required count: {required_count}")

    # Create the output directory
    output_dir = os.path.join(fasta_dir, f"{int(fl * 100)}rep")
    os.makedirs(output_dir, exist_ok=True)

    # Process each FASTA file in the directory
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith("_al.fasta"):
            fasta_path = os.path.join(fasta_dir, fasta_file)
            
            # Count the sequences in the FASTA file
            seq_count = count_sequences_in_fasta(fasta_path)
            print(f"File: {fasta_file}, Sequence count: {seq_count}, Required seqs: {required_count}")

            # Check if the number of sequences meets the threshold
            if int(seq_count) >= int(required_count):
                shutil.copy(fasta_path, output_dir)

os.chdir(diploidclusters)

# filter based on sample representation
filter_fasta(diploidclusters, args.repfilt)

filtfl= float(args.repfilt)
filtstr=f"{int(filtfl * 100)}rep"
filtdir=diploidclusters+filtstr

os.chdir(filtdir)

if args.keepal == 'T':

	os.chdir(filtdir)

	print("Renaming filtered alignment sequence IDs to >ID_species")
	for file in os.listdir(filtdir):
		if file.endswith('al.fasta'):
			with open(file, 'r') as infile:
				for line in infile:
					if '>' in line:
						linspl2=line.split('_')
						name = '>' + linspl2[2] + '_' + linspl2[3]
						replaceAll(file, line, name)


#filter by major clusters
def get_sequence_count(fasta_file):
	"""Counts the number of sequences in a FASTA file."""
	count = 0
	with open(fasta_file, 'r') as f:
		for line in f:
			if line.startswith(">"):
				count += 1
	return count

def filter_clusters_by_sequence_count(input_dir, num_clusters_to_keep):
	"""Keeps only the top N clusters with the highest number of sequences for each locus."""
	odir='MajorClusters_'+num_clusters_to_keep
	os.makedirs(odir)
	numkeep=int(num_clusters_to_keep)
	# Dictionary to store cluster info
	loci_clusters = defaultdict(list)

	# Iterate through the files in the input directory
	for file_name in os.listdir(input_dir):
		if file_name.endswith("duprem_al.fasta"):
			parts = file_name.split('_')
			locus = parts[0]
			cluster = parts[1]
			fasta_file = os.path.join(input_dir, file_name)
			seq_count = get_sequence_count(fasta_file)
			loci_clusters[locus].append((seq_count, file_name))
	
	# Process each locus
	for locus, clusters in loci_clusters.items():
		# Sort clusters by sequence count in descending order
		clusters.sort(reverse=True, key=lambda x: x[0])
		# Keep only the top N clusters
		top_clusters = clusters[:numkeep]
		
		# Move the selected clusters to the output directory
		for _, cluster_file in top_clusters:
			shutil.copy(os.path.join(input_dir, cluster_file), os.path.join(odir, cluster_file))

# get N major clusters from filtered loci
filter_clusters_by_sequence_count(filtdir, args.majorclusters)

clustdir=filtdir+'/MajorClusters_' + args.majorclusters+'/'

os.chdir(clustdir)

#Count filtered sequences
def count_files_in_folder(folder_path):
	try:
		# List all files and directories in the specified folder
		items = os.listdir(folder_path)
		
		# Count only the files (exclude directories)
		file_count = sum(1 for item in items if os.path.isfile(os.path.join(folder_path, item)))
		
		return file_count
	except Exception as e:
		print(f"An error occurred: {e}")
		return 0

outputcount=count_files_in_folder(clustdir)

#Loop to generate consensus references from phased alignments in diploids folder
for align in os.listdir(clustdir):
	if align.endswith('duprem_al.fasta'):
		subprocess.call(["cons %s %sconsensus_reference.fasta" % (align, align[:-15])], shell=True)

#loop throught consensus-references and give them locus cluster reference label
for refc in os.listdir(clustdir):
	if refc.endswith('consensus_reference.fasta'):
		refname=refc.split('_')[0]+'_'+refc.split('_')[1]
		with open(refc, 'r') as refile:
			for line in refile:
				if '>' in line:
					name = '>'+ refname +'\n'
					replaceAll(refc, line, name)


# #concatenate references
subprocess.call(["cat *consensus_reference.fasta > %sLoci_%s_%sMajorClusters_consensusrefs.fasta"% (outputcount, filtstr, args.majorclusters)], shell=True)

#filter phased alignments if present
os.chdir(args.workingdir)

if args.stage2 == 'T':
	#get diploids_phased folder
	phasedir = args.workingdir+'diploids_phased/'
	#make filter folder based on settings in phase directory
	phasefiltdir = phasedir+filtstr+'_MajorClusters_' + args.majorclusters+'/'
	os.makedirs(phasefiltdir)
	for file in os.listdir(clustdir):
		if file.endswith('_al.fasta'):
			clusterid = file.split('_')[0]+'_'+file.split('_')[1]+'_'
			for phasefile in os.listdir(phasedir):
				if phasefile.endswith('_final.fasta'):
					if clusterid in phasefile:
						print('Moving '+phasefile+' to '+phasefiltdir)
						shutil.copy(os.path.join(phasedir, phasefile), os.path.join(phasefiltdir, phasefile))

if args.stage3 == 'T':
	#get diploids_phased folder
	hybdir = args.workingdir+'phaseset/diploidclusters_phased/'
	#make filter folder based on settings in phase directory
	hybfiltdir = hybdir+filtstr+'_MajorClusters_' + args.majorclusters+'/'
	os.makedirs(hybfiltdir)

	for file in os.listdir(clustdir):
		if file.endswith('_al.fasta'):
			clusterid = file.split('_')[0]+'_'+file.split('_')[1]+'_'
			for phasefile in os.listdir(hybdir):
				if phasefile.endswith('_trimmed.fasta'):
					if clusterid in phasefile:
						print('Moving '+phasefile+' to '+hybfiltdir)
						shutil.copy(os.path.join(hybdir, phasefile), os.path.join(hybfiltdir, phasefile))

	for file in os.listdir(clustdir):
		if file.endswith('_al.fasta'):
			clusterid = file.split('_')[0]+'_'+file.split('_')[1]+'_'
			for phasefile in os.listdir(hybdir):
				if phasefile.endswith('differentiated.fasta'):
					if clusterid in phasefile:
						print('Moving '+phasefile+' to '+hybfiltdir)
						shutil.copy(os.path.join(hybdir, phasefile), os.path.join(hybfiltdir, phasefile))

if args.keepal == 'F':

	os.chdir(filtdir)

	for file in os.listdir(filtdir):
		if file.endswith('al.fasta'):
			os.remove(file)

	os.chdir(clustdir)

	for file in os.listdir(clustdir):
		if file.endswith('al.fasta'):
			os.remove(file)

if args.dovcf == 'T':

	os.chdir(args.workingdir)

	#make vcf-folder
	os.makedirs('SORTER2Processor_VCF')
	mcstr=str(args.majorclusters)
	countstr=str(outputcount)
	ref = clustdir+countstr+'Loci_'+filtstr+'_'+mcstr+'MajorClusters_consensusrefs.fasta'
	vcfdir = args.workingdir + 'SORTER2Processor_VCF/'

	#move references to vcf folder
	shutil.copy(ref,vcfdir)

	os.chdir(vcfdir)

	#make reference name strings
	refname = countstr+'Loci_'+filtstr+'_'+mcstr+'MajorClusters_consensusrefs.fasta'

	#Index reference
	subprocess.call(["bwa index %s" % (refname)], shell=True)
	subprocess.call(["samtools faidx %s" % (refname)], shell=True)

	refdir = vcfdir+countstr+'Loci_'+filtstr+'_'+mcstr+'MajorClusters_consensusrefs.fasta'

	#Map Reads to consensus reference for each sample

	for folder in os.listdir(args.workingdir):
		if 'assembly' in folder:
			os.chdir(args.workingdir + folder)
			read = folder[:-8] + 'R1_val_1.fq'
			R2 = folder[:-8] + 'R2_val_2.fq'
			subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (refdir, read, R2, folder[:-8])], shell=True)
			subprocess.call(["samtools sort  %smapreads.sam -o %s" % (folder[:-8], folder[:-8] + 'mapreads.bam')], shell=True)
			os.remove(folder[:-8] + 'mapreads.sam')#Remove sam file
			#Move bam file to workingdir for downstream use in mpileup and call command to generate SNPS
			bampath = args.workingdir + folder + '/' + folder[:-8] + 'mapreads.bam'
			shutil.move(bampath, vcfdir + folder[:-8] + 'mapreads.bam')
			#Get readstats from bam
			os.chdir(args.workingdir)
			statfilename = folder[:-8] + "readstats.txt"
			with open(os.path.join(vcfdir, statfilename), 'a+') as statfile:
				os.chdir(vcfdir)
				statfile.write( folder[:-8] + " Read Statistics" + '\n')
				subprocess.call(["samtools flagstat %s >> %s" % (folder[:-8] + "mapreads.bam", statfilename)], shell=True)
				#get read depth
				subprocess.call(["samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (folder[:-8] + "mapreads.bam", statfilename)], shell=True)
				#get coverage
				subprocess.call(["samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (folder[:-8] + "mapreads.bam", statfilename)], shell=True)
				statfile.close()

	os.chdir(vcfdir)

	#create genotype SNP files
	#Make list of bam files
	subprocess.call(["ls *.bam > bam.filelist"], shell=True)
	subprocess.call(["bcftools mpileup -f %s -b bam.filelist --min-BQ 20 -Ov -d 700  | bcftools call -m -Ov -p .0001 -o %s_allsnps.vcf " % (refname, refname[:-20])], shell=True)
	#subprocess.call(["vcftools --vcf %s_allsnps.vcf --het --out %s_het " % (refname[:-20], refname[:-20])], shell=True)

	#Make dictionary of individuals for readstats and heterozygosity
	HETDICT= {}

	for folder in os.listdir(args.workingdir):
		if folder.endswith("assembly"):
			ind=folder[:-9]
			print(ind)
			if ind not in HETDICT:
				HETDICT[ind]={}

	for readstat in os.listdir(vcfdir):
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

	#get readstat values
	het_values_list = list(HETDICT.values())
	columns = [x for x in het_values_list[0]]
	header = ['Individual']+columns
	os.chdir(vcfdir)
	
	with open('VCF_readstats1.csv', 'w') as outfile:
		writer = csv.writer(outfile)
		writer.writerow(header)
		statlist = list(HETDICT.values())
		stats = statlist[0].keys()
		#samples = columns[0:]
		for ind in HETDICT.keys():
			writer.writerow([ind]+[HETDICT[ind][stat] for stat in stats])

	with open('VCF_readstats1.csv', 'r') as infile:
	    with open('VCF_readstats.csv', 'w') as outfile:
		    data = infile.read()
		    data = data.replace("]", "")
		    data = data.replace("[", "")
		    outfile.write(data)

	for file in os.listdir(vcfdir):
		if file.endswith("readstats1.csv"):
			os.remove(file)

else:
	sys.exit('SORTER2VCF Filtered alignments and generated consensus sequence references for filtered loci, skipped read mapping and vcf')
