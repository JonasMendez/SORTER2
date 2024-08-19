import os
import sys
import subprocess
import shutil
from os import path
from shutil import move
import fileinput
import csv 
import Bio
from os import path
from shutil import move
from shutil import copyfile
from Bio import SeqIO
import re
import itertools
import argparse
import glob
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-cpref", "--cpref")
parser.add_argument("-c", "--coverage")
parser.add_argument("-d", "--depth")
args = parser.parse_args()

os.chdir(args.workingdir)

subprocess.call(["bwa index %s" % (args.cpref)], shell=True)

#Define command to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

os.chdir(args.workingdir)

#remove previous files
for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir+folder)
		for file in os.listdir(args.workingdir+folder):
			if '_cpreads' in file:
				os.remove(file)
			elif file.endswith("_chloroplastreadstats.txt"):
				os.remove(file)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		dirpath = args.workingdir + folder + '/'
		os.chdir(args.workingdir + folder)
		read = folder[:-8] + 'R1_val_1.fq'
		R2 = folder[:-8] + 'R2_val_2.fq'
		sample=folder[:-9]
		subprocess.call(["bwa mem -V %s %s %s > %s_cpreads.sam" % (args.cpref, dirpath+read, dirpath+R2, sample)], shell=True)
		subprocess.call(["samtools sort  %s -o %s" % (dirpath+sample + '_cpreads.sam', dirpath+sample + '_cpreads_srt.bam')], shell=True)
		subprocess.call(["samtools index  %s" % (dirpath+sample + '_cpreads_srt.bam')], shell=True)
		subprocess.call(["bcftools mpileup -B --min-BQ 20 -Ov -d 500 -f %s %s | bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s " % (args.cpref, dirpath+sample + '_cpreads_srt.bam', dirpath+sample + '_cpreads_called.vcf.gz' )], shell=True)
		subprocess.call(["bcftools index  %s" % (dirpath+sample + '_cpreads_called.vcf.gz')], shell=True)
		subprocess.call(["cat %s | bcftools consensus %s > %scp_final.fasta" % (args.cpref, dirpath+sample+'_cpreads_called.vcf.gz', folder[:-8])], shell=True)
		os.remove(sample+"_cpreads.sam")
		for bam in os.listdir(args.workingdir+folder):
			if bam.endswith('_cpreads_srt.bam'):
				statfilename = bam[:-16] + "_chloroplastreadstats.txt"
				with open(os.path.join(args.workingdir + folder, statfilename), 'a+') as statfile:
					statfile.write( folder[:-9] + " Read Statistics" + '\n')
					subprocess.call(["samtools flagstat %s >> %s" % (sample + '_cpreads_srt.bam', statfilename)], shell=True)
					subprocess.call(["samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (sample + '_cpreads_srt.bam', statfilename)], shell=True)
					subprocess.call(["samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (sample + '_cpreads_srt.bam', statfilename)], shell=True)
					statfile.close()

os.chdir(args.workingdir)

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir+folder)
		for chlor in os.listdir(args.workingdir+folder):
				if chlor.endswith("cp_final.fasta"):
					with open(chlor, 'r') as infile:
						for line in infile:
							if '>' in line:
    								print(line)
    								name ='>' + folder[:-9] + '\n'
    								print(name)
    								replaceAll(chlor, line, name)
os.chdir(args.workingdir)

os.mkdir("all_chloroplasts")

for folder in os.listdir(args.workingdir):
	if folder.endswith("assembly"):
		os.chdir(args.workingdir+folder)
		for chlor in os.listdir(args.workingdir+folder):
			if chlor.endswith('cp_final.fasta'):
				chlorsrc = args.workingdir + folder + '/' + chlor
				chlorst = args.workingdir + 'all_chloroplasts/'+chlor
				os.rename(chlorsrc, chlorst)


os.chdir(args.workingdir + 'all_chloroplasts/')

subprocess.call("cat *cp_final.fasta > ALLsamples_chloroplasts.fasta", shell=True)


os.chdir(args.workingdir)

#compile read statistics
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
    with open('readstats_cp.csv', 'w') as outfile:
	    data = infile.read()
	    data = data.replace("]", "")
	    data = data.replace("[", "")
	    outfile.write(data)

for file in os.listdir(args.workingdir):
	if file.endswith("readstats.csv"):
		os.remove(file)

outputdir = str("coverage"+args.coverage+"_depth"+args.depth+"_filtered")

print("Filtering output sequences based on filters to: "+ outputdir)

def filter_fasta_files(csv_file, fasta_dir, min_coverage, min_readdepth, output_dir=outputdir):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    
    # Ensure the necessary columns exist in the CSV file
    if 'coverage' not in df.columns or 'readdepth' not in df.columns:
        raise ValueError("CSV file must contain 'coverage' and 'readdepth' columns")
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        sample_id = row[df.columns[0]]  # Assuming the first column is the ID
        coverage = row['coverage']
        readdepth = row['readdepth']
        
        # Check if the sample meets the criteria
        if coverage >= int(min_coverage) and readdepth >= int(min_readdepth):
            # Construct the expected FASTA file name
            fasta_filename = f"{sample_id}_cp_final.fasta"
            fasta_path = os.path.join(fasta_dir, fasta_filename)
            
            # Check if the FASTA file exists
            if os.path.exists(fasta_path):
                # Copy the file to the filtered directory
                shutil.copy(fasta_path, output_dir)
                print(f"Copied {fasta_filename} to {output_dir}")
            else:
                print(f"FASTA file {fasta_filename} not found in {fasta_dir}")

# Example usage:
csv_file = args.workingdir+'readstats_cp.csv'
fasta_dir = args.workingdir+'all_chloroplasts/'

os.chdir(args.workingdir+'all_chloroplasts/')

filter_fasta_files(csv_file, fasta_dir, args.coverage, args.depth)

os.chdir(args.workingdir+'all_chloroplasts/'+outputdir)

#concatenate filtered sequences
subprocess.call(["cat *final.fasta > Genomes_%s.fasta" % (outputdir)], shell=True)


