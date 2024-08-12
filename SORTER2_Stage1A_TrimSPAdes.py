import os
import sys
import argparse
import subprocess
from datetime import date

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--projname", required=True, help='Name for SORTER Run')
parser.add_argument("-spades", "--spades", required=True, help='Run Spades Assembly? (T/F)')
args = parser.parse_args()

#working directory iterator
rootwd=os.getcwd()+'/'
wdlist = os.listdir(rootwd)

# Get the current date
current_date = date.today()
#Format Date
fdate = current_date.strftime("%d%m%Y")

#Make Project Folder
dst = rootwd + 'SORTER2_'+ args.projname +'_'+fdate+'/'

print('Reads will be processed in:\n' + dst)

#Trim reads, put paired reads in folders corresponding to each sample

for file in wdlist:
	if 'R1' in file:
		filefolder=file.split('_')[0]+'_'+file.split('_')[1]+'_assembly'
		readst= dst + filefolder + '/'
		os.makedirs(os.path.join(readst))
		print('Moving ' + file + ' to: \n' + readst)
		R1_old_path = os.path.join(rootwd, file)
		R1_new_path = os.path.join(readst, file)
		R2_old_path = os.path.join(rootwd, file.replace('_R1.','_R2.'))
		R2_new_path = os.path.join(readst, file.replace('_R1.','_R2.'))
		os.renames(R1_old_path, R1_new_path)
		os.renames(R2_old_path, R2_new_path)
		os.chdir(readst)
		subprocess.call(["trim_galore --quality 20 --length 30 --paired --fastqc %s %s" % (file, file.replace('_R1.','_R2.'))], shell=True)
		os.remove(file)
		os.remove(file.replace('_R1.','_R2.'))
		os.chdir(rootwd)
	else:
		continue

#Assemble Contigs with Spades
if args.spades == 'T':

	for file in os.listdir(dst):
		if 'assembly' in file:
			print("Running Spades on " + file)
			os.chdir(dst + file)
			for read in os.listdir(dst + file):
				if 'R1_val_1.fq' in read:
					R2 = read[:-11] + 'R2_val_2.fq'
					subprocess.call(["spades.py --only-assembler -1 %s -2 %s -o spades_hybrid_assembly" % (read, R2)], shell=True)
					os.chdir(dst)
				else:
					continue
		else:
			continue


sys.exit("Trimgalore and SPADES processing has finished, exiting script")
