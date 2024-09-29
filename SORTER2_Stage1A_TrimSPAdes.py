import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--projname", required=True, help='Name for SORTER Run')
parser.add_argument("-spades", "--spades", required=True, help='Run Spades Assembly? (T/F)')
parser.add_argument("-trim", "--trim", required=True, help='Run Trim Galore? (T/F)')
args = parser.parse_args()


#working directory iterator
rootwd=os.getcwd()+'/'

#Make Project Folder
dst = rootwd + 'SORTER2_'+ args.projname +'_/'
wdlist = os.listdir(rootwd)
print('Reads will be processed in:\n' + dst)

if args.trim == 'T':
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
			subprocess.run(["trim_galore --quality 20 --length 30 --paired --fastqc %s %s" % (file, file.replace('_R1.','_R2.'))], shell=True)
			os.remove(file)
			os.remove(file.replace('_R1.','_R2.'))
			os.chdir(rootwd)
		else:
			continue

#Assemble Contigs with Spades
if args.spades == 'T':

	for file in os.listdir(wdlist):
		if 'assembly' in file:
			print("Running Spades on " + file)
			os.chdir(rootwd + file)
			for read in os.listdir(rootwd + file):
				if 'R1_val_1.fq' in read:
					R2 = read[:-11] + 'R2_val_2.fq'
					subprocess.call(["spades.py --only-assembler -1 %s -2 %s -o spades_hybrid_assembly" % (read, R2)], shell=True)
					os.chdir(rootwd)
				else:
					continue
		else:
			continue


sys.exit("Trimgalore and SPADES processing has finished, exiting script")
