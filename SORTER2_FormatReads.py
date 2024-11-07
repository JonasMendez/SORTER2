import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help='Input table with list of paired fastq files with corresponding IDs and species names (e.g. csv format: Fastq Filename,UniqueID,Species)')
args = parser.parse_args()

direc=os.getcwd()

# #Define command to change sequence IDs
def replaceAll(file,searchExp,replaceExp):
	for line in fileinput.input(file, inplace=1):
		if searchExp in line:
			line = line.replace(searchExp,replaceExp)
		sys.stdout.write(line)

def replace_in_string(file, old_string, new_string):
    # Replace the old_string with new_string in the filename
	new_name = file.replace(old_string, new_string)
	
	return new_name

def changeid(csv_file,direc):
	# Determine the delimiter based on the file extension
	_, file_extension = os.path.splitext(csv_file)

	if file_extension.lower() in ('.tsv', '.txt'):
		delimiter = '\t'
	else:
		delimiter = ','  # Default to comma for CSV and other cases

	#rename files
	with open(os.path.join(direc, csv_file), 'r') as infile:
		for lines in infile:
			for fq in os.listdir(direc):
				if fq in lines:
					line=lines.strip('\n\r')
					#get current R1 filename
					R1=line.split(delimiter)[0]
					#get current R2 filename
					R2=R1.replace('R1','R2')
					#make new R1/R2 filenames
					NewSpeciesID_R1=line.split(delimiter)[1]+'_'+line.split(delimiter)[2]+'_R1.fastq'
					NewSpeciesID_R1b=NewSpeciesID_R1.strip('\n')
					NewSpeciesID_R2=line.split(delimiter)[1]+'_'+line.split(delimiter)[2]+'_R2.fastq'
					NewSpeciesID_R2b=NewSpeciesID_R2.strip('\n')
					print(R1 + '\nNew File Name for SORTER: \n' + NewSpeciesID_R1b)
					print(R2 + '\nNew File Name for SORTER: \n' + NewSpeciesID_R2b)
					#current filepath+name
					old_R1_path = os.path.join(direc, R1)
					old_R2_path = os.path.join(direc, R2)
					#new filepath+name
					new_R1_path = os.path.join(direc, NewSpeciesID_R1b)
					new_R2_path = os.path.join(direc, NewSpeciesID_R2b)
					# Rename the file
					os.rename(old_R1_path, new_R1_path)
					os.rename(old_R2_path, new_R2_path)


#change filenames
changeid(args.input, direc)




