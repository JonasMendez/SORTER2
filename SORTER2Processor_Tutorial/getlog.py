import os
import re
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
args = parser.parse_args()

# Define a regular expression pattern to match log likelihood lines
log_likelihood_pattern = r"Loglikelihood: (-?\d+\.\d+)"

# Create a dictionary to store K values and their associated log likelihoods
results = {}


#Get log likelihoods
for filename in os.listdir(args.workingdir):
    if filename.endswith(".stdout"):
        # Extract the first X value from the filename
        match = re.search(r"\.(\d+)_\d+\.stdout", filename)
        if match:
            x_value = int(match.group(1))
            file_path = os.path.join(args.workingdir, filename)

            try:
                with open(file_path, 'r') as file:
                    lines = file.readlines()

                    # Iterate through the lines in reverse order
                    for line in reversed(lines):
                        match = re.search(log_likelihood_pattern, line)
                        if match:
                            log_likelihood = float(match.group(1))
                            
                            # Store the first found log likelihood
                            if x_value in results:
                                results[x_value].append(log_likelihood)
                            else:
                                results[x_value] = [log_likelihood]
                            
                            # Break after the first match
                            break

            except Exception as e:
                print(f"An error occurred while processing {filename}: {e}")

# Write the results to a tab-delimited text file
output_file = "log_likelihood_results.txt"
with open(output_file, 'w') as out_file:
    out_file.write("K_Value\tLog_Likelihood\n")
    for x, likelihoods in results.items():
        for likelihood in likelihoods:
            out_file.write(f"{x}\t{likelihood}\n")

print(f"LogLikelihood Results saved to {output_file}")


#Make CV data

def extract_values_from_stdout(directory, output_file):
    # Initialize lists to store the K and CVerror values
    k_values = []
    cverror_values = []
    file_names = []
    
    # Regular expression to match the line containing 'CV error (K=X): XVAL'
    pattern = r'CV error \(K=(\d+)\): ([\d.]+)'
    
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.stdout'):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                for line in file:
                    match = re.search(pattern, line)
                    if match:
                        k = int(match.group(1))  # Extract the 'K' value
                        xval = float(match.group(2))  # Extract the 'XVAL' value
                        k_values.append(k)
                        cverror_values.append(xval)
                        file_names.append(filename)
                        break  # If only one match per file is needed, break the loop after finding the first match
    
    # Create a DataFrame with the extracted values
    df = pd.DataFrame({
        'K': k_values,
        'CVerror': cverror_values
    })
    
    # Output the DataFrame to a tab-delimited file
    df.to_csv(output_file, sep='\t', index=False)

# get CV data
output_file2 = args.workingdir+'cv_errors.txt'
extract_values_from_stdout(args.workingdir, output_file2)
print(f"CV Error Results saved to {output_file2}")
