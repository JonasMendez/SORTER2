import os
import csv
import sys
from glob import glob

def read_fasta(filename):
    #Read a FASTA file and return a dict {sample_name: sequence}.
    seqs = {}
    current_name = None
    current_seq = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if any
                if current_name is not None:
                    seqs[current_name] = ''.join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        # Save last
        if current_name is not None:
            seqs[current_name] = ''.join(current_seq)

    return seqs

def count_ATGC(seq):
    #Count the number of A/T/G/C (case-insensitive) in seq.
    count = 0
    for base in seq:
        b = base.upper()
        if b in ['A','T','G','C']:
            count += 1
    return count

def main():
    # Optionally, you could allow an argument for the directory to search.
    # For now, we assume current directory.
    search_dir = '.'

    # Find all folders ending with "_partition"
    partition_dirs = [d for d in os.listdir(search_dir) if os.path.isdir(d) and d.endswith('_partitions')]

    if not partition_dirs:
        print("No directories ending with '_partitions' found.")
        sys.exit(1)

    # region_data structure: {region_filename: {sample_name: sequence}}
    region_data = {}

    # Iterate over each partition folder
    for pdir in partition_dirs:
        fasta_files = glob(os.path.join(pdir, '*.fasta'))
        for ffile in fasta_files:
            region_name = os.path.basename(ffile)  # e.g. geneA.fasta
            seqs = read_fasta(ffile)
            if region_name not in region_data:
                region_data[region_name] = {}
            # Add (or update) sequences
            for sample, seq in seqs.items():
                # If a sample appears multiple times, you can handle duplicates here if needed.
                # For now, we just add/overwrite.
                region_data[region_name][sample] = seq

    # Create combined_partitions directory
    if not os.path.exists('combined_partitions'):
        os.makedirs('combined_partitions')

    # Now output combined FASTA for each region
    # Also prepare data for summary table
    # First gather all samples across all regions
    all_samples = set()
    for rname in region_data:
        for sname in region_data[rname]:
            all_samples.add(sname)
    all_samples = sorted(all_samples)

    # Count dictionary: {sample: {region: atgc_count}}
    counts = {s: {} for s in all_samples}

    for region_name, sample_dict in region_data.items():
        out_path = os.path.join('combined_partitions', region_name)
        with open(out_path, 'w') as out_f:
            # Write all samples for this region
            for sample, seq in sample_dict.items():
                out_f.write(f">{sample}\n{seq}\n")

        # Count ATGC for each sample in this region
        for s in all_samples:
            seq = sample_dict.get(s, "")
            atgc_count = count_ATGC(seq)
            counts[s][region_name] = atgc_count

    # Write combined summary CSV
    region_names = sorted(region_data.keys())
    summary_file = 'combined_partitions_summary.csv'
    with open(summary_file, 'w', newline='') as csv_f:
        writer = csv.writer(csv_f)
        header = ['Sample'] + region_names
        writer.writerow(header)
        for s in all_samples:
            row = [s] + [counts[s].get(r, 0) for r in region_names]
            writer.writerow(row)

    print("Done. Combined FASTA files in 'combined_partitions', summary in 'combined_partitions_summary.csv'.")

if __name__ == "__main__":
    main()
