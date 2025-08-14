import sys
import os
import csv

def parse_attributes(attr_str):
    #Parse GFF3 attributes column into a dictionary
    attrs = {}
    for part in attr_str.split(';'):
        if '=' in part:
            key, val = part.split('=', 1)
            attrs[key] = val
    return attrs

def get_feature_name(attrs):
    #Get a name or ID for the gene from attributes
    if 'Name' in attrs:
        return attrs['Name']
    elif 'ID' in attrs:
        return attrs['ID']
    else:
        return "Unnamed_gene"

def makepartition(gff_file):
    # Store gene coordinates as (start, end, name)
    genes = []
    
    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = parts
            start, end = int(start), int(end)
            
            # Only consider 'gene' coordinates spanning more than 1bp
            if ftype == 'gene' and (end - start >= 1):
                attrs = parse_attributes(attrs_str)
                gname = get_feature_name(attrs)
                genes.append((start, end, gname))
    
    # Sort genes by start coordinate
    genes.sort(key=lambda x: x[0])
    
    # Output genes and intergenic intervals
    output_lines = []
    if genes:
        first_gene_start = genes[0][0]
        if first_gene_start > 1:
            # Leading intergenic region before the first gene
            output_lines.append(f"Start_{genes[0][2]}=1-{first_gene_start-1}")
    
    prev_gene = None
    for i, (gstart, gend, gname) in enumerate(genes):
        # Gene line
        output_lines.append(f"{gname}={gstart}-{gend}")
        
        # If not the first gene, print intergenic interval to previous gene
        if prev_gene is not None:
            pgstart, pgend, pgname = prev_gene
            if gstart > pgend + 1:
                # Intergenic region
                intergenic_start = pgend + 1
                intergenic_end = gstart - 1
                output_lines.insert(-1, f"{pgname}_{gname}={intergenic_start}-{intergenic_end}")
        
        prev_gene = (gstart, gend, gname)
    
    # Write the output to partitions.txt
    with open("partitions.txt", "w") as outfile:
        for line in output_lines:
            outfile.write(line + "\n")

            
def read_partitions(partition_file):

    #Reads the partition file which has lines like:
    #region_name=start-end
    #Returns a list of tuples: [(region_name, start, end)
    partitions = []
    with open(partition_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Split on '='
            # Example: geneA=100-200
            # region_name = 'geneA'
            # interval = '100-200'
            parts = line.split('=')
            if len(parts) != 2:
                continue
            region_name = parts[0]
            interval = parts[1]
            if '-' not in interval:
                continue
            start_str, end_str = interval.split('-')
            start, end = int(start_str), int(end_str)
            partitions.append((region_name, start, end))
    return partitions

def read_fasta_alignment(fasta_file):
    #Reads a FASTA alignment file.
    #Returns a dictionary: {sample_name: sequence}
    seqs = {}
    current_name = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # save previous
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

def clean_region_name(name):
    #Remove '-' and '.' from the region name for file naming.#
    return name.replace('-', '').replace('.', '')

def get_genome_id(gff3_file):
    with open(gff3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('#') and line:
                parts = line.split('\t')
                if len(parts) > 0:
                    genome_id = parts[0]
                    # Remove all '.' characters
                    genome_id = genome_id.replace('.', '')
                    return genome_id
    return None  # If no suitable line is found

def splitfasta(partition_file, fasta_file, genomeid):
    # Read partitions
    partitions = read_partitions(partition_file)
    # Read alignment
    seqs = read_fasta_alignment(fasta_file)

    # Create partitions directory if not exists
    if not os.path.exists(genomeid+'_partitions'):
        os.makedirs(genomeid+'_partitions')

    # We'll collect counts for summary: {sample: {region: count_of_ATGC}}
    samples = list(seqs.keys())
    region_names = [r[0] for r in partitions]
    counts = {sample: {r[0]: 0 for r in partitions} for sample in samples}

    # For each region, extract the relevant slice and write a fasta file
    for region_name, start, end in partitions:
        start_idx = start - 1
        end_idx = end

        # Write the FASTA for this region
        cleaned_name = clean_region_name(region_name)
        output_fasta = os.path.join(genomeid+'_partitions', f"{cleaned_name}.fasta")
        with open(output_fasta, 'w') as out_f:
            for sample in samples:
                full_seq = seqs[sample]
                # Extract region
                region_seq = full_seq[start_idx:end_idx]
                out_f.write(f">{sample}\n{region_seq}\n")

                # Count ATGC
                # Only count characters A,T,G,C (uppercase)
                count_ATGC = sum(base in ['A','T','G','C','a','t','g','c'] for base in region_seq)
                counts[sample][region_name] = count_ATGC

    # create the CSV summary
    summary_file = 'partitions_summary_table.csv'
    with open(summary_file, 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        header = ['Sample'] + region_names
        writer.writerow(header)
        for sample in samples:
            row = [sample] + [counts[sample][r] for r in region_names]
            writer.writerow(row)

    print("Done. Partitioned FASTA files are in 'partitions/' and summary is in 'partitions_summary_table.csv'.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <inputgff3> <fasta_alignment>")
        sys.exit(1)
    makepartition(sys.argv[1])
    genomestring = get_genome_id(sys.argv[1])
    partition_file = 'partitions.txt'
    fasta_file = sys.argv[2]
    splitfasta(partition_file, fasta_file, genomestring)
