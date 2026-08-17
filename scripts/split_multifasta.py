import os
import argparse
from Bio import SeqIO

def split_multifasta(input_file, output_dir, num_files):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the multifasta file
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Calculate the number of sequences per file
    sequences_per_file = len(sequences) // num_files
    remainder = len(sequences) % num_files

    # Split the sequences into N files
    start_idx = 0
    for i in range(num_files):
        end_idx = start_idx + sequences_per_file + (1 if i < remainder else 0)
        output_file = os.path.join(output_dir, f"p{i + 1}.fasta")

        # Write sequences to the output file
        with open(output_file, "w") as output_handle:
            SeqIO.write(sequences[start_idx:end_idx], output_handle, "fasta")

        start_idx = end_idx

def main():
    parser = argparse.ArgumentParser(description="Split a multifasta file into N files with approximately equal number of sequences.")
    parser.add_argument("input_file", help="Path to the input multifasta file.")
    parser.add_argument("output_dir", help="Path to the output directory.")
    parser.add_argument("num_files", type=int, help="Number of output files to create.")

    args = parser.parse_args()
    split_multifasta(args.input_file, args.output_dir, args.num_files)

if __name__ == "__main__":
    main()

