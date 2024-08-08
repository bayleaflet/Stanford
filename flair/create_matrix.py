'''
Author: Baylee Christensen, Intern, Ji Research Group
Date: 8/1/2024

'''
import csv
from collections import defaultdict
import argparse

def parse_line(line):
    """Parse a line from the input file into a header and list of cell barcodes."""
    header, entries = line.strip().split('\t')
    entries = entries.split(',')
    cell_barcodes = [entry.split('#')[0].split('_')[0] for entry in entries]
    return header, cell_barcodes

def build_count_matrix(infile_path):
    """Build a count matrix from the input file."""
    count_matrix = defaultdict(lambda: defaultdict(int))

    with open(infile_path, 'r') as infile:
        for line in infile:
            header, cell_barcodes = parse_line(line)
            for barcode in cell_barcodes:
                count_matrix[header][barcode] += 1

    return count_matrix

def get_all_barcodes(count_matrix):
    """Extract all unique cell barcodes from the count matrix."""
    all_barcodes = sorted({barcode for sub_dict in count_matrix.values() for barcode in sub_dict})
    return all_barcodes

def write_output(count_matrix, all_barcodes, outfile_path):
    """Write the count matrix to a CSV file."""
    with open(outfile_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Gene'] + all_barcodes)

        for header, barcode_counts in count_matrix.items():
            row = [header] + [barcode_counts[barcode] for barcode in all_barcodes]
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Create a count matrix from input file.")
    parser.add_argument('-i', '--input', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')
    args = parser.parse_args()

    count_matrix = build_count_matrix(args.input)
    all_barcodes = get_all_barcodes(count_matrix)
    write_output(count_matrix, all_barcodes, args.output)

if __name__ == "__main__":
    main()


