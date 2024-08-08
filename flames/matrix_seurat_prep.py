'''
Author: Baylee Christensen, Intern, Ji Research Group
Date: 6/26/24
Desc: Takes in data from flames outpu: transcript count, compares them to features to get gene label, and formats output is a specificed criteria
Environment: python-3
'''
import pandas as pd
import gzip
import argparse
import os
import sys

def remove_ensg(value):
    if value.startswith('ENSG'):
        return '_'.join(value.split('_')[1:])
    else:
        return value

def combine_feature_transcript(feature, transcript_id):
    if transcript_id.startswith('ENSG'):
        transcript_id = '_'.join(transcript_id.split('_')[1:])
    return f"{feature}_{transcript_id}"

def parse_commandline():
    parser = argparse.ArgumentParser(description="Process isoform and transcript count data.")
    parser.add_argument('-t', '--transcript_data', required=True, help='Path to the transcript count data file.')
    parser.add_argument('-f', '--features', required=True, help='Path to the features TSV.gz file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file or directory.')
    args = parser.parse_args()
    return args

def main():
    args = parse_commandline()

    # Read transcript count data
    transcript_data = pd.read_csv(args.transcript_data)

    # Read features data
    with gzip.open(args.features, 'rt') as f:
        features = pd.read_csv(f, sep='\t', header=None, names=['gene_id', 'feature', 'other'])

    # Merge the dataframes based on gene_id column
    merged_data = pd.merge(transcript_data, features, on='gene_id', how='left')

    # Combine the feature column with the transcript_id column
    merged_data['transcript_id'] = merged_data.apply(lambda row: combine_feature_transcript(row['feature'], row['transcript_id']), axis =1)

    # Reorder columns to place 'feature' at the front
    cols = ['feature'] + [col for col in merged_data if col != 'feature']
    merged_data = merged_data[cols]

    # Remove 1st and 3rd column
    merged_data.drop(columns=['feature', 'gene_id'], axis =1, inplace=True)

    # Rename columns by adding '-1' to the end of every column name except the first column
    new_columns = [merged_data.columns[0]] + [col + '-1' for col in merged_data.columns[1:]]
    merged_data.columns = new_columns


    # Determine output file path
    output = args.output
    if os.path.isdir(output):
        output_file = os.path.join(output, 'matrix.csv')
    else:
        output_file = output

    # Save the merged dataframe to CSV
    merged_data.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()
