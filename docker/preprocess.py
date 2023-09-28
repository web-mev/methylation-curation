"""
This script creates a simplified version of the probe-mapping file
distributed by Illumina.
"""
import argparse
from tokenize import group
import pandas as pd
import numpy as np


def parse_cl_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_file', required=True)
    parser.add_argument('-o', '--output', dest='output_file', required=True)
    return parser.parse_args()


def parse_delim(s):
    if pd.isna(s):
        return None
    return [x.strip() for x in s.split(';')]


def extract_feature_mapping(row):
    """
    Takes a row from the the dataframe (a pd.Series) and
    returns a dataframe containing only the required information.

    As an example, a row could look like:
    IlmnID                                                          cg27416437
    ...                                                                    ...
    UCSC_RefGene_Name        EWSR1;EWSR1;EWSR1;EWSR1;RHBDD3;EWSR1;EWSR1;EWS...
    UCSC_RefGene_Accession   NM_001163286;NM_013986;NM_001163286;NM_0011632...
    UCSC_RefGene_Group       5'UTR;5'UTR;1stExon;1stExon;TSS200;1stExon;5'U...
    ...                                                                    ...
    Enhancer                                                               NaN
    ...                                                                    ...

    and we return a dataframe that looks like:
    probe          gene       feature        enhancer
    cg27416437     EWSR1      1stExon        NaN
    cg27416437     EWSR1      5'UTR          NaN
    cg27416437     RHBDD3     TSS200         NaN
    """
    probe_id = row['IlmnID']
    gene_list = parse_delim(row['UCSC_RefGene_Name'])
    group_list = parse_delim(row['UCSC_RefGene_Group'])
    is_enhancer = False if pd.isna(row['Enhancer']) else True
    try:
        if gene_list is not None:
            df = pd.DataFrame({'gene_id': gene_list, 'feature': group_list})
        else:
            return None
    except Exception as ex:
        print('x'*200)
        print(f'Unexpected entry for probe {probe_id}:')
        print(row)
        print('x'*200)
        raise ex
    df.drop_duplicates(inplace=True)
    df['enhancer'] = is_enhancer
    df['probe_id'] = probe_id
    return df


def main(input_file, output_file):
    df = pd.read_csv(input_file, skiprows=7)
    reformatted = df.apply(extract_feature_mapping, axis=1)
    reformatted = pd.concat(reformatted.tolist()).reset_index(drop=True)
    reformatted.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    args = parse_cl_args()
    main(args.input_file, args.output_file)
