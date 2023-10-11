import argparse
import sys
import pandas as pd


# name of the column in the probe mapping file which
# dictates where in the gene each probe lands
FEATURE_COL = 'feature'

# name of the column in the probe mapping file which 
# contains a boolean indicating whether the probe is
# located in an enhancer region:
ENHANCER_COL = 'enhancer'

# name of the column in the probe mapping file which
# contains the Illumina cgXXXXXXX probe ID
PROBE_COL = 'probe_id'

# name of the gene symbol column in the probe mapping file.
GENE_COL = 'gene_id'

# key-value pairs for the short name of the platform
# and the associated probe mapping file
PROBE_MAP_FILES = {
    'HM450': '/opt/software/resources/reformatted_probe_mapping.hm450.tsv',
}


def check_features(mapping_df, feature_set):
    """
    Check that the values in feature_set are a subset
    of available 'features' in the mapping file.
    """
    available_features = set(mapping_df[FEATURE_COL].unique())
    diff_set = feature_set.difference(available_features)
    if len(diff_set) > 0:
        raise Exception('The following features are invalid:'
                        f' {", ".join(diff_set)}. Options are'
                        f' {", ".join(available_features)}.')


def parse_cl_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--input',
        help='The full input methylation matrix at the probe level',
        dest='input_file',
        required=True
    )
    parser.add_argument(
        '-o',
        '--output',
        help='Name of the output file',
        dest='output_file',
        required=True
    )
    parser.add_argument(
        '-p',
        '--platform',
        choices=PROBE_MAP_FILES.keys(),
        help='Name of the chip platform',
        dest='platform',
        required=True
    )
    parser.add_argument(
        '-e',
        '--enhancer',
        help='Require to the probe to be annotated in an enhancer region?',
        action='store_true'
    )
    parser.add_argument(
        'features',
        help='Which features to use?',
        nargs='+'
    )
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_cl_args()
    mapping_df = pd.read_table(PROBE_MAP_FILES[args.platform])
    feature_set = set(args.features)
    check_features(mapping_df, feature_set)

    # read the probe-level methylation matrix:
    meth_matrix = pd.read_table(args.input_file, index_col=0)

    # filter the probe mapping file to keep only those features
    # requested:
    mapping_df = mapping_df.loc[mapping_df[FEATURE_COL].isin(feature_set)]

    if args.enhancer:
        mapping_df = mapping_df.loc[mapping_df[ENHANCER_COL]]

    merged_df = pd.merge(mapping_df,
                         meth_matrix,
                         how='inner',
                         left_on=PROBE_COL,
                         right_index=True)

    if merged_df.shape[0] == 0:
        sys.stderr.write('After filtering, the resulting matrix was empty.'
                         ' This can sometimes happen if the probe identifiers'
                         ' are from a platform that we do not support.')
        sys.exit(1)
    else:
        keep_cols = [GENE_COL] + list(meth_matrix.columns)
        merged_df[keep_cols].to_csv(args.output_file, sep='\t', index=False)
