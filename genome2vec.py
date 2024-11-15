import pybedtools
import pandas as pd
import os
from io import StringIO
import argparse

class Genome2Vec:
    def __init__(self, input_file, output_file=None, anno_dir='./anno_data'):
        self.input_file = input_file
        self.output_file = output_file or f"{os.path.splitext(input_file)[0]}_genome2vec.bed"
        self.anno_dir = anno_dir
        self.load_annotations()
    
    def load_annotations(self):
        """Load all necessary annotation files."""
        self.scGPT_emb = pybedtools.BedTool(os.path.join(self.anno_dir, 'gene_embed.bed'))
        self.chromHMM_emb = pybedtools.BedTool(os.path.join(self.anno_dir, 'chromHMM_200bp_UMAPembed.bed'))
        self.is_value = pybedtools.BedTool(os.path.join(self.anno_dir, '40k_is.sort.bed'))
        self.di_value = pybedtools.BedTool(os.path.join(self.anno_dir, '40k_di.sort.bed'))
        self.fi_value = pybedtools.BedTool(os.path.join(self.anno_dir, '40k_fire.sort.bed'))
        self.ab_value = pybedtools.BedTool(os.path.join(self.anno_dir, '250k_hesc_ab.sort.bed'))
        self.hic_value = pybedtools.BedTool(os.path.join(self.anno_dir, '20k_hic.sort.bed'))
        
    def process_bed_files(self):
        """Process BED files and get the closest annotations."""
        input_bed = pybedtools.BedTool(self.input_file)
        # Check if the input file has more than 7 columns
        if input_bed.field_count() < 7:
            raise ValueError(
                "The input file should be a BED6 file format with `.bed` file name. "
                "It should contain at least 7 columns as: 'chr', 'start', 'end', 'name', 'score', "
                "'strand', 'value_1', 'value_2', ..., 'value_n'. "
                "You can use placeholders like '.' for missing columns."
            )

        self.scGPT_bed = input_bed.closest(self.scGPT_emb, t='first', d=False)
        self.chromHMM_bed = input_bed.closest(self.chromHMM_emb, t='first', d=False)
        self.is_bed = input_bed.closest(self.is_value, t='first', d=False)
        self.di_bed = input_bed.closest(self.di_value, t='first', d=False)
        self.fi_bed = input_bed.closest(self.fi_value, t='first', d=False)
        self.ab_bed = input_bed.closest(self.ab_value, t='first', d=False)
        self.hic_bed = input_bed.closest(self.hic_value, t='first', d=False)
    
    @staticmethod
    def read_bedtool_to_df(bedtool_obj):
        """Convert a BedTool object to a pandas DataFrame."""
        bedtool_str = str(bedtool_obj)
        bedtool_df = pd.read_csv(StringIO(bedtool_str), sep="\t", header=None)
        return bedtool_df

    def annotate_data(self):
        """Load BedTool objects, convert them to pandas DataFrames, and annotate data."""
        scGPT_df = self.read_bedtool_to_df(self.scGPT_bed)
        chromHMM_df = self.read_bedtool_to_df(self.chromHMM_bed)
        is_df = self.read_bedtool_to_df(self.is_bed)
        di_df = self.read_bedtool_to_df(self.di_bed)
        fi_df = self.read_bedtool_to_df(self.fi_bed)
        ab_df = self.read_bedtool_to_df(self.ab_bed)
        hic_df = self.read_bedtool_to_df(self.hic_bed)

        # Annotate columns
        scGPT_n_cols = len(scGPT_df.columns)
        scGPT_columns = (
            ['query_chr', 'query_start', 'query_end', 'query_name', 'query_score', 'query_strand'] +
            [f'query_value_{i+1}' for i in range(scGPT_n_cols - 523)] +
            ['near_gene_chr', 'near_gene_start', 'near_gene_end', 'near_gene_name', 'near_gene_strand'] +
            [f"scGPT_emb_{i}" for i in range(1, 513)]
        )
        scGPT_df.columns = scGPT_columns
        scGPT_df['TSS'] = scGPT_df.apply(lambda row: row['near_gene_start'] if row['near_gene_strand'] == '+' else row['near_gene_end'], axis=1)
        scGPT_df['dist_TSS'] = (scGPT_df['query_start'] + scGPT_df['query_end']) / 2 - scGPT_df['TSS']

        chromHMM_df.columns = list(chromHMM_df.columns[:-5]) + ['chromHMM_name', 'chromHMM_UMAPemb_1', 'chromHMM_UMAPemb_2', 'chromHMM_UMAPemb_3', 'chromHMM_UMAPemb_4']
        is_df.columns = list(is_df.columns[:-1]) + ['is_value']
        di_df.columns = list(di_df.columns[:-1]) + ['di_value']
        fi_df.columns = list(fi_df.columns[:-1]) + ['fi_value']
        ab_df.columns = list(ab_df.columns[:-1]) + ['ab_value']
        hic_df.columns = list(hic_df.columns[:-6]) + ['hic_matx', 'hic_maty', 'hic_matz', 'hic_fatx', 'hic_faty', 'hic_fatz']

        # Concatenate data
        out_df = pd.concat([
            scGPT_df, chromHMM_df.iloc[:, -5:], is_df.iloc[:, -1:], di_df.iloc[:, -1:],
            fi_df.iloc[:, -1:], ab_df.iloc[:, -1:], hic_df.iloc[:, -6:]
        ], axis=1)
        
        # Save output
        out_df.to_csv(self.output_file, sep='\t', index=False)
        print(f"Output saved to {self.output_file}")

    def run(self):
        self.process_bed_files()
        self.annotate_data()

def main():
    parser = argparse.ArgumentParser(description="Genome2Vec: Generate genome feature embeddings. -a input_more_than_seven_column.bed -b output_with_header.bed, bedtools > 2.2")
    parser.add_argument("-a", "--input", required=True, help="Input BED file, BED6 format with one or more extra value column ")
    parser.add_argument("-b", "--output", help="Output BED file with header(optional)")
    args = parser.parse_args()
    
    print('preparing...')
    genome2vec = Genome2Vec(input_file=args.input, output_file=args.output)
    genome2vec.run()

if __name__ == "__main__":
    main()
