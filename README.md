![image](https://github.com/user-attachments/assets/db909b0f-d6f1-4880-ba71-f923332bb44c)

# Genome2Vec
A command-line tool for generating genome feature embeddings based on input BED files and pre-defined annotation datasets.

## Installation
copy all file dir with annotation files `./anno_data/`.

make sure your `python` and `bedtools` software are updated and callable.


## Usage
```bash
python genome2vec.py -a input.bed -b output_genome2vec.bed
```

## input
Query region bedgraph, must be sorted using sort -k 1,1 -k2,2n before.

The input file should be a BED6 file format with `.bed` file name. it should contains more than 7 columns as:

'chr', 'start', 'end', 'name', 'score', 'strand', 'value_1', 'value_2',..., 'value_n'

You can use place holder like "." for this columns.

## output
The output file will contain all query columns and then all annotations columns.
