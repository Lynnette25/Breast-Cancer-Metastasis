################################################################################
# Frances Heredia
# Franco Lab
# Description: This script wrangles the counts file obtained from featureCounts to be used for
#          differential expression analysis.

################################################################################

import pandas as pd

# Skip the first line (comment) and read the real header
df = pd.read_csv("/home/fheredia/RNS-seq/Metastasis_Project/RNA_seq/counts_PE.txt", sep="\t",comment='#')

df.columns = list(map(lambda c: c.split('/')[-1].replace('_Aligned.sortedByCoord.out.bam', '') if c.startswith('/home') else c, df.columns))

print(df.columns)

Index(['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'R18LIV_1',
       'R18LIV_2', 'R18LUN_1', 'R18LUN_2', 'R23LIV_1', 'R23LIV_2', 'R23LUN_1',
       'R23LUN_2', 'R33LIV_1', 'R33LIV_2', 'R33LUN_1', 'R33LUN_2', 'R36LIV_1',
       'R36LIV_2', 'R36LUN_1', 'R36LUN_2', 'R39BR_1', 'R39LIV_1', 'R39LIV_2',
       'R49BR_1', 'R49BR_2', 'R49BR_3', 'R49LIV_1', 'R49LIV_2', 'R49LUN_1',
       'R49LUN_2', 'R8BR_1', 'R8BR_2', 'R8LIV_1', 'R8LIV_2', 'R8LUN_1',
       'R8LUN_2'],
      dtype='object')


df['Chr'] = df['Chr'].astype(str).str.split(';').str[0]
df['Strand'] = df['Strand'].astype(str).str.split(';').str[0]
df['Start'] = df['Start'].astype(str).str.split(';').apply(lambda x: min(map(int, filter(None, x))))
df['End'] = df['End'].astype(str).str.split(';').apply(lambda x: max(map(int, filter(None, x))))


gene_names_df = pd.read_csv("Counts.csv")

# Merge on 'Geneid'
merged_df = gene_names_df[['gene_name', 'Geneid']].drop_duplicates().merge(df, on='Geneid', how='right')

merged_df['gene_name'] = merged_df['gene_name'].fillna(merged_df['Geneid'])

merged_df.to_csv("counts_PE.csv", index=False)