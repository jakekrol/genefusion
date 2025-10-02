#!/usr/bin/env bash

# commented section was for logreg top 1000 analysis

# # start w/ top 1000
# head -n 1000 ovary_rand_inferred_sort.tsv > rand_inferred_top1000.tsv
# cut -f 1,2 rand_inferred_top1000.tsv > z && mv z rand_inferred_top1000.tsv
# sed -i 's|\t|--|' rand_inferred_top1000.tsv
# # filter using FusionAnnotator
# FusionAnnotator --annotate rand_inferred_top1000.tsv --genome_lib_dir $GENOME_LIB_DIR --max_neighbor_dist 100000 > rand_inferred_top1000_annotated.txt
# grep -v "NEIGHBOR" rand_inferred_top1000_annotated.txt | grep -v "BLASTPAIR" | grep -v "PARALOG" > rand_inferred_top1000_annotated_filt.txt
# # clean up duplicates (caused by neighbors)
# script='import pandas as pd;f="./rand_inferred_top1000_annotated_filt.txt";df=pd.read_csv(f,sep="\t");df["x"]=df["left--right"].apply(lambda x: set(x.split("--")));df.drop_duplicates(subset=["x"],inplace=True);df.drop(columns=["x"],inplace=True);df.to_csv("rand_inferred_top1000_annotated_filt_no_dups.txt",sep="\t",index=False)'
# python -c "$script"

# # get top tissue expression
# sed 's|--|\t|' rand_inferred_top1000_annotated_filt_no_dups.txt > z
# script='import pandas as pd;\
# from genefusion.genefusion import *;\
# df_g2t = pd.read_csv("../../results/2025_10-hpa_gene_expr/gene_tissue_top3.tsv", sep="\t");\
# df = pd.read_csv("z", sep="\t");\
# df["left_tissues"] = df["left"].swifter.apply(lambda x: gene2tissue(x, df_g2t));\
# df["right_tissues"] = df["right"].swifter.apply(lambda x: gene2tissue(x, df_g2t));\
# df.to_csv("top_ovary_annotated_filt_no_dups_expr.tsv", sep="\t", index=False)'
# python -c "$script"

# # finally rejoin the read data
# join.py \
#     -t left \
#     -x ./top_ovary_annotated_filt_no_dups_expr.tsv \
#     -y ../../logreg/ovary/ovary_rand_all_ftr_filt_no_norm.tsv \
#     -k left,right \
#     -o top_ovary_annotated_filt_no_dups_expr_ftr.tsv

# # subset columns
# cut --complement -f 8,9,12,13,16,17 top_ovary_annotated_filt_no_dups_expr_ftr.tsv \
#     > top_ovary_annotated_filt_no_dups_expr_ftr_slim.tsv

# # sort
# script='import pandas as pd;\
# df = pd.read_csv("top_ovary_annotated_filt_no_dups_expr_ftr_slim.tsv", sep="\t");\
# df.sort_values(by=["onekg_dna_pe_count", "ovary_rna_pe_count_tumor","ovary_dna_pe_count"], ascending=[True, False, False], inplace=True);\
# df.to_csv("top_ovary_annotated_filt_no_dups_expr_ftr_slim_sorted.tsv", sep="\t", index=False)'
# python -c "$script"

### sorting raw ftrs
ovary_rand_ftrs='../../logreg/ovary/ovary_random_all_ftr_filt.tsv'
# rm density and R_trans cols
cut --complement -f 5,6,9,10,13,14,17,18 $ovary_rand_ftrs > ovary_rand_ftrs.tsv
# add tissue expr
script='import duckdb;\
import pandas as pd;\
from genefusion.genefusion import *;\
conn = duckdb.connect();\
df_g2t = conn.execute("SELECT * FROM read_csv_auto(\"../../results/2025_10-hpa_gene_expr/gene_tissue_top3.tsv\", delim=\"\t\")").df();\
df = conn.execute("SELECT * FROM read_csv_auto(\"ovary_rand_ftrs.tsv\", delim=\"\t\")").df();\
print(df.shape);\
print(df_g2t.shape);\
unique_genes = pd.concat([df["left"], df["right"]]).unique();\
gene_dict = {gene: gene2tissue(gene, df_g2t) for gene in unique_genes};\
print("mapping left genes to tissues...");\
df["left_tissues"] = df["left"].map(gene_dict);\
print("mapping right genes to tissues...");\
df["right_tissues"] = df["right"].map(gene_dict);\
conn.execute("COPY (SELECT * FROM df) TO \"ovary_rand_ftrs_expr.tsv\" (DELIMITER \"\t\", HEADER)");\
conn.close()'
python -c "$script"

# add a key column (left--right) for annotation
# write to new file
# use duckdb for this
script='import duckdb;\
import pandas as pd;\
conn = duckdb.connect();\
df = conn.execute("SELECT * FROM read_csv_auto(\"ovary_rand_ftrs_expr.tsv\", delim=\"\t\")").df();\
print("making gene_pair column...");\
df["gene_pair"] = df["left"] + "--" + df["right"];\
conn.execute("COPY (SELECT * FROM df) TO \"ovary_rand_ftrs_expr_key.tsv\" (DELIMITER \"\t\", HEADER)");\
conn.close()'
python -c "$script"

# annotate
FusionAnnotator --annotate ovary_rand_ftrs_expr_key.tsv \
    --genome_lib_dir $GENOME_LIB_DIR \
    --max_neighbor_dist 100000 \
    --fusion_name_col 12 \
    > ovary_rand_ftrs_expr_annotated.tsv




cat > score_batches.py << 'EOF'
import duckdb
import numpy as np
from genefusion.genefusion import *

conn = duckdb.connect()
total_rows = conn.execute('SELECT COUNT(*) FROM read_csv_auto("ovary_rand_ftrs_expr_annotated.txt", delim="\t")').fetchone()[0]
print(f"Total rows: {total_rows}")

batch_size = 500000  # Process in smaller batches
for i in range(0, total_rows, batch_size):
    batch_num = i // batch_size + 1
    total_batches = (total_rows + batch_size - 1) // batch_size
    print(f"Processing batch {batch_num}/{total_batches} (rows {i} to {min(i + batch_size, total_rows)})")
    
    df = conn.execute(f'SELECT * FROM read_csv_auto("ovary_rand_ftrs_expr_annotated.txt", delim="\t") LIMIT {batch_size} OFFSET {i}').df()
    df["fusion_score"] = np.vectorize(score)(
        df["ovary_dna_pe_count_normal"], 
        df["ovary_dna_pe_count"], 
        0, 
        df["ovary_rna_pe_count_tumor"], 
        df["onekg_dna_pe_count"], 
        df["ovary_dna_sample_count_normal"], 
        df["ovary_dna_sample_count"], 
        0, 
        df["ovary_rna_sample_count_tumor"], 
        df["onekg_dna_sample_count"], 
        67, 70, 0, 70, 2536,
        w_dna=0.5, 
        upper_factor=100
    )
    
    # Write each batch with proper append mode
    mode = "w" if i == 0 else "a"
    header = i == 0
    df.to_csv("temp_scored.tsv", sep="\t", index=False, mode=mode, header=header)
    print(f"Batch {batch_num} completed")

# Now sort the complete aggregated file using DuckDB
print("Sorting complete file by fusion_score...")
conn.execute('COPY (SELECT * FROM read_csv_auto("temp_scored.tsv", delim="\t") ORDER BY fusion_score DESC) TO "ovary_rand_ftrs_expr_annotated_scored_sorted.tsv" (DELIMITER "\t", HEADER)')

# Clean up temporary file
print("Cleaning up temporary files...")
import os
os.remove("temp_scored.tsv")

conn.close()
print("Processing complete! Output: ovary_rand_ftrs_expr_annotated_scored_sorted.tsv")
EOF

python score_batches.py


# # score
# # use duckdb for this
# script='import duckdb;\
# import numpy as np;\
# from genefusion.genefusion import *;\
# conn = duckdb.connect();\
# # read in annotated ftrs
# print("reading annotated ftrs...");\
# df = conn.execute("SELECT * FROM read_csv_auto(\"ovary_rand_ftrs_expr_annotated.txt\", delim=\"\t\")").df();\
# # score
# print("scoring...");\
# df["fusion_score"] = np.vectorize(score)(df["ovary_dna_pe_count_normal"], df["ovary_dna_pe_count"], 0, df["ovary_rna_pe_count_tumor"], df["onekg_dna_pe_count"], df["ovary_dna_sample_count_normal"], df["ovary_dna_sample_count"], 0, df["ovary_rna_sample_count_tumor"], df["onekg_dna_sample_count"], 67, 70, 0, 70, w_dna=0.5, upper_factor=100);\
# print("writing scored file...");\
# conn.execute("COPY (SELECT * FROM df ORDER BY fusion_score DESC) TO \"ovary_rand_ftrs_expr_annotated_scored_sorted.tsv\" (DELIMITER \"\t\", HEADER)");\
# conn.close()'
# python -c "$script"




# sort 

# # by tumor dna desc, then onekg dna asc
# script='import duckdb;\
# conn = duckdb.connect();\
# conn.execute("COPY (SELECT * FROM read_csv_auto(\"ovary_rand_ftrs_expr_annotated.tsv\", delim=\"\t\") ORDER BY ovary_dna_pe_count DESC, onekg_dna_pe_count ASC) TO \"ovary_rand_ftrs_expr_dna_tumor_sort.tsv\" (DELIMITER \"\t\", HEADER)");\
# conn.close()'
# python -c "$script"

# # by tumor rna desc, then onekg dna asc
# script='import duckdb;\
# conn = duckdb.connect();\
# conn.execute("COPY (SELECT * FROM read_csv_auto(\"ovary_rand_ftrs_expr_annotated.tsv\", delim=\"\t\") ORDER BY ovary_rna_pe_count_tumor DESC, onekg_dna_pe_count ASC) TO \"ovary_rand_ftrs_expr_rna_tumor_sort.tsv\" (DELIMITER \"\t\", HEADER)");\
# conn.close()'
# python -c "$script"

# # rm permutations where left and right are same but different order
# # use duckdb for this
# script='import duckdb;\
# conn = duckdb.connect();\
# conn.execute("COPY (SELECT * FROM read_csv_auto(\"ovary_rand_ftrs_expr.tsv\", delim=\"\t\") WHERE left < right) TO \"ovary_rand_ftrs_expr_no_dups.tsv\" (DELIMITER \"\t\", HEADER)");\
# conn.close()'
# python -c "$script"

# # annotate


