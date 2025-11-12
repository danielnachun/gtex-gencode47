import pandas as pd
import seaborn as sns
import pyranges as pr
import numpy as np
import bisect

from tqdm import tqdm


tqdm.pandas()


exons = pd.read_csv(
    "/oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/exons_to_match.bed",
    sep="\t",
    header=None,
    names=["chr", "start", "end", "exon_id", "gene_id"],
)
exons["length"] = exons["end"] - exons["start"]
exons = exons[~exons["chr"].isin(["chrM"])]

background = pd.read_csv(
    "/oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/intergenic_regions.bed",
    sep="\t",
    header=None,
    names=["chr", "start", "end"],
)

background = background[background["chr"].isin(exons["chr"].unique())]
background = background.sort_values(["chr", "start"])


gencode_v47_path = "/oak/stanford/groups/smontgom/dnachun/data/dgtex/dGTEx_Consortium_resources/NHP-dGTEx_references/Enhanced_GTFs/Marmosets/v1/Marmoset_mCalJa1.2.pat.X_enhanced_v1_genes_collapsed_stranded.gtf"
gencode_v47 = pr.read_gtf(gencode_v47_path)
gencode_v47_df = gencode_v47.as_df()
gencode_v47_df = gencode_v47_df[~(gencode_v47_df["Chromosome"] == "chrM")]


def get_closest_background_no_overlap(exon, background_starts, background_ends):
    exon_length = exon["end"] - exon["start"]
    idx_right = bisect.bisect_left(background_starts, exon["start"])
    while idx_right < len(background_starts) and background_starts[idx_right] < exon["start"]:
        idx_right += 1

    idx_left = idx_right - 1
    max_iter = len(background_starts)
    iter_count = 0

    while iter_count < max_iter:
        iter_count += 1
        if idx_right < len(background_starts):
            right_start, right_end = background_starts[idx_right], background_ends[idx_right]
            if right_end - right_start >= exon_length:
                background_starts[idx_right] = background_starts[idx_right] + exon_length
                return [right_start, right_start + exon_length]

        if idx_left >= 0:
            left_start, left_end = background_starts[idx_left], background_ends[idx_left]
            if left_end - left_start >= exon_length:
                background_ends[idx_left] = background_ends[idx_left] - exon_length
                return [left_end - exon_length, left_end]

        idx_right += 1
        idx_left -= 1
    
    # Return None if no match found
    return None


null_exons = []
unmatched_exons = []
for chrom, background_chrom in background.groupby("chr"):
    null_exons_chr = []
    print(chrom)
    background_starts = background_chrom["start"].values
    background_ends = background_chrom["end"].values
    exon_chrom = exons[exons["chr"] == chrom]
    print(len(exon_chrom))
    exon_ids = []
    gene_ids = []
    for idx, exon in tqdm(exon_chrom.iterrows(), total=len(exon_chrom)):
        match_result = get_closest_background_no_overlap(exon, background_starts, background_ends)
        if match_result is not None:
            null_exons_chr.append(match_result)
            exon_ids.append(exon["exon_id"])
            gene_ids.append(exon["gene_id"])
        else:
            # Track unmatched exons
            unmatched_exons.append({
                "chr": chrom,
                "start": exon["start"],
                "end": exon["end"],
                "exon_id": exon["exon_id"],
                "gene_id": exon["gene_id"],
                "length": exon["length"]
            })
    
    null_exons_chr_df = pd.DataFrame(null_exons_chr, columns=["start", "end"])
    null_exons_chr_df["exon_id"] = exon_ids
    null_exons_chr_df["gene_id"] = gene_ids
    null_exons_chr_df["chr"] = chrom
    null_exons.append(null_exons_chr_df)

null_exons_all_chrs = pd.concat(null_exons)
null_exons_all_chrs[["chr", "start", "end", "exon_id", "gene_id"]].to_csv(
    "/oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/matched_nongenic_null_intergenic_no_overlap.bed",
    index=False,
    header=None,
    sep="\t",
)

# Track which genes have all exons matched
# Count exons per gene in original data
exon_counts_per_gene = exons.groupby("gene_id").size()
# Count matched exons per gene
matched_exon_counts_per_gene = null_exons_all_chrs.groupby("gene_id").size()

# Only keep genes where all exons were matched
genes_with_all_exons_matched = set(
    gene_id 
    for gene_id in exon_counts_per_gene.index 
    if gene_id in matched_exon_counts_per_gene.index 
    and exon_counts_per_gene[gene_id] == matched_exon_counts_per_gene[gene_id]
)

# Get all genes that were in the original exon list
all_genes_in_exons = set(exons["gene_id"].unique())

# Find genes that were discarded (have at least one unmatched exon)
discarded_genes = all_genes_in_exons - genes_with_all_exons_matched

print(f"Total genes in original exon list: {len(all_genes_in_exons)}")
print(f"Genes with all exons matched: {len(genes_with_all_exons_matched)}")
print(f"Genes discarded (have unmatched exons): {len(discarded_genes)}")
print(f"Total unmatched exons: {len(unmatched_exons)}")

# Save discarded genes to a file
discarded_genes_df = pd.DataFrame({"gene_id": sorted(discarded_genes)})
discarded_genes_df.to_csv(
    "/oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/discarded_genes_with_unmatched_exons.txt",
    index=False,
    header=False,
)
print(f"Discarded genes saved to: /oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/discarded_genes_with_unmatched_exons.txt")

# Save unmatched exons to a file
if unmatched_exons:
    unmatched_exons_df = pd.DataFrame(unmatched_exons)
    unmatched_exons_df[["chr", "start", "end", "exon_id", "gene_id"]].to_csv(
        "/oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/unmatched_exons.bed",
        index=False,
        header=None,
        sep="\t",
    )
    print(f"Unmatched exons saved to: /oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset/unmatched_exons.bed")
else:
    print("No unmatched exons found.")

null_exons_all_chrs["exon_id"] = null_exons_all_chrs["exon_id"].str.split(";").str[0]

null_exons_gtf = pd.merge(
    gencode_v47_df[gencode_v47_df["Feature"] == "exon"].reset_index(),
    null_exons_all_chrs,
    on="exon_id",
    suffixes=["", "_null"],
)

# Filter to only include genes where all exons were matched
null_exons_gtf = null_exons_gtf[null_exons_gtf["gene_id"].isin(genes_with_all_exons_matched)]

null_full_gtf = pd.merge(
    gencode_v47_df.reset_index(),
    null_exons_gtf[["index", "start", "end"]],
    on="index",
    how="left",
    suffixes=["", "_null"],
)

# Filter to only include genes where all exons were matched
null_full_gtf = null_full_gtf[null_full_gtf["gene_id"].isin(genes_with_all_exons_matched)]

gene_starts = null_full_gtf.groupby("gene_id").agg({"start": "min", "end": "max"})
null_full_gtf = pd.merge(
    null_full_gtf, gene_starts, on="gene_id", how="left", suffixes=["", "_gene"]
)

null_full_gtf["Start"] = np.where(
    null_full_gtf["Feature"] == "exon", null_full_gtf["start"], null_full_gtf["start_gene"]
)
null_full_gtf["End"] = np.where(
    null_full_gtf["Feature"] == "exon", null_full_gtf["end"], null_full_gtf["end_gene"]
)


def df_to_gtf(df, output_file):
    with open(output_file, "w") as f:
        f.write(
            "##description: nongenic nulls matched to exons from annotation of the human genome (GRCh38), version 47 (Ensembl 113)\n##provider: Kate Lawrence\n##contact: klawren@stanford.edu\n##format: gtf\n##date: 2025-03-21"
        )
        for index, row in df.iterrows():
            attributes = (
                f'gene_id "{row["gene_id"]}"; transcript_id "{row["transcript_id"]}"; exon_id "{row["exon_id"]}";'
            )
            gtf_line = (
                f'{row["Chromosome"]}\t{row["Source"]}\t{row["Feature"]}\t{row["Start"]}\t{row["End"]}\t{row["Score"]}\t{row["Strand"]}\t{row["Frame"]}\t{attributes}\n'
            )
            f.write(gtf_line)


df_to_gtf(
    null_full_gtf,
    "/oak/stanford/groups/smontgom/klawren/nph_gtex/references/marmoset_null_intergenic.gtf",
)