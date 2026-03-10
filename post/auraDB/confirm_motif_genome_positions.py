import pandas as pd


def load_biomart_mapping(biomart_file):

    mapping = pd.read_csv(biomart_file, sep="\t")
    return mapping


def load_utr_positions(fasta_file, refseq_ids):
    records = []
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            line = line.strip()

            # ignore sequence or empty lines
            if not line.startswith(">"):
                continue

            parts = line.split()
            utr_id = parts[0].lstrip(">").removeprefix("hg38_ncbiRefSeq_").split(".", 1)[0]

            if utr_id in refseq_ids:
                chrom_range = None
                strand = None

                for field in parts[1:]:
                    if "=" not in field:
                        continue

                    key, value = field.split("=", 1)

                    if key == "range":
                        chrom_range = value
                    elif key == "strand":
                        strand = value
                chrom, pos = chrom_range.split(":", 1)
                start_site, end_site = pos.split("-", 1)

                records.append({
                    "utr_id": utr_id,
                    "chr": chrom,
                    "start_site": int(start_site),
                    "end_site": int(end_site),
                    "strand": strand,
                })
    return pd.DataFrame(records)


def load_bed_file(bed_file):
    bed = pd.read_csv(bed_file, sep="\t", skiprows = 1, header = None)
    bed.columns = ["chr", "motif_start", "motif_end"]
    bed["motif_id"] = range(1, len(bed) + 1)
    return bed


def compare_positions(utr_positions, motif_positions, summary_out_file, bed_out_file):
    matches = motif_positions.merge(utr_positions, on="chr", how="inner")

    matches = matches[
        (matches["motif_start"] < matches["end_site"]) &
        (matches["motif_end"] > matches["start_site"])
        ]

    all_plus_mask = matches.groupby("motif_id")["strand"].transform(lambda s: (s == "+").all())

    plus_only = matches[all_plus_mask].drop_duplicates(subset=["motif_id", "motif_start", "motif_end"])
    other = matches[~all_plus_mask].drop_duplicates(subset=["motif_id", "start_site", "end_site"])

    clean_matches = pd.concat([plus_only, other], ignore_index=True)
    clean_matches.to_csv(summary_out_file, index=False)

    minus_mask = clean_matches["strand"] == "-"

    clean_matches.loc[minus_mask, "relative_pos"] = (
            clean_matches.loc[minus_mask, "motif_start"] -
            clean_matches.loc[minus_mask, "start_site"]
    )

    new_start = (
            clean_matches.loc[minus_mask, "end_site"] -
            clean_matches.loc[minus_mask, "relative_pos"]
    )

    clean_matches.loc[minus_mask, "motif_start"] = new_start
    clean_matches.loc[minus_mask, "motif_end"] = new_start + 1
    clean_matches.to_csv(summary_out_file, index=False)

    clean_matches[["chr", "motif_start", "motif_end"]].to_csv(bed_out_file, sep="\t", index=False, header=False)

    return clean_matches

if __name__ == '__main__':

    wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/AuraDB"
    biomart_file = f"{wd}/input/BiomaRt_corrNet2_refSeq_UCSD.tsv"
    fasta_file = f"{wd}/input/human_3UTRsequences.txt"


    # limit to 3'UTRs associated to target genes network (w/ BiomaRt file)
    id_mapping = load_biomart_mapping(biomart_file)
    unique_refseq = set(id_mapping["refseq_mrna"].dropna().unique())

    utr_positions = load_utr_positions(fasta_file, unique_refseq)

    for m in [82,85,94,239]:
        bed_file_prefix = f"{wd}/input/positions/genomePositions_module{m}.bed"
        motif_positions_file = f"{wd}/bed/motif{m}_positions_corrected.csv"
        bed_out_file = f"{wd}/bed/motif{m}_corrected.bed"

        motif_positions = load_bed_file(bed_file_prefix)

        matches = compare_positions(utr_positions, motif_positions, motif_positions_file, bed_out_file)