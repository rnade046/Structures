import pandas as pd


def load_biomart_mapping(biomart_file):

    mapping = pd.read_csv(biomart_file, sep="\t")
    return mapping


def filter_auraDB(auradb_file, hgnc_symbols, fasta_file):

    db = pd.read_csv(auradb_file, sep="\t")

    filtered_db = db[
        db["Target Gene"].isin(hgnc_symbols) &
        db["Target UTR"].str.contains("3UTR", na=False)
    ]

    # get identifiers for fasta file (i.e. remove "hg19_" prefix of UTR ids)
    filtered_db["Target_UTR_id"] = filtered_db["Target UTR"].str.replace(r"^hg19_", "", regex=True)

    # get genome positions for utr in database
    utr_pos = get_utr_genome_positions_fasta(fasta_file, set(filtered_db["Target_UTR_id"].dropna().unique()))
    filtered_db = filtered_db.merge(utr_pos, on="Target_UTR_id", how="inner")

    # get position of binding site from chromosome pos [utr_start/end] +/- [binding site_start/end]
    minus_mask = (
            (filtered_db["strand"] == "-") &
            (filtered_db["Binding site start"] < filtered_db["Target UTR length"])
    )

    filtered_db.loc[minus_mask, "binding_start"] = (
            filtered_db.loc[minus_mask, "end_site"] -
            filtered_db.loc[minus_mask, "Binding site start"]
    )

    filtered_db.loc[minus_mask, "binding_end"] = (
            filtered_db.loc[minus_mask, "end_site"] -
            filtered_db.loc[minus_mask, "Binding site stop"]
    )

    positive_mask = (
            (filtered_db["strand"] == "+") &
            (filtered_db["Binding site start"] < filtered_db["Target UTR length"])
    )

    filtered_db.loc[positive_mask, "binding_start"] = (
            filtered_db.loc[positive_mask, "start_site"] +
            filtered_db.loc[positive_mask, "Binding site start"]
    )

    filtered_db.loc[positive_mask, "binding_end"] = (
            filtered_db.loc[positive_mask, "start_site"] +
            filtered_db.loc[positive_mask, "Binding site stop"]
    )

    # remove motifs that didn't match in the aura db
    filtered_db = filtered_db.dropna(subset=["binding_start", "binding_end"])

    return filtered_db


def get_utr_genome_positions_fasta(fasta_file, utr_ids):
    records = []
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            line = line.strip()

            # ignore sequence or empty lines
            if not line.startswith(">"):
                continue

            parts = line.split()
            utr_id = parts[0].lstrip(">")

            if utr_id in utr_ids:
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
                    "Target_UTR_id": utr_id,
                    "chr": chrom,
                    "start_site": int(start_site),
                    "end_site": int(end_site),
                    "strand": strand,
                })

    return pd.DataFrame(records)


def match_positions(bed_file, aura_db, out_file, clean):
    motif = pd.read_csv(bed_file, sep="\t", header=None)
    motif.columns = ["chr", "motif_start", "motif_end", "pos", "idx"]

    matches = motif.merge(aura_db, on="chr", how="inner")

    matches = matches[
        (matches["motif_start"] < matches["binding_end"]) &
        (matches["motif_end"] > matches["binding_start"])
        ]

    if not matches.empty:
        if not clean:
            matches.to_csv(out_file, sep="\t", index=False)
        else :
            matches2 = clean_matches(matches)
            matches2[["chr", "motif_start","motif_end",  "Regulatory factor",
                "Target Gene", "Target_UTR_ids",   "binding_start", "binding_end"]
            ].to_csv(out_file, sep="\t", index=False)


def clean_matches(match_df):
    summary_df = match_df[
        [
            "chr",
            "motif_start",
            "motif_end",
            "Regulatory factor",
            "Target Gene",
            "Target_UTR_id",
            "strand",
            "binding_start",
            "binding_end",
        ]
    ].copy()

    summary_df = (
        summary_df
        .groupby(
            ["motif_start", "motif_end", "binding_start", "binding_end"],
            as_index=False,
            dropna=False,
        )
        .agg({
            "chr": "first",
            "Regulatory factor": "first",
            "Target Gene": "first",
            "strand": "first",
            "Target_UTR_id": lambda x: ", ".join(sorted(set(x.dropna()))),
        })
        .rename(columns={"Target_UTR_id": "Target_UTR_ids"})
    )

    return summary_df

if __name__ == '__main__':

    wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/AuraDB"
    biomart_file = f"{wd}/input/BiomaRt_corrNet2_refSeq_UCSD.tsv"
    cis_auradb_file = f"{wd}/input/AURAlight_cis-elements.txt"
    trans_auradb_file = f"{wd}/input/AURAlight_trans-factors.txt"
    fasta_file = f"{wd}/input/UTR_hg19.fasta"


    # load cis / trans db : limit to 3'UTRs associated to target genes network (w/ BiomaRt file)
    id_mapping = load_biomart_mapping(biomart_file)
    unique_hgnc = set(id_mapping["hgnc_symbol"].dropna().unique())

    # obtain chromosome positions for target sites (hg19.fasta)
    cisdb = filter_auraDB(cis_auradb_file, unique_hgnc, fasta_file)
    transdb = filter_auraDB(trans_auradb_file, unique_hgnc, fasta_file)

    for i in [82, 85, 94, 239]:
        bed_file = f"{wd}/bed/hglft_genome_{i}.bed"
        cis_out_file = f"{wd}/output/cisdb_motif{i}.tsv"
        trans_out_file = f"{wd}/output/trans_motif{i}.tsv"

        # check for matching chromosome regions
        match_positions(bed_file, cisdb, cis_out_file, False)
        match_positions(bed_file, transdb, trans_out_file, True)
