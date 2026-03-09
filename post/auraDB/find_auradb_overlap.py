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

    coords = (
        filtered_db["Target_UTR_id"]
        .map(utr_pos)
        .str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    )
    coords.columns = ["chr", "utr_start", "utr_end"]

    filtered_db = pd.concat([filtered_db, coords], axis=1)
    filtered_db["utr_start"] = filtered_db["utr_start"].astype("Int64")
    filtered_db["utr_end"] = filtered_db["utr_end"].astype("Int64")

    # get position of binding site from chromosome pos [utr_start/end] +/- [binding site_start/end]

    return filtered_db


def get_utr_genome_positions_fasta(fasta_file, utr_ids):
    utr_to_range = {}
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            line = line.strip()

            # ignore sequence or empty lines
            if not line.startswith(">"):
                continue

            parts = line.split()
            utr_id = parts[0].lstrip(">")
            chrom_range = parts[1].split("=", 1)[1]

            if utr_id in utr_ids:
                utr_to_range[utr_id] = chrom_range
    return utr_to_range


if __name__ == '__main__':

    wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/AuraDB"
    biomart_file = f"{wd}/input/BiomaRt_corrNet2_refSeq_UCSD.tsv"
    cis_auradb_file = f"{wd}/input/AURAlight_cis-elements.txt"
    fasta_file = f"{wd}/input/UTR_hg19.fasta"

    # load cis / trans db : limit to 3'UTRs associated to target genes network (w/ BiomaRt file)
    id_mapping = load_biomart_mapping(biomart_file)
    unique_hgnc = set(id_mapping["hgnc_symbol"].dropna().unique())

    # obtain chromosome positions for target sites (hg19.fasta)
    cisdb = filter_auraDB(cis_auradb_file, unique_hgnc, fasta_file)

    # check for matching chromosome regions
