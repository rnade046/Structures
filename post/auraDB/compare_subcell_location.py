import pandas as pd

def load_atlas(atlas_file, factors):
    fields = ['Gene name', 'Reliability', 'Main location', 'Additional location']
    atlas = pd.read_csv(atlas_file, sep="\t", usecols=fields)

    # remove rows with Reliability='uncertain'
    atlas = atlas[atlas['Reliability'] != 'Uncertain'].copy()

    # keep atlas rows for RBP or target proteins
    genes_to_keep = set(factors["Regulatory factor"]).union(factors["Target Gene"])
    atlas = atlas[atlas["Gene name"].isin(genes_to_keep)].copy()

    # parse locations
    atlas["all_locations"] = atlas.apply(parse_locations, axis=1)
    atlas = atlas.drop(columns=["Main location", "Additional location"])
    return atlas


def parse_locations(row):
    main = [] if pd.isna(row["Main location"]) else row["Main location"].split(";")
    add = [] if pd.isna(row["Additional location"]) else row["Additional location"].split(";")

    # combine main and additional locations
    return {
        loc.strip()
        for loc in main + add
        if loc and loc.strip()
    }

def load_factors(trans_file):
    #fields = ['Regulatory factor', 'Target Gene']
    factors = pd.read_csv(trans_file, sep="\t").drop_duplicates()
    return factors


def assess_shared_locations(factors, atlas, out_file, summary_out_file):
    location_map = atlas.set_index("Gene name")["all_locations"]
    out = factors.copy()
    out["rf_locations"] = out["Regulatory factor"].map(location_map)
    out["tg_locations"] = out["Target Gene"].map(location_map)

    out["matching_locations"] = out.apply(get_overlap, axis=1)
    out.to_csv(out_file, sep="\t", index=False)

    out["co_localized"] = out["matching_locations"].notna()

    cols_to_keep = ["chr", "motif_start", "motif_end", "Regulatory factor", "Target Gene",
                    "Target_UTR_ids", "binding_start", "binding_end", "co_localized"]

    out[cols_to_keep].to_csv(summary_out_file, sep="\t", index=False)


def get_overlap(row):
    if not isinstance(row["rf_locations"], set) or not isinstance(row["tg_locations"], set):
        return pd.NA

    overlap = row["rf_locations"] & row["tg_locations"]
    return overlap if overlap else pd.NA

if __name__ == '__main__':
    wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/AuraDB"

    atlas_file = f"{wd}/protein_atlas/subcellular_location.tsv"

    for m in [82, 85, 94, 239]:
        motif = f"{wd}/output/trans_motif{m}.tsv"
        matched_location_file = f"{wd}/protein_atlas/atlas_subcell_location_auradb_transfactors_motif{m}.tsv"
        summary_file = f"{wd}/output_summary/auradb_transfactors_location_motif{m}.tsv"

        trans = load_factors(motif)
        atlas = load_atlas(atlas_file, trans)
        assess_shared_locations(trans, atlas, matched_location_file, summary_file)

