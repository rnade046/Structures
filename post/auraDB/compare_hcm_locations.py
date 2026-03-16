import pandas as pd

def load_hcm_locations(hcm_file, factors):
    fields = ['symbol', 'MMF localization', 'SAFE domain']
    atlas = pd.read_csv(hcm_file, sep="\t", usecols=fields)

    # keep atlas rows for RBP or target proteins
    genes_to_keep = set(factors["Regulatory factor"]).union(factors["Target Gene"])
    atlas = atlas[atlas["symbol"].isin(genes_to_keep)].copy()

    return atlas


def load_factors(trans_file):
    fields = ['Regulatory factor', 'Target Gene']
    factors = pd.read_csv(trans_file, sep="\t", usecols=fields).drop_duplicates()
    return factors


def assess_shared_nmf_locations(factors, atlas, out_file):
    location_map = atlas.set_index("symbol")["MMF localization"]
    out = factors.copy()
    out["rf_locations"] = out["Regulatory factor"].map(location_map)
    out["tg_locations"] = out["Target Gene"].map(location_map)

    out["matching_locations"] = out.apply(get_overlap, axis=1)
    out.to_csv(out_file, sep="\t", index=False)


def assess_shared_safe_locations(factors, atlas, out_file):
    location_map = atlas.set_index("symbol")["MMF localization"]
    out = factors.copy()
    out["rf_locations"] = out["Regulatory factor"].map(location_map)
    out["tg_locations"] = out["Target Gene"].map(location_map)

    out["matching_locations"] = out.apply(get_overlap, axis=1)
    out.to_csv(out_file, sep="\t", index=False)

def get_overlap(row):
    if not isinstance(row["rf_locations"], set) or not isinstance(row["tg_locations"], set):
        return pd.NA

    overlap = row["rf_locations"] & row["tg_locations"]
    return overlap if overlap else pd.NA


if __name__ == '__main__':
    wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/AuraDB"

    hcm_location_file = f"{wd}/hcm/preys-latest.txt"

    for m in [82, 85, 94, 239]:
        motif = f"{wd}/output/trans_motif{m}.tsv"
        nmf_out_file = f"{wd}/hcm/hcm_NMF_location_auradb_transfactors_motif{m}.tsv"
        safe_out_file = f"{wd}/hcm/hcm_SAFE_location_auradb_transfactors_motif{m}.tsv"

        transfactors = load_factors(motif)
        locations = load_hcm_locations(hcm_location_file, transfactors)
        assess_shared_nmf_locations(transfactors, locations, nmf_out_file)
        assess_shared_safe_locations(transfactors, locations, safe_out_file)