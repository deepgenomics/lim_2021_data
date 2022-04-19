import pandas as pd
import re
from pathlib import Path


def _get_data_fn(fn):
    return Path(__file__).parent / "data" / fn


def _process_counts():
    df = [
        pd.read_csv(_get_data_fn(fn))
        for fn in ("lim_2021_totalrna_dna.csv", "lim_2021_totalrna_polysome.csv")
    ]

    combined = pd.merge(df[0], df[1], on=("description", "barcode"), validate="1:1")
    overlapping = [col for col in combined.columns if col.endswith("_x")]
    for col in overlapping:
        base = col[:-2]
        other = base + "_y"
        assert other in combined.columns
        assert (combined[col] == combined[other]).all()
        combined.drop(other, axis=1, inplace=True)
        combined.rename({col: base}, axis=1, inplace=True)

    mutant_counts = {}
    wt_counts = {}
    bad_genes = set()

    for desc, group in combined.groupby("description"):
        match = re.search(
            r"^([A-Za-z0-9]+)_([ACGT])_([ACGT])_chr[0-9XY]+_\d+_\d+$", desc
        )
        if match:
            gene, ref, mut = match.groups()
            # We should only expect a single set of counts to be provided for a
            # given (gene, ref, mut) tuple. This assumption can be violated if
            # the same type of mutation (e.g., G -> A) occurs at multiple
            # positions in the UTR sequence. Determining the meaning of these
            # positions, and/or matching them to the positions given in the
            # sequences dataframe, has been difficult, so we ignore them.
            # Consequently, when this assumption is violated, discard the gene
            # and all its mutations to ensure we're not providing incorrectly
            # matched data.
            if (gene, ref, mut) in mutant_counts:
                bad_genes.add(gene)
            else:
                mutant_counts[(gene, ref, mut)] = group
            continue

        # We can match multiple line endings with this. Sometimes a position is
        # provided with the wildtype (after the chromosome), but other times it
        # isn't.
        match = re.search(r"^([A-Za-z0-9]+)_WT_chr[0-9XY]", desc)
        if match:
            gene = match.group(1)
            # Again, expect only one wildtype sequenec to be provided for a gene.
            if gene in wt_counts:
                bad_genes.add(gene)
            else:
                wt_counts[gene] = group
            continue

        # print("No match for {}".format(desc))

    # Every mutant should have a corresponding wildtype, and each wildtype
    # should have one or more mutants.
    good = (
        set(wt_counts.keys()) & set([k[0] for k in mutant_counts.keys()])
    ) - bad_genes
    mutant_counts = {K: V for (K, V) in mutant_counts.items() if K[0] in good}
    wt_counts = {K: V for (K, V) in wt_counts.items() if K in good}
    # print(len(bad_genes), len(wt_counts), len(mutant_counts))
    # print(
    #    len(combined),
    #    sum(len(g) for g in mutant_counts.values()),
    #    sum(len(g) for g in wt_counts.values()),
    # )
    return (wt_counts, mutant_counts)


def _process_seqs():
    df = pd.read_csv(_get_data_fn("lim_2021_fig6a_utr5_seqs.csv"))

    wt_seqs = {}
    mutant_seqs = {}
    bad_genes = set()

    for idx in range(len(df)):
        name = df.iloc[idx]["5' UTR genomic coordinate"]
        seq = df.iloc[idx]["sequence of 5' UTR"]

        # Handle SNVs
        match = re.search(
            r"^([A-Za-z0-9]+)_chr[0-9XY]+_\d+_([ACGT])_([ACGT])_UTR5$", name
        )
        if match:
            gene, ref, mut = match.groups()
            # As with the count data, we expect only one sequence to be
            # provided for a given (gene, ref, mut). This assumption can be
            # violated if the same mutation occurs at different positions of
            # the 5' UTR, but we can't make sense of the genomic coordinates
            # that are provided, so we throw out all mutations when this case
            # occurs.
            if (gene, ref, mut) in mutant_seqs:
                bad_genes.add(gene)
            else:
                mutant_seqs[(gene, ref, mut)] = seq
            continue

        # Handle deletions (though we don't actually use these, I think, since
        # no count data seem to be provided for these sequences)
        match = re.search(r"^([A-Za-z0-9]+)_chr[0-9XY]+_\d+_del([ACGT]+)_UTR5$", name)
        if match:
            gene, deletion = match.groups()
            ref, mut = deletion, ""
            if (gene, ref, mut) in mutant_seqs:
                bad_genes.add(gene)
            else:
                mutant_seqs[(gene, ref, mut)] = seq
            continue

        # Handle wildtype
        match = re.search(r"^([A-Za-z0-9]+)_chr[0-9XY]+_\d+_WT_UTR5$", name)
        if match:
            gene = match.group(1)
            if gene in wt_seqs:
                bad_genes.add(gene)
            else:
                wt_seqs[gene] = seq
            continue

        raise Exception("No match for {}".format(name))

    # Every wildtype seq should have one or more mutant seqs, and every mutant
    # seq should have exactly one wildtype seq.
    good = (set(wt_seqs.keys()) & set([k[0] for k in mutant_seqs.keys()])) - bad_genes
    mutant_seqs = {K: V for (K, V) in mutant_seqs.items() if K[0] in good}
    wt_seqs = {K: V for (K, V) in wt_seqs.items() if K in good}
    # print(len(bad_genes), len(wt_seqs), len(mutant_seqs))
    return (wt_seqs, mutant_seqs)


def _retain_only_intersection(A, B):
    common = set(A.keys()) & set(B.keys())
    A = {K: V for (K, V) in A.items() if K in common}
    B = {K: V for (K, V) in B.items() if K in common}
    return (A, B)


def _add_seqs(wt_seqs, mutant_seqs, wt_counts, mutant_counts):
    mutants = pd.DataFrame()
    for (gene, ref, mut), group in mutant_counts.items():
        group["gene"] = gene
        group["ref"] = ref
        group["mut"] = mut
        group["seq"] = mutant_seqs[(gene, ref, mut)]
        # Reset the index via `ignore_index` so that it takes continguous
        # values.
        mutants = pd.concat((mutants, group), ignore_index=True)

    wt = pd.DataFrame()
    for gene, group in wt_counts.items():
        group["gene"] = gene
        group["seq"] = wt_seqs[gene]
        wt = pd.concat((wt, group), ignore_index=True)

    return (wt, mutants)


def _process():
    wt_seqs, mutant_seqs = _process_seqs()
    wt_counts, mutant_counts = _process_counts()

    # Require that we have sequence and count data for every entity, either mutant or wildtype.
    wt_seqs, wt_counts = _retain_only_intersection(wt_seqs, wt_counts)
    mutant_seqs, mutant_counts = _retain_only_intersection(mutant_seqs, mutant_counts)

    # Generate new dataframes
    return _add_seqs(wt_seqs, mutant_seqs, wt_counts, mutant_counts)


def _aggregate_measurements(df):
    """
    Aggregate replicates (across columns) and barcodes (across rows) for measurements of same sequence.
    """
    assert df.isna().any(axis=1).sum() == 0

    rep_group_cols = set()
    rep_cols = set(col for col in df.columns if re.search("_rep\d+$", col))
    while len(rep_cols) > 0:
        col = rep_cols.pop()
        prefix = re.search("^(.+)_rep\d+$", col).group(1)
        rep_group_cols.add(prefix)
        measure_cols = set((col,)) | set(c for c in rep_cols if c.startswith(prefix))
        assert col in measure_cols
        df[prefix] = df[list(measure_cols)].mean(axis=1)
        rep_cols -= measure_cols
        df.drop(measure_cols, axis=1, inplace=True)

    df.drop("barcode", axis=1, inplace=True)
    other_cols = list(set(df.columns) - rep_group_cols)
    df = df.groupby(other_cols, as_index=False)[list(rep_group_cols)].mean()
    return df


def main():
    outputs = {}
    outputs["wildtype"], outputs["mutants"] = _process()
    for name, df in outputs.items():
        df = _aggregate_measurements(df)
        df.to_csv(_get_data_fn(name + ".csv"), index=False)


main()
