# Lim et al. (2021) processing code
## Summary
This directory contains data-munging code to combine the DNA, total RNA, and
polysome counts with corresponding 5' UTR sequence data from [Lim et al.
(2021)](https://dx.doi.org/10.1038/s41467-021-24445-6) ("Multiplexed
functional genomic analysis of 5’ untranslated region mutations across the
spectrum of prostate cancer").

## Input
The code takes three inputs:

1. `lim_2021_fig6a_utr5_seqs.csv`: lightly munged CSV export of [Supplementary
   Data
   6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-24445-6/MediaObjects/41467_2021_24445_MOESM8_ESM.xlsx)
   from the paper. After the CSV export, I removed a few extraneous rows toward
   the end and filtered out invalid rows with this command:

    grep -v -i -e "Incomplete set" -e "Not in pool" -e "Not in pacbio" | grep -v '^,,,$'

2. `lim_2021_totalrna_dna.csv`: total RNA and DNA counts provided via e-mail by
   the authors on 2022.04.11, then exported to CSV.

3. `lim_2021_totalrna_polysome.csv`: total RNA and polysome counts. For a given
   barcode, the total RNA counts in this file should exactly match those in
   `lim_2021_totalrna_dna.csv` (though this assumption was not validated beyond
   spot-checking a few cases by eye).

Note that `lim_2021_totalrna_dna.csv` and `lim_2021_totalrna_polysome.csv` were
extracted from a single Excel file, where they were present in the same sheet
as distinct column sets.

## Outputs
The code produces two output files.

1. `mutants.csv`: mutated 5' UTR sequences where, for a given barcode, all of
   (polysome, total RNA, DNA) read-outs were available. Only SNVs should be
   present.

2. `wildtype.csv`: 5' UTR wildtype sequences where, for a given barcode, all of
   (polysome, total RNA, DNA) read-outs were available

For both files, the full 5' UTR sequences for the genes (with or without
SNV) should be provided. Records were retained only when the 5' UTR
sequence could be resolved unambiguously from the gene (and, for mutant
records, the mutation) listed.

Read-outs were combined (by taking the mean) across multiple barcodes (i.e.,
collapsing the dataframes row-wise), and across multiple replicates (i.e.,
collapsing the dataframes column-wise).

Exactly the same genes should be represented in `mutants.csv`, with a
many-to-one relationship between mutant and wildtype records (i.e., every
wildtype sequence should have one or more mutated sequences in correspondence).

Values for DNA, RNA, and polysome are in units of counts-per-million (CPM).
