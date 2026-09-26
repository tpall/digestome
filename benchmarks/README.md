# Benchmark inputs

`gtdb_r226_sample.tsv` is the 3,077-genome GTDB r226 species-representative sample behind the GTDB
benchmark in the manuscript (3,043 of them yielded a proteome). It is the output of
`scripts/sample_gtdb_reps.py` as first released (commit `48e9bae`) run on GTDB's
`{ar53,bac120}_taxonomy_r226.tsv` restricted to the representatives in `sp_clusters_r226.tsv`.
The current sampler also counts *Methanonatronarchaeales* as a methanogen order, which moves five
background genomes, so this list is the fixed reference rather than the sampler's current output.

Columns: accession, phylum, order, GTDB lineage, stratum (`methanogen` or `background`, before the
two sourced corrections described in the manuscript's Methods).
