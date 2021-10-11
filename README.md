## MBG

<strong>M</strong>inimizer based sparse de <strong>B</strong>ruijn <strong>G</strong>raph constructor. Homopolymer compress input sequences, pick syncmers from hpc-compressed sequences, connect syncmers with an edge if they are adjacent in a read, unitigify and homopolymer decompress. Suggested input is PacBio HiFi/CCS reads. May or may not work with Illumina reads. Algorithmic details and citation: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab004/6104877

#### Installation

Bioconda: `conda install -c bioconda mbg`

#### Compilation

- `git clone https://github.com/maickrau/MBG.git`
- `cd MBG`
- `git submodule update --init --recursive`
- `make bin/MBG`

#### Usage

`MBG -i input_reads.fa -o output_graph.gfa -k kmer_size -w window_size -a kmer_min_abundance -u unitig_min_abundance`

eg `MBG -i reads.fa -o graph.gfa -k 1501 -w 1450 -a 1 -u 3`

Multiple read files can be inputted with "-i file1.fa -i file2.fa" etc. Input read type can be .fa / .fq / .fa.gz / .fq.gz.

#### Parameters

- `-k`: k-mer size. Must be odd and at least 31
- `-w`: window size. Cannot be greater than k-30. Default k-30
- `-a`: minimum k-mer abundance. Discard k-mers whose coverage is less than this. Default 1
- `-u`: minimum unitig abundance. Discard unitigs whose average coverage is less than this. Default 2
- `-t`: number of threads. Default 1
- `-i`: input read files. Can use multiple times to input multiple files. Format .fa/.fasta/.fa.gz/.fasta.gz/.fq/.fastq/.fq.gz/.fastq.gz
- `-o`: output graph file

Other options:
- `-r`: assemble using a multiplex DBG and increase k-mer size up to `r`. See [Bankevich et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.10.420448v2) for a description of the multiplex DBG algorithm.
- `-h`: print help
- `-v`: print version
- `--blunt`: output a graph without edge overlaps. Cannot be combined with `-r` or `--output-sequence-paths`
- `--error-masking`: mask errors in the input reads. Options are `hpc` (default) to mask homopolymer errors and call consensus on homopolymer run lengths, `no` to disable error masking, `collapse` to collapse homopolymer runs to one base pair (not recommended outside testing), `collapse-dinuc` to collapse homopolymer runs and call consensus on dinucleotide runs, `collapse-msat` to collapse homopolymer runs and call consensus on microsatellite where the repeat motif is up to 6bp long, `dinuc` to call consensus on dinucleotide runs, `msat` to call consensus on microsatellites where the repeat motif is up to 6bp long.
- `--include-end-kmers`: force k-mers at the ends of input sequences to be picked. Recommended if building from a reference of a linear genome, not recommended when building from reads or from a reference of a circular genome. Uses significantly more time and memory but prevents the last up to `w` base pairs at chromosome ends from being clipped.
- `--output-sequence-paths`: output the paths of the input sequences as alignments in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) to the given file. The alignments are exact in homopolymer compressed space but might differ in homopolymer run lengths.

k and w can be arbitrarily large but at some point the error rate and limited read length will cause the graph to be fragmented. Runtime stays approximately the same if the ratio k/w is kept constant. All repeats shorter than k are separated, all repeats longer than k+w are collapsed, and repeats in between may be separated or collapsed depending on if a k-mer was selected from within the repeat. When using `--blunt`, you should clean the graph afterwards with [vg](https://github.com/vgteam/vg). `--blunt` uses an extension of an algorithm invented by Hassan Nikaein (personal communication).
