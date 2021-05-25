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
- `-w`: window size. Cannot be greater than k-30
- `-a`: minimum k-mer abundance. Discard k-mers whose coverage is less than this
- `-u`: minimum unitig abundance. Discard unitigs whose average coverage is less than this
- `-t`: number of threads

Other options:
- `-i`: input read files. Can use multiple times to input multiple files
- `-o`: output graph file
- `-h`: print help
- `-v`: print version
- `--blunt`: output a graph without edge overlaps
- `--no-hpc`: don't homopolymer compress the reads. Not recommended outside testing purposes
- `--collapse-hpc`: disable homopolymer run length consensus and collapse runs to one bp. Recommended if the reads are already homopolymer compressed, otherwise not
- `--include-end-kmers`: force k-mers at the ends of input sequences to be picked. Recommended if building from a reference of a linear genome, not recommended when building from reads or from a reference of a circular genome. Uses significantly more time and memory but prevents the last up to `w` base pairs at chromosome ends from being clipped.

k and w can be arbitrarily large but at some point the error rate and limited read length will cause the graph to be fragmented. Runtime stays approximately the same if the ratio k/w is kept constant. All repeats shorter than k are separated, all repeats longer than k+w are collapsed, and repeats in between may be separated or collapsed depending on if a k-mer was selected from within the repeat. When using `--blunt`, you should clean the graph afterwards with [vg](https://github.com/vgteam/vg). `--blunt` uses an extension of an algorithm invented by Hassan Nikaein (personal communication).
