## MBG

<strong>M</strong>inimizer based sparse de <strong>B</strong>ruijn <strong>G</strong>raph constructor. Homopolymer compress input sequences, winnow minimizers from hpc-compressed sequences, connect minimizers with an edge if they are adjacent in a read, unitigify. Suggested input is PacBio HiFi/CCS reads. May or may not work with Illumina reads. Algorithmic details: https://www.biorxiv.org/content/10.1101/2020.09.18.303156v1

#### Installation

Bioconda: `conda install -c bioconda mbg`

#### Compilation

- `git clone https://github.com/maickrau/MBG.git`
- `cd MBG`
- `git submodule update --init --recursive`
- `make bin/MBG`

#### Usage

`MBG -i input_reads.fa -o output_graph.gfa -k kmer_size -w window_size -a kmer_min_abundance -u unitig_min_abundance`

eg `MBG -i reads.fa -o graph.gfa -k 1501 -w 1500 -a 1 -u 3`

Multiple read files can be inputted with "-i file1.fa -i file2.fa" etc. Input read type can be .fa / .fq / .fa.gz / .fq.gz.

#### Parameters

- `-k`: k-mer size for minimizer winnowing
- `-w`: window size for minimizer winnowing. Cannot be greater than k
- `-a`: minimum k-mer abundance. Discard minimizers whose coverage is less than this
- `-u`: minimum unitig abundance. Discard unitigs whose average coverage is less than this, and discard edges whose coverage is less than this

Other options:
- `-i`: input read files. Can use multiple times to input multiple files
- `-o`: output graph file
- `-h`: print help
- `-v`: print version
- `--no-hpc`: don't homopolymer compress the reads. Not recommended outside testing purposes
- `--collapse-hpc`: disable homopolymer run length consensus and collapse runs to one bp. Recommended if the reads are already homopolymer compressed, otherwise not

k and w can be arbitrarily large but at some point the error rate and limited read length will cause the graph to be fragmented. Runtime stays approximately the same if the ratio k/w is kept constant. All repeats shorter than k are separated, all repeats longer than k+w are collapsed, and repeats in between may be separated or collapsed depending on if a minimizer was selected from within the repeat.
