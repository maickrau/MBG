## MBG

Minimizer based de Bruijn-like Graph constructor. Homopolymer compress input sequences, winnow minimizers from hpc-compressed sequences, connect minimizers with an edge if they are adjacent in a read, unitigify. Suggested input is PacBio HiFi/CCS reads. May or may not work with Illumina reads.

#### Compilation

- `git clone https://github.com/maickrau/MBG.git`
- `cd MBG`
- `git submodule update --init --recursive`
- `make bin/MBG`

#### Usage

`bin/MBG -i input_reads.fa -o output_graph.gfa -k k -w w -a kmerabundance -u unitigabundance`

eg `bin/MBG -i reads.fa -g graph.gfa -k 2501 -w 2000 -a 1 -u 3`

#### Parameters

- k: k-mer size for minimizer winnowing
- w: window size for minimizer winnowing. Cannot be greater than k
- kmerabundance: discard minimizers whose coverage is less than this
- unitigabundance: discard unitigs whose average coverage is less than this, and discard edges whose coverage is less than this

k and w can be arbitrarily large but at some point the error rate and limited read length will cause the graph to be fragmented. Runtime stays approximately the same if the ratio k/w is kept constant. All repeats shorter than k are separated, all repeats longer than k+w are collapsed, and repeats in between may be separated or collapsed depending on if a minimizer was selected from within the repeat.

#### Runtime

Approximately O(s k / w + s / w) where s is amount of sequence in base pairs, k is k-mer size and w is window size. Runtime is dominated by reading & winnowing sequences and the graph construction is trivial.
