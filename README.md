# HerbChomper
A bioinformatic tool for trimming poorly-aligned ends from DNA sequences.

DNA sequences assembled from herbarium and other museum specimens often have error-rich ends that do not align well. Poorly-aligned ends on otherwise reliable sequences can lead to unnecessary data loss when filters are applied only to whole sequences or whole columns in an alignment. HerbChomper trims the ends of a target sequence until identity to a reference within the same alignment calculated within a sliding window reaches a user-defined threshold. 

NEW UPDATE: June 7, 2020 - Beta version 0.3 
Automated reference choosing is much faster in this version.

Installation: download herbchomper.R or clone this repository.

Dependencies: R and the SeqinR package. The latter will be automatically installed if it is not already present.

Usage: Rscript herbchomper.R -a [alignment in] -o [alignment out] -t [sequence to trim] -r [reference sequence or "auto"] -w [size of sliding window] -i [identity cutoff] -g [gap size necessary to restart trimming]

The reference can be user specified or chosen automatically by specifying "-a auto". In the latter case, the sequence with the highest similarity to the target will be chosen, ignoring gaps and undetermined characters.

Example command line using the provided data (exampleData.zip):
Rscript herbchomper.R -a gene038.single.trimmed -o gene038.single.chomped -t A_jarrettiae_SAN120933 -r auto -w 10 -i 0.8 -g 20
