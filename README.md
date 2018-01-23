# idinf
Tools for disambiguating metagenomic reads

**phylowalk_m8parser.py**

Implements taxonomic tree reversal for .m8 files.

`python3 phylowalk_m8parser.py -i [input .m8 file] -o [output file]`


**align_and_convert_reference_microbe.sh**

For any input .fasta file, runs GSNAPL allowing 1000 matches and subsequently launches the phylowalk_m8parser.py script

`bash align_and_conver_reference_microbe.sh [fasta file name, not including .fasta]`

ie: `bash align_and_convert_reference_microbe.sh UGANDA/UGANDA_S9_NP_orthopneumovirus-hits`

