# cleave2text
In silico protein digests in *.txt format

## About
cleave2text is an R/Bioconductor script that performs in silico digestion of a FASTA protein sequence and exports the resulting peptide sequences (as a non-redundant list) in text format. Optionally, the protein sequence can be truncated before it is digested, and peptides can be filtered to remove those that are shorter than a user-defined minimum length (the default behaviour is to remove single amino acids). For further details on usage, please see the header of the script file (**cleave2text.R**).

The script uses CRAN package [seqinR](https://cran.r-project.org/web/packages/seqinr/index.html) for the handling of protein sequence data, and Bioconductor package [cleaver](https://bioconductor.org/packages/release/bioc/html/cleaver.html) for the in silico digestion.

## Getting started
Section under construction..
