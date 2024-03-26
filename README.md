# cleave2text
In silico protein digests in *.txt format

## About
cleave2text is an R/Bioconductor script that performs in silico digestion of a FASTA protein sequence and exports the resulting peptide sequences (as a non-redundant list) in text format. Optionally, the protein sequence can be truncated* and/or exported before it is digested, and peptides can be filtered to remove those that are shorter than a user-defined minimum length (the default behaviour is to remove single amino acids). For further details on usage, please see the header of the script file (**cleave2text.R**).

The script uses CRAN package [seqinR](https://cran.r-project.org/web/packages/seqinr/index.html) for the handling of protein sequence data, and Bioconductor package [cleaver](https://bioconductor.org/packages/release/bioc/html/cleaver.html) for the in silico digestion. The software and package versions under which the script was developed and tested can be found in **v1.0.0_testrun_sessionInfo.txt**.

\* I would avoid using the truncation function on anything other than original full-length sequences from UniProt, otherwise things could get very confusing. Specifically, you could end up with amino-acid index numbers that no longer relate to the UniProt accession number in your FASTA headers/file names. 

## Getting started
1. Download the script file (**cleave2text.R**) and the example input file (**uniprotkb_2024_02_19_test-set.fasta**).
2. Open R and, if necessary, install the required packages.
3. Transfer the example input file to your R working directory [if you don't know where this is, use the `getwd()` command to find out].
4. Open the script file in R (File > Open script).
5. Run the script (Edit > Run all).
6. Review the output in the R console, noting any warnings.
7. Check your working directory for output files. You should find a text file containing an in silico digest of the SARS-CoV-2 spike glycoprotein.
8. Try re-running the script with different user-defined parameters (see script header for details).

## Acknowledgements
The sequences in the example input file were obtained from [UniProtKB](https://www.uniprot.org/).
