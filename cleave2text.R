# +-----------------------------------------------------------------------------------+
#  cleave2text v1.0.0
# +-----------------------------------------------------------------------------------+

# Notes:
# 1.  This script performs in silico digestion of a protein sequence and filters the
#      resulting peptides, removing those that are shorter than a user-defined minimum
#      length. Peptides passing through the filter are exported in one of two possible
#      text formats.
# 2.  The script requires a *.fasta file containing one or more protein sequences, 
#      plus the index number of the sequence to be processed (1 = first or only 
#      sequence, 2 = second sequence, and so on). FASTA headers should be UniProtKB-
#      style (i.e., of the form
#                              >databaseID|proteinID|description
#      where 'databaseID', 'proteinID' and 'description' are variable strings).
# 3.  Optionally, the sequence can be truncated prior to digestion (e.g., to remove a 
#      signal peptide). Truncation is achieved by defining the part of the sequence
#      that you want to keep (and be digested) rather than the part that you want to
#      remove.
# 4.  Optionally, the sequence - with any truncation applied - can be saved in FASTA
#      format before it is digested.
# 5.  CRAN package 'seqinr' is required for the handling of protein sequences. 
#      Bioconductor package 'cleaver' is required for the in silico digestion.
# 6.  The following user-defined parameters are set within the body of the script:
#       fasta.file.path     The path to a *.fasta file.
#       fasta.element       The index number of the sequence to be digested (default:
#                            1).
#       residues            Either 'all' (the default) or a range (e.g., '2-375').
#       save.residues       Logical, indicating whether or not the segment defined by
#                            'residues' should be exported as a new *.fasta file 
#                            (default: FALSE).
#       rule                The cleavage rule, as defined in the 'cleaver' manual
#                            (default: 'trypsin').
#       missed              A set of values for the 'missedCleavages' argument. 
#                            Examples:
#                             -------------------------------------------------------
#                              Value(s)            Products
#                             -------------------------------------------------------
#                              0 (= the default)   Peptides with 0 missed cleavages
#                              1                   Peptides with 1 and only 1 missed
#                                                   cleavage (i.e., none with 0 or >1
#                                                   missed cleavages).
#                              c(0, 1, 2)          Peptides with <=2 missed
#                                                   cleavages.
#                             --------------------------------------------------------
#
#       min.peptide.length  The minimum length of peptide to be included in the output 
#                            files. The default value of 2 excludes only single amino 
#                            acids.
#       output.folder.path  The path to a folder. The default is the current working 
#                            directory path, which is stored in object 'wd'. Must 
#                            include a trailing path separator ('/' or '\\').
#       collapse.isd        Logical, enabling output format to be switched as follows:
#                             ---------------------------------------
#                              Value      Format
#                             ---------------------------------------
#                              TRUE       PEPTIDE1 PEPTIDE2 PEPTIDE3
# 
#                              FALSE      PEPTIDE1
#                              (= the     PEPTIDE2
#                              default)   PEPTIDE3
#                             ---------------------------------------
#
# 7.  Specifically, 'collapse.isd' controls whether or not the in silico digest 
#      ('isd') should be collapsed into a single string of space-separated values 
#      before it is written to file.
# 8.  A final parameter, 'extended.filename', is for development purposes only and can
#      be ignored.
# 9.  Digest output files are named like this 
#                       isd_P60709_2-375_trypsin(0,1,2)_20240219162446.txt
#      or this
#                          isd_P60709_trypsin(0,1,2)_20240219162446.txt
#      where:
#       ------------------------------------------------------------------------------
#        Substring           Meaning
#       ------------------------------------------------------------------------------
#        isd                 Denotes an in silico digest.
#        P60709              Protein ID, as extracted from the FASTA header.
#        2-375               Range of indices passed from 'residues' (omitted if no
#                             range was specified).
#        trypsin(0,1,2)      Rule used for digestion, followed by the set of values
#                             used in the 'missedCleavages' argument.                 
#        20240219162446      Time stamp.
#       ------------------------------------------------------------------------------
#
# 10. Protein output file names are named like this 
#                               P60709_2-375_20240219162446.fasta
#      or this 
#                                  P60709_20240219162446.fasta
#      with substrings as defined in the previous note.
# 11. The lists of peptides generated by this script are non-redundant, on account of
#      the 'unique = TRUE' argument in the 'cleave' step.
# 12. Package 'cleaver' masks functions 'grepl', 'paste' and 'strsplit' from package 
#      'base', hence why calls to the respective 'base' functions have been made 
#      explicit.   

# Start time
start.time <- as.POSIXlt(Sys.time())
print(start.time)

# Options
options(warnPartialMatchAttr = TRUE)
options(warnPartialMatchDollar = TRUE)

# Get working directory path and add trailing separator if necessary
wd <- getwd()
if(substr(x = wd, start = nchar(wd), stop = nchar(wd)) != "/")
wd <- base::paste(wd, "/", sep = "")
print(base::paste("The current working directory is", wd), quote = FALSE)

# Set user-defined parameters
fasta.file.path <- base::paste(wd, "uniprotkb_2024_02_19_test-set.fasta", sep = "")
fasta.element <- 1 # index number
residues <- "all" # a string: either 'all' or a range such as '2-375'
save.residues <- FALSE # Should the (possibly truncated) sequence be exported?
rule <- "trypsin" # for 'enzym' argument
missed <- 0 # e.g., 'c(0, 1, 2)' for <=2 missed cleavages
min.peptide.length <- 2 # number of amino acid residues
output.folder.path <- wd # must include trailing path separator ('/' or '\\')
collapse.isd <- FALSE # should digest be collapsed into single string?

# Set development parameters
extend.isd.filename <- FALSE # Should 'collapse.isd' be represented in file name?
# | if you're not sure, go with 'FALSE'

# Load packages
library(seqinr)
library(cleaver)

# Hour, minutes and seconds for timestamp
if(nchar(as.character(start.time$hour)) == 1)
hour <- base::paste("0", as.character(start.time$hour), sep = "") else
hour <- as.character(start.time$hour)
if(nchar(as.character(start.time$min)) == 1)
minutes <- base::paste("0", as.character(start.time$min), sep = "") else
minutes <- as.character(start.time$min)
if(nchar(as.character(floor(start.time$sec))) == 1)
seconds <- base::paste("0", as.character(floor(start.time$sec)), sep = "") else
seconds <- as.character(floor(start.time$sec))

# Timestamp based on start time
time.stamp <- base::paste(
	gsub(pattern = "-", 
		replacement = "", 
		x = as.Date(as.POSIXlt(start.time))),
	hour, minutes, seconds, sep = "")
print(base::paste("Timestamp for this run:", time.stamp), quote = FALSE)

# Make substring (possibly 0 characters long!) via which to extend filename 
if(extend.isd.filename)
extension <- base::paste(substr(x = as.character(collapse.isd), start = 1, stop = 1), 
	"_", sep = "") else
extension <- ""

# Load *.fasta file
fasta.object <- read.fasta(fasta.file.path,
	seqtype = "AA",
	as.string = TRUE, # want string rather than vector
	whole.header = TRUE) # default is FALSE
print(base::paste("FASTA object contains ", 
	length(fasta.object), " sequences.", sep = ""),	quote = FALSE)

# Extract sequence to be digested and print its length
print(base::paste("Attempting to extract sequence number ", 
	fasta.element, ".", sep = ""), quote = FALSE)
protein <- fasta.object[[fasta.element]]
print(base::paste("A ", getLength(protein), 
	"-residue sequence was extracted.", sep = ""), quote = FALSE)

# Extract relevant metadata from FASTA header
fasta.header <- attr(x = protein, which = "name")
cat("FASTA header of extracted sequence:", "\n", fasta.header, "\n", sep = "")
bits <- base::strsplit(x = fasta.header, split = "|", fixed = TRUE)
bits <- unlist(bits)

# Different courses of action depending on the nature of 'residues'
if(residues == "all")
{termini <- c(1, nchar(protein)) # full-length
residues <- "" # omit range if protein is full-length
} else
if(base::grepl(pattern = "-", x = residues)) # hyphen as diagnostic of a range
{bits[3] <- residues # in preparation for 'write.fasta' step
termini <- as.integer(unlist(base::strsplit(residues, split = "-")))
residues <- base::paste(residues, "_", sep = "") # for correct file name format
} else
print("Error: couldn't understand which residues you wanted.", quote = FALSE)

# Process the protein sequence
protein <- substr(x = as.character(protein), 
	start = termini[1], 
	stop = termini[2])
print(base::paste("A ", nchar(protein), 
	"-residue sequence will be digested.", sep = ""), quote = FALSE)

# Save the protein sequence in its current form if instructed to do so
if(save.residues)
write.fasta(sequences = base::strsplit(protein, split = ""), # list is OK
	as.string = FALSE, # the default?
	names = base::paste(bits, collapse = "|"), # 'residues'-dependent!
	nbchar = 60, # the default
	file.out = base::paste(output.folder.path,
		bits[2], 
		"_", residues,
		time.stamp, ".fasta", sep = ""),
	open = "w") # the default

# Perform in silico digestion
isd <- cleave(x = protein, 
	enzym = rule, 
	missedCleavages = missed, 
	unique = TRUE)

# Filter the digest and print diagnostic counts
isd <- as.character(unlist(isd)) # prepare character vector
counts <- integer(2) # make empty vector in which to enter counts
names(counts) <- c("Pre", "Post") # add names
counts[1] <- length(isd) # enter initial unique-peptide count
isd <- isd[nchar(isd) >= min.peptide.length] # do the filtering
counts[2] <- length(isd) # enter final unique-peptide count
{print("Numbers of unique peptides, pre and post filtering:", quote = FALSE)
print(counts)} # print the vector

# Collapse 'isd' if instructed to do so
if(collapse.isd)
isd <- base::paste(isd, collapse = " ")

# Write 'isd' to file
write.table(x = isd, 
	file = base::paste(output.folder.path,
		"isd_", 
		bits[2], 
		"_", residues, rule, "(", 
		base::paste(missed, collapse = ","), ")_",
		extension,
		time.stamp, ".txt", sep = ""), 
	row.names = FALSE, 
	col.names = FALSE,
	quote = FALSE)

# Print session info
print(sessionInfo())

# Print system info
print(Sys.info())

# Print run end time
print(Sys.time())

# +-----------------------------------------------------------------------------------+
