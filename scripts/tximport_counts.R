suppressPackageStartupMessages({
	    library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)

tx2gene_file <- args[1]
output_file  <- args[2]

# Read transcript -> gene map
tx2gene <- read.delim(
	tx2gene_file,
	header = FALSE,
	col.names = c("TXNAME", "GENEID")
)

# Each staged Kallisto directory contains abundance.tsv
files <- list.files(
	path = ".",
	pattern = "abundance.tsv$",
	recursive = TRUE,
	full.names = TRUE
)

if (length(files) == 0) {
	    stop("No Kallisto abundance.tsv files found")
}

# Sample name = directory containing abundance.tsv
sample_names <- basename(dirname(files))
names(files) <- sample_names

# Import and summarize transcript estimates to genes
txi <- tximport(
	files,
	type = "kallisto",
	tx2gene = tx2gene,
	countsFromAbundance = "no",
	ignoreAfterBar = TRUE,
	dropInfReps = TRUE
)

# Write genes x samples matrix
write.table(txi$counts,
	file = output_file,
	sep = "\t",
	quote = FALSE,
	col.names = NA
)
