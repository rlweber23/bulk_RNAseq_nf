suppressPackageStartupMessages({
    library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: tximport_counts.R <tx2gene_file> ",
        "<countsFromAbundance> <output_file>"
    )
}

tx2gene_file <- args[1]
counts_from_abundance <- args[2]
output_file <- args[3]

valid_counts_from_abundance <- c(
    "no",
    "scaledTPM",
    "lengthScaledTPM"
)

if (!counts_from_abundance %in% valid_counts_from_abundance) {
    stop(
        "Invalid countsFromAbundance: ", counts_from_abundance,
        "\nValid options: ",
        paste(valid_counts_from_abundance, collapse = ", ")
    )
}

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

cat(
    "Running tximport with countsFromAbundance =",
    counts_from_abundance,
    "\n"
)

cat(
    "Found",
    length(files),
    "Kallisto samples\n"
)

# Import and summarize transcript estimates to genes
txi <- tximport(
    files,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = counts_from_abundance,
    dropInfReps = TRUE,
    ignoreAfterBar = TRUE
)

# Write genes x samples matrix
write.table(
    txi$counts,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

