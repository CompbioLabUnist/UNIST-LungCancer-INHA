rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Input TSV data"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output data"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input, output)
{
    library(data.table)
    library(deconstructSigs)
    library(BSgenome.Hsapiens.UCSC.hg38)

    input_data <- read.table(file=input, sep="\t", header=TRUE)
    sigs_input <- mut.to.sigs.input(input_data, bsg=BSgenome.Hsapiens.UCSC.hg38)

    sigs_output <- data.frame()
    for(i in 1:nrow(sigs_input)){
        sample_info <- rownames(sigs_input)[i]
        sample_signature <- whichSignatures(tumor.ref=sigs_input, signatures.ref=signatures.cosmic, sample.id=sample_info, contexts.needed=TRUE, tri.counts.method="default")
        sigs_output <- rbind(sigs_output, sample_signature$weights)
    }

    write.table(sigs_output, file=output, quote=FALSE, sep="\t", col.names=NA)
}

if (length(opt) == 3)
{
    main(opt$input, opt$output)
} else {
    print_help(opt_parser)
}
