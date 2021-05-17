rm(list = ls())

library(optparse)
option_list = list(make_option(c("-r", "--readcount"), type="character", default=NULL, help="Read count data"),
                   make_option(c("-c", "--coldata"), type="character", default=NULL, help="Column data"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output data"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(readcount, cdata, output)
{
    library(DESeq2)
    require(data.table)

    coldata <- data.frame(fread(cdata, fill=TRUE, header=TRUE), row.names=1)
    coldata[["condition"]] <- factor(coldata[["condition"]])

    cts <- data.frame(fread(readcount), row.names=1)[, colnames(data.table(t(coldata)))]
    cts <- round(cts)
    cts <- subset(cts, rowSums(cts) != 0)
    cts <- cts + 1

    dds <- DESeq(DESeqDataSetFromMatrix(countData=cts, colData=coldata, design=~condition))
    res <- results(dds)

    write.table(res, output, sep="\t", quote=FALSE, append=FALSE)
}

if (length(opt) == 4)
{
    main(opt$readcount, opt$coldata, opt$output)
} else {
    print_help(opt_parser)
}
