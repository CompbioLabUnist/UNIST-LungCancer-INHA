rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Input TSV file"),
                   make_option(c("-a", "--annotation"), type="character", default=NULL, help="Reference annotation TXT file"),
                   make_option(c("-r", "--reference"), type="character", default=NULL, help="Reference TSV file"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output TSV file"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input_file, annotation_file, reference_file, output_file)
{
    install.packages("BisqueRNA", repos="http://cran.us.r-project.org")
    library(Biobase)
    library(BisqueRNA)
    require(data.table)

    input_data <- read.table(file=input_file, sep="\t", header=TRUE, row.names=1)
    input_eset <- new("ExpressionSet", exprs=as.matrix(input_data))

    annotation_data <- na.omit(read.table(file=annotation_file, sep="\t", header=TRUE, row.names=1))
    annotation_data <- annotation_data[annotation_data$Sample_Origin %in% c("nLung", "tLung"), ]
    barcodes <- order(rownames(annotation_data))
    annotation_data <- annotation_data[barcodes, ]
    print(head(annotation_data))

    reference_data <- as.data.frame(readRDS(reference_file))
    reference_data <- reference_data[, barcodes]

    reference_eset <- ExpressionSet(as.matrix(reference_data), phenoData=AnnotatedDataFrame(data.frame(row.names=colnames(reference_data), samples=annotation_data[["Sample"]], clusters=annotation_data[["Cell_subtype"]])))
    print(reference_eset)

    estimation_proportion <- ReferenceBasedDecomposition(bulk.eset=input_eset, sc.eset=reference_eset, cell.types="clusters", subject.names="samples", verbose=TRUE, use.overlap=FALSE)
    output_data <- estimation_proportion$bulk.props

    print(head(output_data))
    write.table(output_data, file=output_file, quote=FALSE, sep="\t", col.names=NA)
}

if (length(opt) == 5)
{
    main(opt$input, opt$annotation, opt$reference, opt$output)
} else {
    print_help(opt_parser)
}
