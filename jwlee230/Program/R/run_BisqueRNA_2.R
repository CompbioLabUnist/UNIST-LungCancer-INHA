rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Input TSV file"),
                   make_option(c("-r", "--reference"), type="character", default=NULL, help="Reference RDS file"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output TSV file"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input_file, reference_file, output_file)
{
    library(Biobase)
    library(BisqueRNA)
    require(data.table)
    library(Seurat)

    input_data <- read.table(file=input_file, sep="\t", header=TRUE, row.names=1)
    input_eset <- new("ExpressionSet", exprs=as.matrix(input_data))

    load(reference_file)
    reference_data <- as.data.frame(eleven.tils.cd3.integrated@assays$RNA@counts)
    print(head(reference_data))

    reference_eset <- ExpressionSet(as.matrix(reference_data), phenoData=AnnotatedDataFrame(data.frame(row.names=colnames(reference_data), samples=eleven.tils.cd3.integrated@meta.data[["orig.ident"]], clusters=Idents(eleven.tils.cd3.integrated))))
    print(reference_eset)

    estimation_proportion <- ReferenceBasedDecomposition(bulk.eset=input_eset, sc.eset=reference_eset, cell.types="clusters", subject.names="samples", verbose=TRUE, use.overlap=FALSE)
    output_data <- estimation_proportion$bulk.props

    print(head(output_data))
    write.table(t(output_data), file=output_file, quote=FALSE, sep="\t", col.names=NA)
}

if (length(opt) == 4)
{
    main(opt$input, opt$reference, opt$output)
} else {
    print_help(opt_parser)
}
