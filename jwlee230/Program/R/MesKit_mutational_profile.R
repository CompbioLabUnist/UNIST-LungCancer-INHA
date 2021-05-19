rm(list = ls())

library(optparse)
option_list = list(make_option(c("-m", "--maf"), type="character", default=NULL, help="Read count data"),
                   make_option(c("-c", "--clinical"), type="character", default=NULL, help="Column data"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output data"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(maf, clinical, output)
{
    library(MesKit)

    maf <- readMaf(mafFile=maf, clinicalFile=clinical, refBuild="hg38")
}

if (length(opt) == 4)
{
    main(opt$maf, opt$clinical, opt$output)
} else {
    print_help(opt_parser)
}
