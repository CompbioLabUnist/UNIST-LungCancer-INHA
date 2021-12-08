rm(list = ls())

library(optparse)
option_list = list(make_option(c("-i", "--input"), type="character", default=NULL, help="Input RDS file"),
                   make_option(c("-o", "--output"), type="character", default=NULL, help="Output TSV file"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input_file, output_file)
{
    input_data <- readRDS(input_file)
    write.table(input_data$results[[1]]$seg, file=output_file, quote=FALSE, sep="\t", row.names=FALSE)
}

if (length(opt) == 3)
{
    main(opt$input, opt$output)
} else {
    print_help(opt_parser)
}
