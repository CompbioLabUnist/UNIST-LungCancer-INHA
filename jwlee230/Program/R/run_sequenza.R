rm(list = ls())

library(optparse)
option_list = list(make_option(c("-s", "--sampleid"), type="character", default=NULL, help="sample name", metavar="character"), make_option(c("-i", "--input"), type="character", default=NULL, help="input .seqz.gz file", metavar="character"), make_option(c("-o", "--outdir"), type="character", default=NULL, help="output directory", metavar="character"), make_option(c("-p", "--parallel"), type="integer", default=NULL, help="output directory", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(sampleid, input, outdir, parallel)
{
    library(sequenza)

    test <- sequenza.extract(input, assembly="hg19", parallel=parallel)
    CP <- sequenza.fit(test, mc.cores=parallel)
    sequenza.results(sequenza.extract=test, cp.table=CP, sample.id=sampleid, out.dir=outdir)
}

if (length(opt) == 5)
{
    main(opt$sampleid, opt$input, opt$outdir, opt$parallel);
}
else
{
    print_help(opt_parser)
}
