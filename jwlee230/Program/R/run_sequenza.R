rm(list = ls())

library(optparse)
option_list = list(make_option(c("-s", "--sampleid"), type="character", default=NULL, help="sample name", metavar="character"),
                   make_option(c("-i", "--input"), type="character", default=NULL, help="input .seqz.gz file", metavar="character"),
                   make_option(c("-o", "--outdir"), type="character", default=NULL, help="output directory", metavar="character"),
                   make_option(c("-p", "--parallel"), type="integer", default=NULL, help="parallel number of CPUs", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(sampleid, input, outdir, parallel)
{
    library(sequenza)

    chrs = c(1:22, "X", "Y")
    chrs <- paste0("chr", chrs)

    test <- sequenza.extract(file=input, chromosome.list=chrs, parallel=parallel)
    CP <- sequenza.fit(sequenza.extract=test, mc.cores=parallel, chromosome.list=chrs, XY=c(X="chrX", Y="chrY"))
    sequenza.results(sequenza.extract=test, cp.table=CP, sample.id=sampleid, out.dir=outdir, chromosome.list=chrs, XY=c(X="chrX", Y="chrY"))
}

if (length(opt) == 5)
{
    main(opt$sampleid, opt$input, opt$outdir, opt$parallel);
} else {
    print_help(opt_parser)
}
