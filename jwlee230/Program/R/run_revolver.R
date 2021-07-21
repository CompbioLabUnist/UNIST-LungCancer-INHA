# https://caravagnalab.github.io/revolver/articles/Cross_CRC_MSeq.html
rm(list = ls())

library(optparse)
option_list = list(make_option(c("--input"), type="character", default=NULL, help="Input TSV data"),
                    make_option(c("--cluster"), type="character", default=NULL, help="Input TSV data"),
                    make_option(c("--driver"), type="character", default=NULL, help="Input TSV data"),
                    make_option(c("--dendrogram"), type="character", default=NULL, help="Input TSV data"),
                    make_option(c("--det"), type="character", default=NULL, help="Input TSV data"),
                    make_option(c("--clonality"), type="character", default=NULL, help="Input TSV data"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

main <- function(input_file, cluster_file, driver_file, dendrogram_file, det_file, clonality_file)
{
    library(revolver)
    require(data.table)

    input_data <- read.table(file=input_file, sep="\t", header=TRUE, as.is=TRUE)

    revolver_input_data <- revolver_cohort(dataset=input_data, ONLY.DRIVER=FALSE, MIN.CLUSTER.SIZE=0, annotation="")
    revolver_check_cohort(revolver_input_data)

    non_recurrent <-  Stats_drivers(revolver_input_data) %>% filter(N_tot == 1) %>% pull(variantID)
    revolver_input_data <- remove_drivers(revolver_input_data, non_recurrent)

    revolver_input_data <- compute_mutation_trees(revolver_input_data)

    revolver_input_data <- revolver_fit(revolver_input_data, parallel=FALSE, initial.solution=NA)

    revolver_input_data <- revolver_cluster(revolver_input_data, min.group.size=1)

    plot_clusters(revolver_input_data, cutoff_trajectories=1, cutoff_drivers=0, arrow.symbol=" >> ")
    ggsave(filename=cluster_file, device="pdf", width=9, height=20, units="in")

    plot_drivers_graph(revolver_input_data)
    ggsave(filename=driver_file, device="pdf", width=9, height=9, units="in")

    plot_dendrogram(revolver_input_data)
    ggsave(filename=dendrogram_file, device="pdf")

    plot_DET_index(revolver_input_data)
    ggsave(filename=det_file, device="pdf", width=9, height=16, units="in")

    plot_drivers_clonality(revolver_input_data)
    ggsave(filename=clonality_file, device="pdf", width=9, height=16, units="in")
}

if (length(opt) == 7)
{
    main(opt$input, opt$cluster, opt$driver, opt$dendrogram, opt$det, opt$clonality)
} else {
    print_help(opt_parser)
}
