### DESeq2 for lncRNA
library(DESeq2)
require(data.table)
library(rtracklayer)

#### rsem result merging
rsem.results = Sys.glob("/BiO/Research/UNIST-NSMF-2021-02/QC/CRSM/genome_bam/*.genes.results")
rsem.results = rsem.results[-grep("201209SHIN#18", rsem.results)] #remove 201209SHIN#18
rsem.table <- fread(rsem.results[1], header = T, sep = "\t")
merge.df = data.frame(rsem.table[, 1])
colnames(merge.df) <- "gene_id"


# merging gene expected_read_count of each sample
for(rsem.result in rsem.results){
    sample.id <- strsplit(basename(rsem.result), ".", fixed = TRUE)[[1]][1]
    rsem.df <- fread(rsem.result, header = T)
    #rsem.df <- rsem.df[, -"transcript_id(s)"]
    sub.df <- rsem.df[, c("gene_id", "expected_count")]
    colnames(sub.df) <- c("gene_id", sample.id)
    merge.df <- cbind(merge.df, sub.df[, 2])
}
merge.df <- unique(merge.df)
input.table = as.data.frame(merge.df)
#input.table <- as.data.frame(merge.df, row.names = merge.df$gene_id)
#input.table <- input.table[, 2:20]
input.table <- subset(input.table, rowSums(input.table[,2:20]) != 0)


#i = import("/BiO/Share/Database/GENCODE/M26/gencode.vM26.chr_patch_hapl_scaff.annotation.gtf")
#a = mcols(i)[, c("gene_id", "transcript_id", "gene_name")]
#write.table(a, "/BiO/Share/Database/GENCODE/M26/gencode.vM26.gene_id.transcript_id.gene_name", quote = FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)


gencode.file <- file.path("/BiO/Share/Database/GENCODE/M26/gencode.vM26.gene_id.transcript_id.gene_name")
gencode.file.df <- fread(file = gencode.file, sep = "\t")
colnames(gencode.file.df) <- c("gene_id", "transcript_id", "gene_name")
gencode.file.gene_id_name.df <- subset(gencode.file.df, select = c("gene_id", "gene_name"))
gencode.file.gene_id_name.df <- unique(gencode.file.gene_id_name.df)

cts <- merge(x = input.table, y = gencode.file.gene_id_name.df, by = "gene_id", all.x = TRUE)
cts$gene_id <- NULL
cts <- unique(cts)
cts.df <- cts[, c(20, 1:19)]
cts.df <- subset(cts.df, rowSums(cts.df[, 2:20]) != 0)
cts.df <- na.omit(cts.df)


unique.ids <- unique(cts.df$gene_name)
f.x <- function(x){
    unique.df <- subset(cts.df, gene_name == x)
    mean.row <- as.data.frame(t(colMeans(unique.df[, -1])))
    rownames(mean.row) <- x
    return(mean.row)
}
x <- do.call(rbind, lapply(unique.ids, f.x))

write.table(x, "/BiO/Research/UNIST-NSMF-2021-02/script/final.gene_name.read_count.txt", sep = "\t", quote = FALSE, append = FALSE)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------


#### loading input.table
#coldata = data.frame(fread("/BiO/Research/UNIST-NSMF-2021-02/script/coldata.NSMF_WT_KO.adenomas.txt", fill = T, header = T), row.names = 1)
coldata <- data.frame(fread("/BiO/Research/UNIST-NSMF-2021-02/script/coldata.adeno_vs_adjnor.NSMF_KO.txt", fill = T, header = T), row.names = 1)
coldata = subset(coldata, select = "condition")
coldata.df <- data.table(t(coldata))


cts <- data.frame(fread("/BiO/Research/UNIST-NSMF-2021-02/script/final.gene_name.read_count.txt"), row.names = 1)
colnames(cts) = gsub("X", "", colnames(cts))
colnames(cts) = gsub("N.", "N#", colnames(cts))
cts <- cts[, colnames(coldata.df)]
cts <- round(cts)
cts <- subset(cts, rowSums(cts) != 0)
cts <- cts + 1


#### deseq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

out_file = "/BiO/Research/UNIST-NSMF-2021-02/DEG/DESeq2/DESeq2_result_rsem.adeno_vs_adjnor.NSMF_KO"
write.table(res, out_file, sep = "\t", quote = FALSE, append = FALSE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

#### extract UP and DOWN expression
DESeq2.file = "/BiO/Research/UNIST-NSMF-2021-02/DEG/DESeq2/DESeq2_result_rsem.adeno_vs_adjnor.NSMF_KO"
DESeq2.df <- read.table(DESeq2.file, header = T, row.names = 1)

#DESeq2.df.pvalue.cut <- subset(DESeq2.df, pvalue != 'NA' & padj != 'NA' & pvalue < 0.01 & padj < 0.01)
DESeq2.df.pvalue.cut <- subset(DESeq2.df, pvalue != 'NA' & padj != 'NA' & pvalue < 0.05 & padj < 0.05)

#UP <- subset(DESeq2.df.pvalue.cut, log2FoldChange >= log2(3))
#DOWN <- subset(DESeq2.df.pvalue.cut, log2FoldChange <= -log2(3))
UP <- subset(DESeq2.df.pvalue.cut, log2FoldChange >= 1)
DOWN <- subset(DESeq2.df.pvalue.cut, log2FoldChange <= -1)

UP <- UP[order(-(UP$log2FoldChange)), ]
DOWN <- DOWN[order(DOWN$log2FoldChange), ]


## NSMF.KO_vs_WT.adenomas
outfile_UP = "/BiO/Research/UNIST-NSMF-2021-02/DEG/DESeq2/DESeq2_result_rsem.adeno_vs_adjnor.NSMF_KO.pval_0.05.l2fc2_1.UP"
write.table(UP, outfile_UP, sep = "\t", quote = FALSE, append = FALSE)
outfile_DOWN = "/BiO/Research/UNIST-NSMF-2021-02/DEG/DESeq2/DESeq2_result_rsem.adeno_vs_adjnor.NSMF_KO.pval_0.05.l2fc2_1.DOWN"
write.table(DOWN, outfile_DOWN, sep = "\t", quote = FALSE, append = FALSE)
