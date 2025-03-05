setwd("/home/yuq22/ihb-intestine-evo/examples/genomic_features")
data <- read.table("/home/yuq22/ihb-intestine-evo/Annotation/UCSC/hg38/phyloP100way/LCT_MCM6.hg38.phyloP100way.bedGraph", sep="\t", stringsAsFactors=F)
colnames(data) <- c("chr","start","end","score")

library(GenomicRanges)
phyloP_gr <- makeGRangesFromDataFrame(data, keep.extra.columns=TRUE)

# Create a range covering the entire set of regions
total_range <- range(phyloP_gr)
# Define your sliding window parameters
n1 <- 500
window_size <- round(width(total_range)/n1)
step_size <- round(window_size/2)

# Create sliding windows across the total range
sliding_windows <- slidingWindows(total_range, width = window_size, step = step_size)[[1]]

# Find overlaps between the phyloP scores and the sliding windows
hits <- findOverlaps(query=phyloP_gr, subject=sliding_windows)
p <- Pairs(phyloP_gr, sliding_windows, hits = hits)
i <- pintersect(p)
df <- as.data.frame(i)
df$window_idx <- subjectHits(hits)

# Calculate the average phyloP score for each window
# and get the genomic coordinate of center of each window
window_averages <- t(sapply(sort(unique(df$window_idx)), function(i) {
  idx <- which(df$window_idx==i)
  score <- sum(df$score[idx] * df$width[idx], na.rm=T)/sum(df$width[idx])
  range <- c(df$start[idx], df$end[idx])
  coor <- (min(range)+max(range))/2
  return(c(score, coor))
}))
colnames(window_averages) <- c("score", "coor")

data <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Table_S6_hg38_coor.rds")
aa <- findOverlaps(subject=data, phyloP_gr)



pdf("Plot_scatterPlot_average_phyloP_around_LCT_and_MCM6.pdf")
plot(window_averages[,"coor"], window_averages[,"score"], pch=16, cex=0.5, col="#696969", ylab="Window average phyloP", xlab="hg38 coordinate", bty="n")
abline(v=135851076, lty=2)
text(x=135851076, y = max(window_averages[,"score"])*0.9, labels="rs4988235 / MCM6" )
dev.off()

# visualize the CMS data with something like Mahanton plot
setwd("/home/yuq22/ihb-intestine-evo/examples/genomic_features")
data <- readRDS("/home/yuq22/ihb-intestine-evo/evo_signature/CMS_positive_selection_Grosmann_2010_Science/Table_S6_hg38_coor.rds")
df <- as.data.frame(data, stringsAsFactors=F)
input <- df[,c("SNP.ID", "seqnames", "start", "Normalized.CMS")]
names(input)[1:3] <- c("SNP", "CHR", "BP")
input$CHR <- as.numeric(sub("chr", "", input$CHR))
pdf("Plot_manhattan_plot_all_chr_CMS.pdf", height=5)
manhattan(input, p = "Normalized.CMS", logp = FALSE, ylab = "Normalized CMS", genomewideline = FALSE,
    suggestiveline = FALSE, main = "Positive selection")
dev.off()

pdf("Plot_manhattan_plot_chr2_CMS.pdf", height=5)
manhattan(subset(input,CHR == 2), p = "Normalized.CMS", logp = FALSE, ylab = "Normalized CMS", genomewideline = FALSE,
    suggestiveline = FALSE, main = "Positive selection", highlight = "rs4988235")
dev.off()
