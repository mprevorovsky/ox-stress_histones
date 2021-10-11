# This script produces Figures 4A, B
#####################################################

files <- read.delim('./files', header = FALSE, stringsAsFactors = FALSE)

data <- read.delim(files[1, 1], header = FALSE, stringsAsFactors = FALSE, skip = 2)
for (i in 2:nrow(files)){
  temp <- read.delim(files[i, 1], header = FALSE, stringsAsFactors = FALSE, skip = 2)
  data <- cbind(data, temp[, 2])
#  print(length(which(data[, 1] != temp[, 1]))) # the order of genes is the same in all files
}

colnames(data) <- c('gene', 
                    'WT_0', 'WT_15', 'WT_60',
                    'cbf11_0', 'cbf11_15', 'cbf11_60',
                    'cbf12_0', 'cbf12_15', 'cbf12_60')

# normalize all samples to WT t=0min
data.norm <- data
for (i in 2:ncol(data.norm)){
  data.norm[, i] <- data.norm[, i] / data[, 2]
}

# filter out genes with too many missing data points that would otherwise impair heatmap construction
x <- data.norm[, 2:ncol(data.norm)]
rows <- nrow(x)
cols <- ncol(x)
max.NA <- 0.4 # maximum ratio of missing values per row
keep <- NA
for (i in 1:rows)
{
  NAs <- 0
  for (j in 1:cols)
  {
    if (is.na(x[i, j]) == TRUE)
    {
      NAs <- NAs + 1
    }
  }
  NAs.rat <- NAs / cols
  if (NAs.rat <= max.NA)
  {
    keep[length(keep) + 1] <- i
  }
}
keep <- keep[2:length(keep)]
data.norm <- data.norm[keep, ]
rm(x, keep)

# filter for DEGs
threshold.down <- 0.5
threshold.up <- 2

keep <- NULL
for (i in 1:nrow(data.norm)){
  if (min(data.norm[i, 2:ncol(data.norm)], na.rm = TRUE) <= threshold.down){
    keep <- c(keep, data.norm[i, 1])
  }
}
for (i in 1:nrow(data.norm)){
  if (max(data.norm[i, 2:ncol(data.norm)], na.rm = TRUE) >= threshold.up){
    keep <- c(keep, data.norm[i, 1])
  }
}
keep <- unique(keep)
data.diff.norm <- data.norm[data.norm[, 1] %in% keep, ]
rm(keep)

# plot Figure 4A
library(gplots)
pdf('peroxide_timecourse_DEGs.pdf', height = 7, width = 3)
heatmap.2(log2(as.matrix(data.diff.norm[, 2:ncol(data.diff.norm)])), 
          dendrogram = 'row', 
          Colv = FALSE, 
          trace = 'none',
          scale = 'none', 
          symbreaks = TRUE, 
          density.info = 'none',
          breaks = seq(-4, 4, 0.4),
          margins = c(7, 3),
          colsep = c(3, 6),
          sepwidth = c(0.05, 0.05),
          keysize = 2,
          col = colorRampPalette(c("blue", 'black', 'yellow'))(20))
dev.off()

# plot Figure 4B
length(which(data.norm[, 5] >= threshold.up)) # highest hit is kanMX
length(which(data.norm[, 5] <= threshold.down)) # lowest hit is cbf11
data.norm <- data.norm[!data.norm$gene %in% c('c736.08', 'kan-mx'), ]

pdf('peroxide_timecourse_Cbf11_DEGs_stripchart.pdf', height = 2.5, width = 5)
stripchart(c(log2(data.norm[which(data.norm[, 5] >= threshold.down & data.norm[, 5] <= 1), 5]), 
             log2(data.norm[which(data.norm[, 5] <= threshold.up & data.norm[, 5] > 1), 5])), 
           method = 'jitter', xlim = c(-5, 6), pch = 21, bg = rgb(0.5, 0.5, 0.5, 1), 
           xlab = 'log2FCE cbf11 vs WT time 0',
           las = 1, col = 'grey50', cex= 0.5)
stripchart(log2(data.norm[which(data.norm[, 5] <= threshold.down), 5]), 
           method = 'jitter', pch = 21, bg = rgb(0, 0, 1, 1), 
           add = TRUE, col = 'grey50', cex= 0.5)
stripchart(log2(data.norm[which(data.norm[, 5] >= threshold.up), 5]), 
           method = 'jitter', pch = 21, bg = rgb(1, 1, 0, 1), 
           add = TRUE, col = 'grey50', cex= 0.5)
abline(v = c(-1, 0, 1), col = 'grey50', lty = 2)

stripchart(c(log2(data.norm[which(data.norm[, 4] >= threshold.down & data.norm[, 4] <= 1), 4]), 
             log2(data.norm[which(data.norm[, 4] <= threshold.up & data.norm[, 4] > 1), 4])), 
           method = 'jitter', xlim = c(-5, 6), pch = 21, bg = rgb(0.5, 0.5, 0.5, 1), 
           xlab = 'log2FCE WT H2O2 vs WT time 0',
           las = 1, col = 'grey50', cex= 0.5)
stripchart(log2(data.norm[which(data.norm[, 4] <= threshold.down), 4]), 
           method = 'jitter', pch = 21, bg = rgb(0, 0, 1, 1), 
           add = TRUE, col = 'grey50', cex= 0.5)
stripchart(log2(data.norm[which(data.norm[, 4] >= threshold.up), 4]), 
           method = 'jitter', pch = 21, bg = rgb(1, 1, 0, 1), 
           add = TRUE, col = 'grey50', cex= 0.5)
abline(v = c(-1, 0, 1), col = 'grey50', lty = 2)
dev.off()

