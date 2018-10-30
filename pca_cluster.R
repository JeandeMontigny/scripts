library(ggplot2)

# ---------------- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2){
	stop("arg error - need 2 arg: [path to measures_output files] [nb of cluster]")
}
setwd(args[1])

# ---------------- #
#nbOfLine <- system(" cat measures_output.txt | wc -l")
df <- read.table("measures_output.txt", header=T) # column separated by " "
#df
pca <- prcomp(df[,2:5]) # scaled pca excluding cell name and cell type column
scores <- pca$x[,1:3] # scores for first three PC's

# k-means clustering; centers=nb of cluster
cluster_nb <- as.numeric(args[2])
km <- kmeans(scores, centers=cluster_nb, nstart=10) #algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
ggdata <- data.frame(scores, Cluster=km$cluster, cell_type=df$cell_type)
	ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2, color=factor(Cluster)), size=5, shape=20) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster"))

# ---------------- #
message("done")

