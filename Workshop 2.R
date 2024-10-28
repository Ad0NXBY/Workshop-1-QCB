#what type of elements are in the Celltype and Donor columns
class(sample.meta.data$Celltype)
class(sample.meta.data$Donor)

#Create separate vectors for Donor and Cell type and convert the elements to be factors
Donor <- factor(sample.meta.data$Donor)
Celltype <- factor(sample.meta.data$Celltype)

class(Celltype)
class(Donor)

#Create a model matrix called "design" for these samples using a single factor experimental design for Celltype
design <- model.matrix(~Celltype) #~ notation for the explanatory vairables you want to include in the model
design

#Update the dgelist.filtered.norm object with dispersion estimates. Use the experimental design from Q3
dgelist.filtered.norm <- estimateDisp(dgelist.filtered.norm, design = design)

#Fit a linear model to the data and dispersion estimates in dgelist.filtered.norm, using experimental design from Q3
#Use function glmQLFit, assign the output to an object called "fit"
fit <- glmQLFit(dgelist.filtered.norm, design)

#Take a look at the coefficients dataframe in the fit object
head(fit$coefficients)

#Carry out statistical testing using glmQLFTest() to assess if there is differential expression between M2 and M1 cell types
#Assign to output called "M2.v.M1" and look at its structure
M2.v.M1 <- glmQLFTest(fit, coef = 2) #2 is used as we are choosing the 2nd column and referencing it to the first column
head(M2.v.M1$table)

#Use topTags() function to order genes by p-value and assign table output to an object called DE
DE <- topTags(M2.v.M1, n = Inf)
head(DE)

#Join the gene meta data object to the DE table to add in the gene names info
DE <- left_join(rownames_to_column(DE$table, "gene_id"), gene.meta.data)
head(DE)

#Build a volcano plot
#log2 fold change v p-value
ggplot() + geom_point(data = DE, aes(x= logFC, y= PValue))

#log2 fold chanve v -log(10)p-value
ggplot() + geom_point(data = DE, aes(x= logFC, y= -log10(PValue)))

#Plot another volcano plot where genes with an FDR-adjusted p-value less than 0.05 are black and more than 0.05 are grey
ggplot() +
  geom_point(data = filter(DE, FDR <0.05), aes(x= logFC, y= -log10(PValue))) +
  geom_point(data = filter(DE, FDR >0.05), aes(x=logFC, y= -log10(PValue)), color = "grey") +
  theme_bw()

#Add another layer to the volcano plot from Q9 to label top 20 differentially expressed genes
library(ggrepel)
ggplot() +
  geom_point(data = filter(DE, FDR <0.05), aes(x= logFC, y= -log10(PValue))) +
  geom_point(data = filter(DE, FDR >0.05), aes(x=logFC, y= -log10(PValue)), color = "grey") +
  geom_text_repel(data = slice_head(DE, n=20), aes(x= logFC, y= -log10(PValue), label = gene_name))
  theme_bw()

#Rename cpm.filtered.norm to cpm
cpm <- cpm.dgelist.filtered.norm
head(cpm)

#Rename the cpm columns so that they provide information on both the Donor and Cell type for that sample
sample.meta.data <- sample.meta.data %>% mutate("sample_name"= paste0(Celltype, "_Donor", Donor))
head(sample.meta.data)
colnames(cpm) <- sample.meta.data$sample_name
head(cpm)

#Rename the rows of CPM to be gene symbols rather than gene ids
dim(gene.meta.data)
dim(cpm)
head(gene.meta.data)
head(rownames(cpm))

filtered.gene.meta.data <- left_join(data.frame("gene_id" = rownames(cpm)), gene.meta.data)
head(filtered.gene.meta.data)

rownames(cpm) <- filtered.gene.meta.data$gene_name

#Use Z-score scaling on the cpm data and assign this to an object called z.scaled.genes
z.scaled.genes <- t(cpm) %>% scale() %>% t()
head(z.scaled.genes)

#Find the Euclidean distance between samples and assign this to sample.scaled_distances
sample.scaled_distance <- dist(t(z.scaled.genes), method = "euclidean")
sample.scaled_distance

#Cluster samples using hierarchical clustering and function hclust(), use the "complete"linkage mothod and plot the output object
sample.scale_hclust <- hclust(sample.scaled_distance, method = "complete")
plot(sample.scale_hclust)

#Cluster the gene-wise scaled CPM values using a Euclidean distance measure and avearge linkage
gene_distance <- dist(z.scaled.genes, method = "euclidean")
gene_hclust <- hclust(gene_distance, method = "average")
plot(gene_hclust, labels = FALSE) #labs still give x row labels

#use cutree functionto cut the clustering tree from Q19 into 8 clusters
clusters.gene.k8 <- cutree(gene_hclust, k=8)
head(clusters.gene.k8)
table(clusters.gene.k8)

#subset the z-score scaled data for the genes that are in the gene cluster 3.
z.scaled.genes.cluster3 <- z.scaled.genes[clusters.gene.k8==3, ]

#check the dimensions of the subsetted data 
dim(z.scaled.genes)
dim(z.scaled.genes.cluster3)

#write the rownames to a csv file
write.csv(rownames(z.scaled.genes.cluster3), file = "cluster3_genenames.csv")

#Install ComplexHeatmap package using BiocManager
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#Plot a heatmap of gene cluster 3 values WITHOUT clustering rows or columns, do not plot gene names
Heatmap(z.scaled.genes.cluster3, cluster_rows= FALSE, cluster_columns = FALSE, show_row_names = FALSE)

#Plot the same data again, clustering columns and rows with the default linkage methods for sample and genes, DO NOT plot gene names
Heatmap(z.scaled.genes.cluster3, cluster_rows= TRUE, cluster_columns = TRUE, show_row_names = FALSE)

#Save the Heatmap as an image file
png(filename = "Cluster3.z.score_heatmap.png", height = 30, width = 10, units = "cm", res = 200)
ht <- Heatmap((z.scaled.genes.cluster3, 
               cluster_rows= TRUE,
               cluster_columns= TRUE,
               show_row_names= FALSE,
               height= unit(20,"cm"))
draw(ht)
dev.off()