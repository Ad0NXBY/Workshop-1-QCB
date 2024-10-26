if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("recount3")
library(edgeR)
library(recount3)
#Following the workshop
#Make the tidyverse collection of packages available in your R session.
install.packages("tidyverse")
library(tidyverse)

#Read in the data from Session1_GSE81518_NormalizedCounts.txt and assign this to an object called indata. 
indata <- read.delim("Session1_GSE81518_NormalizedCounts.txt")

#Take a look at the first few rows
head(indata)

#How many genes are there?
dimension <- dim(indata)
nrow(indata)
#19284

#Extract the part of the sample name that categorizes the sample as originating from bone marrow or cerebrospinal fluid
sample_names <- colnames(indata[-1])
split_sample_names <- str_split_fixed(sample_names, "_", 2)

#Assign this vector to object called tissue
tissue <- split_sample_names[,2]

#Extract the part of the sample name that categorizes the patient
patient_number <- str_sub(sample_names,9,9)

#What are the total number of counts per sample?
sample_sums <- colSums(indata[ ,-1])

#How many genes have all zeros for the bone marrow samples or the cerebrospinal fluid samples?
counts_per_gene_sum <- data.frame("Gene" = indata[ ,1],
                                  "BM" = rowSums(indata[ , c(2,4,6)]),
                                  "CNS"= rowSums(indata[ , c(3,5,7)]))

#Try to make use of the pipe to pipe output of one command into the next, filter dataframe to find all zeros from BM and CNS sums
filter(counts_per_gene_sum, BM == 0 | CNS == 0) %>% nrow()

#Using the data in the current format, use base R graphics to plot boxplots of the count distribution per sample
boxplot(indata[, -1])

#Which gene/genes correspond to the very highly expressed outlier?
highly_expressed_genes <- which(indata[,-1] > 3000000, arr.ind = TRUE)
indata[highly_expressed_genes[,1],1]

#Transform the data matrix from wide format to a long format where every observation has its own row
long.indata <- indata %>%
  as.data.frame() %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Counts") %>%
  mutate("Patients" = str_sub(Sample,9,9),
         "Tissue" = str_split_fixed(Sample, "_",2)[,2])
head(long.indata)                                

#Replot the boxplot distribution from question 7 using the long format data and ggplot2 package
library(ggplot2)
ggplot(data = long.indata) + geom_boxplot(aes(x = Sample, y = Counts))

#Use appropriate data transformation to improve visualization
long.indata <- long.indata %>%
  mutate("Log10_counts" = log10(Counts))
ggplot(data = long.indata) + geom_boxplot(aes(x = Sample, y =Log10_counts))

#Filter the data for the gene VEGFA and plot a scatter plot of the counts per sample, splt by tissue(BM or CNS) and colored by patient
VEGFA_data <- filter(long.indata, Gene == "VEGFA")
ggplot(VEGFA_data,aes(x=Tissue, y=Counts, color = Patients)) + 
  geom_point()


#Load recount3 and edgeR library
library(recount3)
library(edgeR)

#Find SRA project ID from GEO accession number GSE36952 
#Use code snippet at the bottom of the page to retrieve data set
SRP012015.recount3 <- recount3::create_rse_manual(
  project = "SRP012015",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

#What class of R object is the data?
class(SRP012015.recount3)

#use functions Assay() and rowData() to look at the first 6 rows of the data matrix and the gene meta data
assay(SRP012015.recount3) %>% head()
rowData(SRP012015.recount3) %>% head()

#What are the column names in the sample meta data? why do you add colData()??
colnames(colData(SRP012015.recount3))

#What is the name of the assay for the data matrix that's been retrieved, and what are these values?
assayNames(SRP012015.recount3)
assay(SRP012015.recount3, "counts") <- compute_read_counts(SRP012015.recount3, round = FALSE)
assayNames(SRP012015.recount3)

#Create a dataframe called sample.meta.data containing only the SRR id and the sra.sample_attributes column
sample.meta.data <- colData(SRP012015.recount3) %>% as.data.frame() %>% select(external_id, sra.sample_attributes)
head(sample.meta.data)

#Extract the cell type category (M1 or M2) and donor category (1,2 or 3) from the sra.sample_attributes column
sample.meta.data<-sample.meta.data %>% 
  mutate("Celltype"=str_split_fixed(sra.sample_attributes, ";;", 3)[,2] %>% str_sub(start=1, end=2),
         "Donor"=str_split_fixed(sra.sample_attributes, ";;", 3)[,3] %>% str_sub(start=1, end=1)) %>%
  select(-sra.sample_attributes)

head(sample.meta.data)

#Create a dataframe called gene.meta.data containing only the gene id and the gene name
gene.meta.data <- rowData(SRP012015.recount3) %>% as.data.frame() %>% select(gene_id, gene_name)
head(gene.meta.data)

#Create a DGEList object called dgelist containing the count matrix data #don't understand this part
dgelist <- DGEList(assay(SRP012015.recount3, "counts"))
class(dgelist)
head(dgelist$counts)
head(dgelist$samples)

#how many genes have all zeros across the samples?
number.of.all.zero <- sum(rowSums(dgelist$counts)== 0)

#what fraction of the total number of genes is this
fraction.all.zero <- number.of.all.zero/sum(rowSum(dgelist$counts)) #why does the sum(rowSum) not work?
fraction.all.zero<- number.of.all.zero/dim(dgelist)[1]

#Calculate counts per million (cpm) and Log2 counts per million (lcpm) data matrices
cpm <- cpm(dgelist)
lcpm <- cpm(log(dgelist)) ??? non-numeric argument
lcpm <- cpm(dgelist, log = TRUE)
head(cpm)
head(lcpm)

#create a long (tidy) format fata frame for LogCPM values
long.lcpm <- as.data.frame(lcpm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogPCM", names_to = "SRR_ID")

head(long.lcpm)

#Plot sample distribution by SAMPLE
ggplot(long.lcpm) + geom_density(aes(LogPCM, color = SRR_ID)) #why is it like this can you explain this to me

#Obtain an indec of those genes that pass the expression filter using default parameters, using filterByExpr()
index.keep.expr <- filterByExpr(dgelist)

#Specifiy the group parameters in filterByExpr(), set this as the cell type from the sample.meta.data dataframe
index.keep.expr <- filterByExpr(dgelist, group = sample.meta.data$Celltype)
head(index.keep.expr)  

#Filter dgelist using filtering indx, assign it to a new object called dgelist.filtered
dgelist.filtered <- dgelist[index.keep.expr, , keep.lib.sizes = FALSE]
dim(dgelist.filtered)

#Plot sample distribitions of LogCPM for the filtered data
lcpm.filtered <- cpm(dgelist.filtered, log = TRUE)
long.lcpm.filtered <- as.data.frame(lcpm.filtered) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")
ggplot(long.lcpm.filtered) + geom_density(aes(LogCPM,color = SRR_ID))

#use calcNormFactors() functions to reate new DGEList object called dgelist.filtered.norm
dgelist.filtered.norm <- calcNormFactors(dgelist.filtered, method = "TMM")
dgelist.filtered.norm$samples$norm.factors

#Calculate CPM and Log2 CPM values from the dgelist.filtered.norm object
cpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm)
lcpm.dgelist.filtered.norm <- cpm(dgelist.filtered.norm, log = TRUE)

#Create long format dataframe called df.plotting containing these LogCPM and CPM values
long.cpm.dgelist.filtered.norm <- as.data.frame(cpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "CPM", names_to = "SRR_ID")

long.lcpm.dgelist.filtered.norm <- as.data.frame(lcpm.dgelist.filtered.norm) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, values_to = "LogCPM", names_to = "SRR_ID")

df.plotting <- full_join(long.cpm.dgelist.filtered.norm, long.lcpm.dgelist.filtered.norm)

dim(long.cpm.dgelist.filtered.norm)
dim(long.lcpm.dgelist.filtered.norm)
dim(df.plotting)

#Join the sample.meta.data and gene.meta.data dataframes to df.plotting
df.plotting <- left_join(df.plotting, sample.meta.data, by = c("SRR_ID"= "external_id"))
#head(df.plotting)
df.plotting <- left_join(df.plotting, gene.meta.data, by = c("Gene_ID" = "gene_id"))
head(df.plotting)

#Filter for the gene FCGR1A, plot the CPM values split by Celltype and colored by Donor
ggplot(df.plotting %>% filter(gene_name== "FCGR1A")) +
  geom_point(aes(x= Celltype, y= CPM, color = Donor)) + 
  ggtitle("FCGR1A")

#save following objects in an Rdata file 
save(dgelist.filtered.norm, sample.meta.data, gene.meta.data, cpm.dgelist.filtered.norm, lcpm.dgelist.filtered.norm, file = "Workshop1_output.Rdata")
