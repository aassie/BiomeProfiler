#Wormbiome R script
#Last update 14/01/25 -- V1.0

message(sprintf("Running Wormbiom R script\n"))
options(error = quote({dump.frames(to.file=TRUE); q()}))
set.seed(123)

## Get script arguments
args <- commandArgs(trailingOnly=TRUE)

## Set data location for loading
Data.location <- args[1]      # ASVTable
Meta.location <- args[2]      # ASVMeta  
Bin.location <- args[3]       # BlastResult
out.folder <- paste("./",args[4],"/",sep="")  # Output folder
GenomeDB.location <- args[5]  # Genome Metadata
KO.meta.location <- args[6]   # KeggMeta
WB.location <- args[7]        # Genome Database
selecteur1 <- args[8]         # Column to use as comparator
SSU_Count.location <- args[9] # SSUCount
VERBOSE <- args[10]           # DEBUG
STATS <- args[11]             # Performing stats on table
count.threshold <- as.numeric(args[12]) # Count Threshold for removing a sample
Tableopt1 <- args[13]         # ASV Count table type
SourceDir <- args[14]         # Directory with accessory scripts

## Load extra scripts
source(SourceDir)

# Check for required packages
has_deseq2 <- requireNamespace("DESeq2", quietly = TRUE)
has_phyloseq <- requireNamespace("phyloseq", quietly = TRUE)
has_biocyc <- requireNamespace("BioCyc", quietly = TRUE)

# Check for required packages and stop if missing
missing_pkgs <- c()
if(!has_deseq2) missing_pkgs <- c(missing_pkgs, "DESeq2")  
if(!has_phyloseq) missing_pkgs <- c(missing_pkgs, "phyloseq")

if(length(missing_pkgs) > 0) {
    stop(paste0("Required packages missing. Please install:\n",
                paste(" -", missing_pkgs, collapse="\n"), 
                "\n\nUse BiocManager::install() to install these packages."))
}

# Warn about optional packages
if(!has_biocyc && VERBOSE==1) {
    stop(cat("Note: BioCyc package not installed - please install the package\n"))
}

# Load the core libraries needed
suppressMessages(require(readr))
suppressMessages(require(tidyverse))
suppressMessages(require(reshape2))
suppressMessages(require(vegan))
suppressMessages(require(RColorBrewer))
suppressMessages(require(labdsv))

message(sprintf("Libraries loaded\n"))

# Setting up analysis pipeline

cat("------------------\n")
cat("Preparing Analysis\n")
cat("------------------\n")
cat("\n")

# Update debug message to match new argument order
if(VERBOSE==1){
    message(paste("c(\"",
        args[1],"\",\"",  # ASVTable
        args[2],"\",\"",  # ASVMeta
        args[3],"\",\"",  # BlastResult
        args[4],"\",\"",  # Output folder
        args[5],"\",\"",  # Genome Metadata
        args[6],"\",\"",  # KeggMeta
        args[7],"\",\"",  # Genome Database
        args[8],"\",\"",  # Sel1
        args[9],"\",\"",  # SSUCount
        args[10],"\",\"", # DEBUG
        args[11],"\",\"", # Performing stats on table
        args[12],"\",\"", # Count Threshold for removing a sample
        args[13],"\",\"", # ASV Count table type
        args[14],"\")",   # Directory with accessory scripts
        sep=""))
}

## Create output folder
system(paste("mkdir -p", out.folder))
system(paste("mkdir -p ", out.folder,"/log", sep=""))
# Core script
# Let's go!

## Load the data
### ASV

if(VERBOSE==1){cat("Loading MetaData\n")}
M<-read.csv(Meta.location, header = T, sep="\t")

### Verify Count table
if(VERBOSE==1){cat("Loading Data\n")}

### Check input table type
if(Tableopt1=="BIOM"){
D<-read.csv(Data.location, header=T, sep="\t", row.names = 1, skip=1)
}

if(Tableopt1=="TSV"){
  D<-read.csv(Data.location, header=T, sep="\t", row.names = 1)
}

if(Tableopt1=="CSV"){
  D<-read.csv(Data.location, header=T, row.names = 1)
}

colnames(M)[1]<-"X.SampleID"
M$X.SampleID<-make.names(M$X.SampleID)

### Checking if samples are present in meta data and vice versa

if(VERBOSE==1){cat("\tCross checking data and metadata\n")}

Missing.sample.1<-M[!M$X.SampleID %in% colnames(D),1]
M<-M[M$X.SampleID %in% colnames(D),]

if(length(Missing.sample.1)>0){
  cat("* The following sample are in the metadata but not in the dataset\n")
  print(Missing.sample.1)
  write.table(Missing.sample.1, paste(out.folder,"/log/Sample.No.Data.txt", sep=""))
}

Missing.sample.2<-colnames(D[,!colnames(D) %in% M$X.SampleID])
D<-D[,colnames(D) %in% M$X.SampleID]

if(length(Missing.sample.2)>0){
  cat("* The following sample are in the dataset but not in the metadata\n")
  print(Missing.sample.2)
  write.table(Missing.sample.2, paste(out.folder,"/log/Sample.No.Metadata.txt", sep=""))
}

### Filter samples

#Removing empty row and column

null.row<-rownames(D[rowSums(D)==0,])
null.col<-colnames(D[,colSums(D)==0])

if(length(null.row)>0){
  cat(paste("\t/!\\ The following ASV sequences have no count:\n"))
  print(null.row)
  cat(paste("They will be removed now\n"))
  write.table(null.row, paste(out.folder,"/log/ASV.0.count.txt", sep=""))
}

if(length(null.col)>0){
  cat(paste("\t/!\\ The following samples are empty\n"))
  print(null.col)
  cat(paste("They will be removed now\n"))
  write.table(null.col, paste(out.folder,"/log/Sample.0.count.txt", sep=""))
}

D<-D[!rownames(D) %in% null.row,]
D<-D[,!colnames(D) %in% null.col]

##Removing sample below threshold

D.tmp<-D[,colSums(D)>count.threshold]
Sample.filteredout.read<-M[!M$X.SampleID %in% colnames(D.tmp),1]

if(length(Sample.filteredout.read)>0){
  cat(paste("\t/!\ The following were below the count threshold of", count.threshold, "reads\n"))
  print((Sample.filteredout.read))
  cat(paste("They will be removed now\n"))
  write.table(Sample.filteredout.read, paste(out.folder,"/log/ASV.under.count.threshold.txt", sep=""))
}

D<-D.tmp
M<-M[M$X.SampleID %in% colnames(D.tmp),]

#Checking if there are still samples
if(length(ncol(D))==0){
  cat(paste("\t/!\\ There are no sample in your dataset\n"))
  cat(paste("I will stop now, please double check your data\n\n"))
  Gameover()
  quit(save="ask")
}

### Load Adjust SSU count
if(VERBOSE==1){cat("SSU Count module, currently deactivated\n")}
SSU_Count<-read.csv(SSU_Count.location, sep="\t", header = T, row.names = NULL)
#SSU_Count$SSU_Count<-gsub(0,1,SSU_Count$SSU_Count)
SSU_Count$SSU_Count<-1

### Load Genome database
if(VERBOSE==1){cat("Loading Genome Database\n")}

WB<-read_tsv(WB.location)
KO.meta<-read_tsv(KO.meta.location)

### Genome metadata
if(VERBOSE==1){cat("Loading Genome Metadata\n")}
GenomeDB<-read_tsv(GenomeDB.location)

unique.db.genome<-unique(WB$Genome)

## Adjustments
### ASV to relative abundance
cat("----------------\n")
cat("ASV Summary Step\n")
cat("----------------\n")
cat("\n")
if(VERBOSE==1){cat("Calculating Relative Abundance of ASV table\n")}
if(VERBOSE==1){cat("----------------\n")}
D<-D[complete.cases(D),]
D1<-as.data.frame(t(prop.table(as.matrix(D), margin=2)))

### Selection of correct samples
if(VERBOSE==1){cat("Merging Metadata\n")}

selectors<-c("X.SampleID",selecteur1)
D2<- as.data.frame(t(D)) %>% 
  rownames_to_column("X.SampleID") %>% 
  left_join(M %>% select(all_of(selectors))) 
## Plots

### Taxonomy plots
if(VERBOSE==1){cat("Plotting ASV\n")}
if(VERBOSE==1){cat("----------------\n")}

if(VERBOSE==1){
  print("Sanity check")
  print(head(D1))
}

Bplot<-D1 %>% 
  rownames_to_column("X.SampleID") %>%
  pivot_longer(!X.SampleID, names_to = "variable", values_to = "value") %>% 
  left_join(M) %>% 
  mutate(variable=factor(as.character(variable)))

colourCount = length(unique(Bplot$variable))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

ASV.barplot<-ggplot(Bplot, aes(x=X.SampleID, y=value, fill=variable)) +
  theme_bw() +
  geom_bar(stat = "identity", position="fill") +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  facet_wrap(as.formula(paste("~",selecteur1)), scale="free",nrow=1) + 
  scale_fill_manual(values = getPalette(colourCount)) +
  guides(fill=guide_legend(ncol=1))+
  xlab("Samples")+
  ylab("Relative ASV abundance in %") +
  ggtitle(paste("Relative ASV abundance split by",selecteur1))

ggsave(paste(out.folder, "ASV.png", sep=""), ASV.barplot, device = "png", width = 20, height = 10, dpi=150)

### PCoA
if(VERBOSE==1){cat("Plotting ASV PCoA\n")}
if(VERBOSE==1){cat("-----------------\n")}
selectors2<-c(colnames(D2)[1],selecteur1)
Fdist<-vegdist(D2[,!colnames(D2) %in% selectors2], method="jaccard")
Fpco<-pco(Fdist, k=10)
eig1<- round(Fpco$eig[1]/sum(Fpco$eig)*100,2)
eig2<- round(Fpco$eig[2]/sum(Fpco$eig)*100,2)
PCOplot<-as.data.frame(Fpco$points) %>% 
  mutate(X.SampleID=D2$X.SampleID) %>% 
  left_join(M)

colourCount = length(unique(PCOplot[[selecteur1]]))
ASV.PCoa<-ggplot(PCOplot, aes(x=V1, y=V2, col=get(selecteur1))) +
  geom_point() +
  theme_bw()+
  xlab(paste("Dimension 1", eig1, "%",sep=" "))+
  ylab(paste("Dimension 2", eig2, "%",sep=" "))+
  ggtitle("Principal Coordinates Analysis of ASV table")+
  scale_color_manual(values = getPalette(colourCount), name=paste(selecteur1))

ggsave(paste(out.folder, "PCOA.ASV.results.png", sep=""), ASV.PCoa, device = "png", width = 10, height = 10, dpi=150)

## ASV to Genomes

### Bin
if(VERBOSE==1){cat("Loading Blast Results\n")}
Bin<-read.csv(Bin.location,header = F, sep="\t")
colnames(Bin)<-c("ASV","Genome","pID")
Bin<-Bin[,c(1,2,3)]
Bin$Genome<-gsub("BIGb0243","BIGb0244",Bin$Genome)

### List of missing genomes

# genome present
Genome.present<-unique(Bin[Bin$Genome %in% unique(unique.db.genome),2])

# genome missing
Genome.absent<-unique(Bin[!Bin$Genome %in% unique(unique.db.genome),2])

if(length(Genome.absent)>=1){
  write.table(Genome.absent, paste(out.folder,"/log/Genome.not.present.in.database.txt", sep=""))
}

Genome.ratio<-paste(length(Genome.absent)," out of ",length(Genome.present)," genomes from your collection are missing from our database", sep="")
if(VERBOSE==1){cat(Genome.ratio,"\n")}
if(VERBOSE==1){cat("\n")}

### Selecting bins

Bin<-unique(Bin)
Bin.select<-tibble()
Bin.log<-tibble() #Add a bin log

for (i in unique(Bin$Genome)){
  cat(i,"\n")
  if(length(Bin[Bin$Genome==i&Bin$pID==100,1])==0){
    cat(" - No direct Match\n")
    a=max(Bin[Bin$Genome==i,3])
    ASV.temp.table<-Bin[Bin$Genome==i & Bin$pID==a,]
    if (nrow(ASV.temp.table)==1){
      b<-ASV.temp.table[1,1]
      cat("\t Unique ASV/Genome\n")
      cat("Closest hit: ",a,"% identity\n")
      cat(" Matching ASV: ")
      cat(b,"\n")
      Bin.select<-rbind(Bin.select, Bin[(Bin$ASV %in% b) & (Bin$Genome %in% i),])
    }else{
      cat("\t Multiple Match, skipping\n")
    }
  }else{
    cat("\t Direct Match\n")
    Bin.select<-rbind(Bin.select, Bin[Bin$Genome==i&Bin$pID==100,])
  }
}

Bin<-Bin.select

#Which ASV don't map to anything and how much they represent
D.RS<-as.data.frame(rowSums(D))
colnames(D.RS)<-"Sum"
D.RS$ASV<-rownames(D.RS)
D.RS$Per<-(D.RS$Sum/sum(D.RS$Sum))*100
D.RS$Per<-round(D.RS$Per, digits = 2)

SubRS<-D.RS[!D.RS$ASV %in% Bin$ASV,c(2,3)]
ASV.Missing<-paste("The ASV not mapping to any selected genomes represent ",sum(SubRS[,2]),"% of the dataset", sep="")

if(nrow(SubRS)>=1){
  write.table(SubRS, paste(out.folder,"/log/ASV.not.mapping.to.genome.txt", sep=""))
}

if(VERBOSE==1){cat(ASV.Missing,"\n")}
if(VERBOSE==1){cat("----------------\n")}
if(VERBOSE==1){cat("\n")}

### Keep the genome we know we have
if(VERBOSE==1){cat("Filtering with available genomes\n")}
if(VERBOSE==1){cat("----------------\n")}
if(VERBOSE==1){cat("\n")}
Bin.backup<-Bin
Bin<-Bin[Bin$Genome %in% unique(unique.db.genome),]

### How much what we have represent:
D1b<-t(D1[,-ncol(D1)])
MeanASV<-(sum(colSums(D[rownames(D) %in% Bin$ASV,]))/sum(colSums(D[complete.cases(D),])))*100

Genome.div<-paste("The sequenced data set represents ", round(MeanASV,digits = 1),"% of the diversity")
if(VERBOSE==1){cat(Genome.div,"\n")}
if(VERBOSE==1){cat("----------------\n")}
if(VERBOSE==1){cat("")}
### Reploting with genomes and corrected ASV
if(VERBOSE==1){cat("Plotting with Genome Information\n")}
if(VERBOSE==1){cat("----------------\n")}
Da<-D
Da<-Da[rowSums(Da)>0,]
ASV.noCount<-rownames(Da[!rowSums(Da)>0,])
Single<-tibble()
All<-tibble()
ASV.totmissing<-tibble()

for (i in unique(rownames(Da))) {
  Bin.length<-length(Bin[Bin$ASV==i,2])
  if(Bin.length>1){
    a<-Bin.length
    b<-Da[rownames(Da)==i,]/a
    for (j in Bin[Bin$ASV==i,2]){
      Single<-as.data.frame(b)
      Single$ID<-j
      All<-rbind(All,Single)
      Single<-tibble()
    }
  }
  else if(Bin.length==1){
    a<-1
    Single<-as.data.frame(Da[rownames(Da)==i,]/a)
    Single$ID<-Bin[Bin$ASV==i,2]
    All<-rbind(All,Single)
    Single<-tibble()
  }
  else if(Bin.length==0){
    ASV.missing<-as.data.frame(Da[rownames(Da)==i,])
    ASV.totmissing<-rbind(ASV.totmissing,ASV.missing)
    #cat(i,"is missing and has", sum(Da[i,]), "sequences in the dataset\n")
  }
}

ASV.totmissing<-as.data.frame(t(colSums(ASV.totmissing)))
ASV.totmissing$ID<-"Missing"
All1<-All %>% 
  rbind(ASV.totmissing) %>% 
  group_by(ID) %>%
  summarise_all(sum)
D3<-as.data.frame(t(All1[,-1]))
colnames(D3)<-All1$ID

cat("-------------------\n")
cat(" Save Genome Table\n")
cat("-------------------\n")
cat("\n")

if(VERBOSE==1){cat("Writing ASV table to:\n")}
write.csv(All1, paste(out.folder, "NEW.ASV.to.Genome.table.csv", sep=""))

BplotD4<-D3 %>% 
  rownames_to_column("X.SampleID") %>% 
  pivot_longer(!X.SampleID, names_to = "Genome", values_to = "value") %>% 
  left_join(M)

colourCount = length(unique(BplotD4$Genome))

Genome.barplot<-ggplot(BplotD4, aes(x=X.SampleID, y=value, fill=Genome)) +
  theme_bw() +
  geom_bar(stat = "identity", position="fill") +
  theme(axis.text.x = element_blank(),
        legend.position = "right") +
  facet_wrap(as.formula(paste("~",selecteur1)), scale="free", nrow=1) + 
  scale_fill_manual(values = getPalette(colourCount))+
  guides(fill=guide_legend(ncol=2)) +
  ggtitle(paste("Relative Genome abundance split by",selecteur1))

ggsave(paste(out.folder, "ASV.genome.png", sep=""), Genome.barplot, device = "png",width = 20, height = 10, dpi=150)

## ASV to genome table
cat("----------------\n")
cat("Pathways section\n")
cat("----------------\n")
cat("\n")

D5<-t(D3)
D5<-D5[!rownames(D5)=="Missing",]


#Making sure now row or column empty

null.row<-rownames(D5[rowSums(D5)==0,])
null.col<-colnames(D5[,colSums(D5)==0])

if(length(null.row)>0){
  cat(paste("The following genomes were present as ASV sequences but are not present in the count table:\n"))
  cat(null.row,"\n")
  cat(paste("They will be removed now\n"))
}

if(length(null.col)>0){
  cat(paste("The following samples are empty\n"))
  cat(null.col,"\n")
  cat(paste("They will be removed now\n"))
}

D5<-D5[!rownames(D5) %in% null.row,]
D5<-D5[,!colnames(D5) %in% null.col]

### Calculating how much of each sample we have
D4<-prop.table(as.matrix(D3),margin=1)
D4<-D4[,!colnames(D4)=="Missing"]
D4<-D4[rownames(D4) %in% colnames(D5),]
D4<-D4[,colnames(D4) %in% rownames(D5)]
D4<-rowSums(D4)

### PCoA
D6<-as.data.frame(t(D5)) %>% rownames_to_column("ID")
if(VERBOSE==1){cat("Plotting Genome PCoA\n")}
if(VERBOSE==1){cat("-----------------\n")}
selectors2<-c(colnames(D6)[1],selecteur1,selecteur1)
GeDist<-vegdist(D6[,!colnames(D6) %in% selectors2], method="jaccard")
GePco<-pco(GeDist, k=10)
Geig1<- round(GePco$eig[1]/sum(GePco$eig)*100,2)
Geig2<- round(GePco$eig[2]/sum(GePco$eig)*100,2)
GePCOplot<-as.data.frame(GePco$points)
GePCOplot$ID<-D6$ID
GePCOplot<-merge(GePCOplot,M, by.x="ID", by.y="X.SampleID", all.x=T)

colourCount = length(unique(GePCOplot[[selecteur1]]))
Genome.PCoa<-ggplot(GePCOplot, aes(x=V1, y=V2, col=get(selecteur1))) +
  geom_point() +
  theme_bw()+
  xlab(paste("Dimension 1", eig1, "%",sep=" "))+
  ylab(paste("Dimension 2", eig2, "%",sep=" "))+
  ggtitle("Principal Coordinates Analysis of Genome table")+
  scale_color_manual(values = getPalette(colourCount), name=paste(selecteur1))

ggsave(paste(out.folder, "PCOA.Genome.results.png", sep=""), Genome.PCoa, device = "png", width = 10, height = 10, dpi=150)

# Genome to Gene profiles

## Table generation profiles

PanGen<-WB %>%
  select(Genome, gene_cluster_id) %>%
  filter(gene_cluster_id!="NA") %>%
  mutate(gene_cluster_id=strsplit(gene_cluster_id,"\\|")) %>%
  unnest(gene_cluster_id) %>%
  group_by(Genome,gene_cluster_id) %>%
  summarise(panGene=n()) %>%
  pivot_wider(id_cols=gene_cluster_id,names_from = Genome,values_from = panGene, values_fill = 0)
Pan.db<-PanGen %>% column_to_rownames("gene_cluster_id") %>% t()
PGdb2<-Pan.db[rownames(Pan.db) %in% rownames(D5),]

###MOST IMPORTANT STEP OF THE SCRIPT
PG.results<-ASVtoFUNCTION(D5,PGdb2,out.folder,"Pangene")
if(VERBOSE==1){cat("Saving table to:\n")}
if(VERBOSE==1){cat(paste0(out.folder, "PanGene.profile.table.csv"),"\n")}
if(VERBOSE==1){cat("----------------\n")}
write.csv(x = PG.results,file = paste0(out.folder, "PanGene.profile.table.csv"))

# PCoA based on pan gene distribution

#PG.results=read.csv(paste0(out.folder, "PanGene.profile.table.csv"))
PG.results=PG.results %>% column_to_rownames("X")
## Generate PCoA plot of genomic composition
if(VERBOSE==1){cat("Plotting PCoA from Genome Data\n")}
if(VERBOSE==1){cat("----------------\n")}
PGRB<-PG.results
#K.results<-t(prop.table(as.matrix(K.results), margin = 2))
PG.results<-t(PG.results)*D4
PG.clr<-compositions::clr(PG.results)
PG.dist<-vegdist(PG.clr, method = "euclidian")
PGres2<-pco(PG.dist, k=2)
PGplot<-as.data.frame(PGres2$points)
PGplot$ID<-rownames(PGplot)
PGplot<-merge(PGplot,M[,colnames(M)], by.x="ID", by.y="X.SampleID", all.x=T)

PGeig1<- round(PGres2$eig[1]/sum(PGres2$eig)*100,2)
PGeig2<- round(PGres2$eig[2]/sum(PGres2$eig)*100,2)

PanGene.pcoa<-ggplot(PGplot, aes_string(x="V1", y="V2", col=selecteur1)) +
  geom_point()+
  theme_bw()+
  xlab(paste("Dimension 1", PGeig1, "%",sep=" "))+
  ylab(paste("Dimension 2", PGeig2, "%",sep=" "))+
  ggtitle("Principal Coordinates Analysis of Pangene Profiles Predictions")+
  scale_color_manual(values = getPalette(colourCount), name=paste(selecteur1))

ggsave(paste(out.folder,"PanGen.pcoa.png", sep=""), PanGene.pcoa, device = "png",width = 5, height = 5, dpi=150)

# Generate HTML report
if(VERBOSE==1){cat("Generating HTML report...\n")}
PWD <- getwd() # Save current working directory
rmarkdown::render(input = file.path(dirname(SourceDir), "Report.Rmd"),
                 output_file = file.path(out.folder, "WormBiome_Report.html"),
                 params = list(
                   PWD = PWD,
                   M = M,
                   selecteur1 = selecteur1,
                   Missing.sample.1 = Missing.sample.1,
                   Missing.sample.2 = Missing.sample.2,
                   count.threshold = count.threshold,
                   Sample.filteredout.read = Sample.filteredout.read,
                   Genome.present = Genome.present,
                   Genome.absent = Genome.absent,
                   null.row = null.row,
                   null.col = null.col,
                   SubRS = SubRS,
                   MeanASV = MeanASV,
                   ASV.barplot = ASV.barplot,
                   ASV.PCoa = ASV.PCoa,
                   out.folder = out.folder
                 ),
                 envir = environment())

if(VERBOSE==1){cat("Report generated successfully!\n")}

