#Function ASV to gene profile
ASVtoFUNCTION <- function(ASV.table, Function.table, out.folder, pref) {
  if(VERBOSE==1){
    cat("Processing", ncol(ASV.table), "samples\n")
    cat("ASV table dimensions:", dim(ASV.table)[1], "x", dim(ASV.table)[2], "\n")
    cat("Function table dimensions:", dim(Function.table)[1], "x", dim(Function.table)[2], "\n")
  }
  
  # Convert to matrices and ensure proper orientation
  ASV.matrix <- as.matrix(ASV.table)
  Function.matrix <- as.matrix(Function.table)
  
  # Pre-calculate the multiplication once
  # Function.matrix should be genes x genomes, ASV.matrix should be genomes x samples
  if(VERBOSE==1) cat("Calculating function profiles...\n")
  result <- t(Function.matrix) %*% ASV.matrix
  
  # Convert to data frame and handle column names
  PathwayList <- as.data.frame(result)
  colnames(PathwayList) <- colnames(ASV.table)
  rownames(PathwayList) <- colnames(Function.matrix)
  
  # Filter zero rows and columns
  P2 <- PathwayList[rowSums(PathwayList) > 0, , drop = FALSE]
  P2 <- P2[, colSums(P2) > 0, drop = FALSE]
  
  if(VERBOSE==1){
    cat("Final matrix dimensions:", dim(P2)[1], "x", dim(P2)[2], "\n")
    cat("Writing", pref, "table to:", paste0(out.folder, pref, ".profile.table.csv"), "\n")
  }
  write.csv(P2, paste0(out.folder, pref, ".profile.table.csv"), row.names = TRUE)
  
  return(P2)
}

Deseq.wrapper<-function(RTable,MetaT,selector){  
  rownames(MetaT)<-MetaT$X.SampleID
  phy<-phyloseq::otu_table(round(RTable), taxa_are_rows = T)
  MM<-phyloseq::sample_data(MetaT)
  DES<-phyloseq::phyloseq(phy,MM)
  phydeseqB<-suppressMessages(phyloseq::phyloseq_to_deseq2(physeq = DES, eval(parse(text = paste("~", selector)))))
  DEphy<-suppressMessages(DESeq2::DESeq(phydeseqB))
  return(DEphy)
}

Extract.results<-function(deseqobj, MetaT, type){
  Deseq.result<-list()
  DS<-list()
  j=0
  for (i in 1:ncol(comb.comp)){
    j=j+1
    a=comb.comp[1,i]
    b=comb.comp[2,i]
    if(VERBOSE==1){cat("Comparing ",a," vs ", b, "\n")}
    Deseq.result[[j]]<-suppressMessages(DESeq2::results(deseqobj, alpha=0.01, contrast=c("Selecteur", a, b)))
    DS[[j]] <- Deseq.result[[j]][which(Deseq.result[[j]]$padj < 0.01), ]
    DS[[j]] <- as(DS[[j]], "data.frame")
    DS[[j]]$sel <- rownames(DS[[j]])
    colnames(DS[[j]]) <- gsub("sel",type,colnames(DS[[j]]))
    #colnames(DS[[j]]) <- paste("C_",a,"vs",b,"_",colnames(DS[[j]]), sep="")
    colnames(DS[[j]])[7]<-type
    if(nrow(DS[[j]])>0){
    DS[[j]]$Comparison<-paste(a,b,sep=" - ")
    }
    if(j==1){
      Full.data<-DS[[j]]
    }
    
    if(j>1 & nrow(DS[[j]])>0){
      Full.data<- rbind(Full.data,DS[[j]])
    }
  }
  Full.data =Full.data %>%
    inner_join(MetaT, by=type)
  return(Full.data)
}


Gameover<-function(x){                                                              
  message(paste(" ____                                            "))
  message(paste("/ ___|  __ _ _ __ __    ____    ____     _____ _ __ "))
  message(paste("| |  _ / _` | '_ ` _ \\ / _  \\  / __ \\ \\ / / _ \\ '__|"))
  message(paste("| |_| | (_| | | | | | |  ___/ | |__| \\ V /  __/ |   "))
  message(paste("\\_____|\\__,_|_| |_| |_|\\___|   \\____/ \\_/ \\___|_|   "))
}
