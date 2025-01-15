#Function ASV to gene profile
ASVtoFUNCTION<-function(ASV.table, Function.table, out.folder,pref){
  Single<-tibble()
  All<-list()
  Function.table<-as.data.frame(Function.table)
  ASV.col<-ncol(ASV.table)
  cur.col<-0
  for (sample in colnames(ASV.table)){
    cur.col<-cur.col+1
    cat("Sample: ", cur.col, "/",ASV.col,"\n")
    if(VERBOSE==1){cat(paste("Processing sample", sample,"\n"))}
    Single=Function.table*ASV.table[,colnames(ASV.table)==sample]
    tmp<-as.data.frame(colSums(Single))
    Single<-tibble()
    colnames(tmp)<-sample
    All<-c(All,list(tmp))
    if(VERBOSE==1){cat("---\n\n")
      cat("")}
  }
  PathwayList<-as.data.frame(All)
  PathwayList[is.na(PathwayList)]<-0
  P2<-PathwayList[!rowSums(PathwayList)==0,]
  P2<-P2[,!colSums(P2)==0]
  
  if(VERBOSE==1){cat("Writing ",pref," table to:")}
  write.csv(P2, paste(out.folder, pref,".profile.table.csv", sep=""))
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
