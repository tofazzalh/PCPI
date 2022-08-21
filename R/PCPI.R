#' Prediction of interaction between circRNA and protein
#'
#' This function takes the circRNA and protein sequences as input and predicts
#' the interaction between them.
#' @param circRNA_seq A dataframe containing the names and sequences of circRNAs
#' @param protein_seq A dataframe containing the names and sequences of proteins
#' @param output_folder The name of the output folder
#'
#' @return The interaction between the circRNA and protein will written in
#' output_folder
#'
#' @importFrom utils write.table
#' @importFrom stats predict
#' @importFrom Biostrings readDNAStringSet
#' @importFrom seqinr read.fasta
#' @importFrom e1071 svm
#' @export
#'
#' @examples
#' #Loading an example data for circRNA names and sequences.
#' a<-data("circRNA_seq_df")
#' circ_df<-circRNA_seq_df
#'
#' #Loading an example data for protein names and sequences.
#' b<-data("protein_seq_df")
#' protein_df<-protein_seq_df
#'
#' #predicting circRNA and protein interaction. Here, the output will be
#' #written in file interaction_data.txt in out_dir directory
#' out_dir<-tempdir()
#' circ_protein_interaction(circ_df, protein_df, out_dir)
#'
#'
#'
circ_protein_interaction<-function(circRNA_seq, protein_seq, output_folder){
  #Read circRNA sequence
  circ_seq_name = circRNA_seq[,1]
  circ_sequence = circRNA_seq[,2]
  #Read protein sequence
  df_protein<-protein_seq
  row.names(df_protein)<-1:dim(df_protein)[1]
  colnames(df_protein)<-c("protein_name", "protein_seq")
  #All possible combination of feature 1
  P_wi_1 <- expand.grid(rep(list(1:8), 3))
  F_1<-paste0(P_wi_1$Var1, P_wi_1$Var2, P_wi_1$Var3, collaspe="")
  #All possible combination of feature 2
  P_wi_2 <- expand.grid(rep(list(1:7), 3))
  F_2<-paste0(P_wi_2$Var1, P_wi_2$Var2, P_wi_2$Var3, collaspe="")
  # start of function for counting k-mer
  k_mer_count<-function(s,k){
    kmer=NULL
    for(i in 1:nchar(s)-k+1){
      kmer[i]=substr(s,i,i+k-1)
    }
    kmer
  }
  # end of function for counting k-mer

  ### start of generating feature from circRNA
  Feature1_value=NULL
  #a<-load(file='F:/Third_project/Third_project/index_zero_mean_N.rda')
  Feature_1_colmean_zero_index<-index_zero_mean_N[which(index_zero_mean_N<=512)]
  for(k in 1:length(circ_seq_name)){
    d<-circ_sequence[k] #assign kth circRNA seq
    D<-toupper(d) #convert the nuecliotides into upper case
    D_R<-chartr("T","U",D) # convert "T" into "U"
    kmer<-k_mer_count(D_R,3) # counting 3-mers of the circRNA seq
    # start of codon transfer: from RNA to protein
    for(i in 1:length(kmer)){
      if(kmer[i]=="UUU"|kmer[i]=="UUC")kmer[i]="F"
      else if (kmer[i]=="UUA"|kmer[i]=="UUG"|kmer[i]=="CUU"|kmer[i]=="CUC"|kmer[i]=="CUA"|kmer[i]=="CUG") kmer[i]="L"
      else if (kmer[i]=="AUU"|kmer[i]=="AUC"|kmer[i]=="AUA")kmer[i]="I"
      else if (kmer[i]=="AUG")kmer[i]="M"
      else if (kmer[i]=="GUU"|kmer[i]=="GUC"|kmer[i]=="GUA"|kmer[i]=="GUG")kmer[i]="V"
      else if (kmer[i]=="UCU"|kmer[i]=="UCC"|kmer[i]=="UCA"|kmer[i]=="UCG"|kmer[i]=="AGU"|kmer[i]=="AGC")kmer[i]="S"
      else if (kmer[i]=="CCU"|kmer[i]=="CCC"|kmer[i]=="CCA"|kmer[i]=="CCG")kmer[i]="P"
      else if (kmer[i]=="ACU"|kmer[i]=="ACC"|kmer[i]=="ACA"|kmer[i]=="ACG")kmer[i]="T"
      else if (kmer[i]=="GCU"|kmer[i]=="GCC"|kmer[i]=="GCA"|kmer[i]=="GCG")kmer[i]="A"
      else if (kmer[i]=="UAU"|kmer[i]=="UAC")kmer[i]="Y"
      else if (kmer[i]=="UAA"|kmer[i]=="UAG"|kmer[i]=="UGA")kmer[i]="Z"
      else if (kmer[i]=="CAU"|kmer[i]=="CAC")kmer[i]="H"
      else if (kmer[i]=="CAA"|kmer[i]=="CAG")kmer[i]="Q"
      else if (kmer[i]=="AAU"|kmer[i]=="AAC")kmer[i]="N"
      else if (kmer[i]=="AAA"|kmer[i]=="AAG")kmer[i]="K"
      else if (kmer[i]=="GAU"|kmer[i]=="GAC")kmer[i]="D"
      else if (kmer[i]=="GAA"|kmer[i]=="GAG")kmer[i]="E"
      else if (kmer[i]=="UGU"|kmer[i]=="UGC")kmer[i]="C"
      else if (kmer[i]=="UGG")kmer[i]="W"
      else if (kmer[i]=="CGU"|kmer[i]=="CGC"|kmer[i]=="CGA"|kmer[i]=="CGG"|kmer[i]=="AGA"|kmer[i]=="AGG")kmer[i]="R"
      else kmer[i]="G"
    }
    Pro<-paste(kmer,collapse="")
    Pro<-paste(kmer,collapse="")
    # end of codon transfer: from RNA to protein

    # start of converting protein seq into 8 letter characters
    Pro_8=NULL
    for (j in 1: nchar(as.character(Pro))){
      if (substring(Pro, j,j)=="A"|substring(Pro, j,j)=="G"|substring(Pro, j,j)=="V")Pro_8[j]="1"
      else if (substring(Pro, j,j)=="I"|substring(Pro, j,j)=="L"|substring(Pro, j,j)=="F"|substring(Pro, j,j)=="P")Pro_8[j]="2"
      else if (substring(Pro, j,j)=="Y"|substring(Pro, j,j)=="M"|substring(Pro, j,j)=="T"|substring(Pro, j,j)=="S")Pro_8[j]="3"
      else if (substring(Pro, j,j)=="H"|substring(Pro, j,j)=="N"|substring(Pro, j,j)=="Q"|substring(Pro, j,j)=="W")Pro_8[j]="4"
      else if (substring(Pro, j,j)=="R"|substring(Pro, j,j)=="K")Pro_8[j]="5"
      else if (substring(Pro, j,j)=="D"|substring(Pro, j,j)=="E")Pro_8[j]="6"
      else if (substring(Pro, j,j)=="C")Pro_8[j]="7"
      else Pro_8[j]="8"
    }
    Pro_8<-paste(Pro_8, collapse="")
    # end of converting protein seq into 8 letter characters

    Pro_8_kmer<-k_mer_count(Pro_8,3) # counting 3-mers of the 8 letter protein seq
    A<-table(Pro_8_kmer)# counting the frequency of each 3-mer

    # start of counting circRNA feature value
    F_1v=rep(0,512)
    for( m in 1:dim(A)){
      v<-which(F_1==names(A)[m])
      F_1v[v]=A[[m]]
    }
    # End of counting circRNA feature value
    Feature1_value<-rbind(Feature1_value, F_1v) # accumulating circRNA feature value
  }
  Feature1_value_modified<-Feature1_value[,-Feature_1_colmean_zero_index]
  ### end of generating feature from circRNA

  #### start of generating feature from protein
  Feature2_value=NULL
  Feature_2_colmean_zero_index<-index_zero_mean_N[which(index_zero_mean_N>512)]-512
  for(k in 1:dim(df_protein)[1]){
    d<-as.character(df_protein[k,2]) #assign kth protein seq
    Pro<-d
    # start of converting protein seq into 7 letter characters
    Pro_7=NULL
    for (j in 1: nchar(as.character(Pro))){
      if (substring(Pro, j,j)=="A"|substring(Pro, j,j)=="G"|substring(Pro, j,j)=="V")Pro_7[j]="1"
      else if (substring(Pro, j,j)=="I"|substring(Pro, j,j)=="L"|substring(Pro, j,j)=="F"|substring(Pro, j,j)=="P")Pro_7[j]="2"
      else if (substring(Pro, j,j)=="Y"|substring(Pro, j,j)=="M"|substring(Pro, j,j)=="T"|substring(Pro, j,j)=="S")Pro_7[j]="3"
      else if (substring(Pro, j,j)=="H"|substring(Pro, j,j)=="N"|substring(Pro, j,j)=="Q"|substring(Pro, j,j)=="W")Pro_7[j]="4"
      else if (substring(Pro, j,j)=="R"|substring(Pro, j,j)=="K")Pro_7[j]="5"
      else if (substring(Pro, j,j)=="D"|substring(Pro, j,j)=="E")Pro_7[j]="6"
      else Pro_7[j]="7"
    }
    Pro_7<-paste(Pro_7, collapse="")
    # end of converting protein seq into 7 letter characters

    Pro_7_kmer<-k_mer_count(Pro_7,3) # counting 3-mers of the 7 letter protein seq
    B<-table(Pro_7_kmer)# counting the frequency of each 3-mer

    # start of counting protein feature value
    F_2v=rep(0,343)
    for( n in 1:dim(B)){
      v<-which(F_2==names(B)[n])
      F_2v[v]=B[[n]]
    }
    # end of counting protein feature value
    Feature2_value<-rbind(Feature2_value, F_2v) # accumulating protein feature value
  }
  Feature2_value_modified<-Feature2_value[,-Feature_2_colmean_zero_index]
  #### end of generating feature from protein

  ######claculating interaction
  interaction<-data.frame(circRNA=character(), protein=character(), interaction_probability=numeric(), interaction_status=character())
  write.table(interaction, file.path(output_folder,"interaction_data.txt"),sep="\t", row.names=F, quote=F)
  #my_model<-readRDS(file='F:/Third_project/Third_project/Example_data/my_model_final.rda')
  for(circRNA in 1:length(circ_seq_name)){
    for (protein in 1:dim(df_protein)[1]){
      data<-c(Feature1_value_modified[circRNA,],Feature2_value_modified[protein,])
      data_1<-data.frame(t(data))
      colnames(data_1)<-paste("V",1:500,sep="")
      pred<- predict(my_model,data_1, probability=T)
      probability<-attr(pred, 'probabilities')[1]
      if (probability>=0.5)Int_status="Yes" else Int_status="No"
      Int<-data.frame(circ_seq_name[circRNA], df_protein[protein,1], probability, Int_status)
      write.table(Int, file.path(output_folder,"interaction_data.txt"),sep="\t", row.names=F, col.names=F, quote=F, append=T)
    }
  }
}
