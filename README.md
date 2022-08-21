
# PCPI: an R package for predicting circRNA and protein interaction

## Requirements
* R (>= 4.2.0)
* Biostrings
* seqinr
* e1071

## Installation
### From github
To install the package from github first you need to install the package “devtools” using the following command:

    install.packages("devtools", dep=T)

The package "PCPI" depends on one bioconductor packages "Biostrings" which cannot be installed automatically while installing "PCPI" using "devtools". So, you need to install this packages manually using the following way:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("Biostrings")

Finally, install “circRNAFull” by the following command:

    devtools::install_github("tofazzalh/PCPI", dep = T)

Start analysis by loading the package with the following command:

    library("PCPI")

## Prediction of interaction between circRNA and protein

#### Description
This function takes the circRNA and protein sequences as input and predicts the interaction between them.

#### Usage
    
    circ_protein_interaction(circRNA_seq, protein_seq, output_folder)
    
#### Arguments

`circRNA_seq` A dataframe containing the names and sequences of circRNAs.

`protein_seq` A dataframe containing the names and sequences of proteins.

`output_folder` The name of the output folder.

#### Value
The interaction between the circRNA and protein will written in *output_folder*

#### Example
    
    #Loading an example data for circRNA names and sequences.
    a<-data("circRNA_seq_df")
    circ_df<-circRNA_seq_df

    #Loading an example data for protein names and sequences.
    b<-data("protein_seq_df")
    protein_df<-protein_seq_df

    #predicting circRNA and protein interaction. Here, the output 
    #will be written in file interaction_data.txt in out_dir directory
    out_dir<-tempdir()
    circ_protein_interaction(circ_df, protein_df, out_dir)
