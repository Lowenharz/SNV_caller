#! /illumina/scratch/K2-I/Users/qliu2/JIRA/Code/R-3.4.1/bin/Rscript --vanilla
VERSION = "1.0.0.1"
arg=commandArgs()
craft_run = unlist(strsplit(arg[ pmatch("--Folder",arg)], "="))[2]
outputFolder = unlist(strsplit(arg[ pmatch("--outputFolder",arg)], "="))[2]
cellline = unlist(strsplit(arg[ pmatch("--cell_lines",arg)],"="))[2]
cellline_list = unlist(strsplit(cellline,","))
cellline_list = cellline_list[ cellline_list != "" ]

if (length(cellline_list)==0) {
  stop("Please specific your cell line")
}


options(warn=-1)

## Check if file exists
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder)
}

## Load Libraries ## 
if (!require("stringr")) install.packages("stringr",repos = "http://cran.us.r-project.org")
library(stringr)

#### Start of the Code ### 
craft_run_list=list.files(craft_run)
folder_index=grep(pattern="fold.change",craft_run_list)

if (length(folder_index)==0){
  stop("Folder does not contain fold change files. Check again")
}
############ PRINT VERSION and script name ###########
print(paste("####### CNV Caller", VERSION, "#######"))
print("Loading CCLE CNV database")
CCLE_complete=read.table("/illumina/scratch/K2-I/software/SNV_auto_caller/database/CCLE_copynumber_byGene_2013-12-03.txt",header=T)



######################## Function to return cell line in CCLE table ##########################
lite_database_gen<-function(CCLE_name) {
  location=c()
  CCLE_header=read.table("/illumina/scratch/K2-I/software/SNV_auto_caller/resource/lite_CCLE_copynumber_byGene_2013-12-03.txt",header=F,stringsAsFactors = FALSE)
  for(k in 1:length(CCLE_header)){
    check_name = unlist(strsplit(as.character(CCLE_header[k]),"_"))[1]
    if (check_name==CCLE_name){
      location=cbind(k)
    }
  }
  return(location)
}

filename=craft_run_list[folder_index]
######################## Generate CCLE_index.txt ######################
i = 1
for(i in 1:length(filename)){
  print(paste("Processing ", filename[i], "..."))
  data_name = unlist(strsplit(filename[i],".",fixed=TRUE))[3]
  for(j in 1: length(cellline_list)) {
    sample_name=cellline_list[j]
    cell_line_name=toupper(unlist(strsplit(sample_name,"_",fixed=TRUE))[1])
    CCLE_compatible_name=gsub("[^[:alnum:]]", "", cell_line_name)
    
    CCLE_position=lite_database_gen(CCLE_compatible_name)
    if(length(CCLE_position)>1){ 
      stop(paste("There is ambiguity in cell line",CCLE_compatible_name,"at position", CCLE_position))
    }
    
    if (length(CCLE_position)!=0){
      mapping_index=cbind(CCLE_compatible_name, CCLE_position)
      # map=rbind(map,mapping_index)
      CCLE_index=as.numeric(mapping_index[1,2])
      lite=CCLE_complete[,cbind(2,CCLE_index)]
      lite$CCLE_log_2_power=with(lite,2^lite[,2])
      orgin=read.table(paste(craft_run,filename[i],sep=''),header=T,stringsAsFactors=F)
      # orgin$log_2_power=with(orgin,2^orgin[,2])
      match_table=c()
      for(k in 1:nrow(orgin)){
        gene=orgin[k,1]
        matched_gene=lite[grep(paste("^",gene,"$", sep=""),lite$SYMBOL),]
        if (length(matched_gene) > 0) {
          match_table=rbind(match_table,matched_gene)
        }
      }
      combined_table=merge(orgin,match_table,by.x=c("gene"),by.y=c("SYMBOL"))
      print(paste("Writing CCLE_combined_", filename[i], ' with ',sample_name, " cell_line", " ...",sep=""))
      colnames(combined_table)= c("gene","fold_change", "t_stat", "q_score", "num_targets", "HS675T_LARGE_INTESTINE","CCLE_fold_change")
      write.table(combined_table,file=paste(outputFolder,"/CCLE_combined_",data_name,"_",sample_name,".txt",sep=""),
                  col.names = TRUE,row.names = FALSE,quote = FALSE)
    } else { 
      print(paste(cell_line_name,"does not exist in CCLE Database"))
    }
  } ########### Loop thourgh each cell line ###########  
}
