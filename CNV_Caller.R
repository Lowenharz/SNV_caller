#! /illumina/scratch/K2-I/Users/qliu2/JIRA/Code/R-3.4.1/bin/Rscript --vanilla

arg=commandArgs()
craft_run = unlist(strsplit(arg[ pmatch("--Folder",arg)], "="))[2]
outputFolder = unlist(strsplit(arg[ pmatch("--outputFolder",arg)], "="))[2]
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
folder_index=grep(pattern="foldChange",craft_run_list)

if (length(folder_index)==0){
  stop("Folder does not contain fold change files. Check again")
}

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
for(i in 1:length(filename)){
  sample_name=unlist(strsplit(filename[i],".",fixed=TRUE))[1]
  cell_line_name=toupper(unlist(strsplit(sample_name,"_",fixed=TRUE))[1])
  CCLE_compatible_name=gsub("[^[:alnum:]]", "", cell_line_name)
  print(paste("Processing ", filename[i], "..."))
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
    orgin$log_2_power=with(orgin,2^orgin[,2])
    match_table=c()
    for(j in 1:nrow(orgin)){
      gene=orgin[j,1]
      matched_gene=lite[grep(paste("^",gene,"$", sep=""),lite$SYMBOL),]
      if (length(matched_gene) > 0) {
        match_table=rbind(match_table,matched_gene)
      }
    }
    combined_table=merge(orgin,match_table,by.x=c("GeneName"),by.y=c("SYMBOL"))
    print(paste("Writing CCLE_combined_", filename[i], " ...",sep=""))
    write.table(combined_table,file=paste(outputFolder,"/CCLE_combined_",sample_name,".txt",sep=""),
                col.names = TRUE,row.names = FALSE,quote = FALSE)
  } else { 
    print(paste(cell_line_name,"does not exist in CCLE Database"))
  }
  
}
