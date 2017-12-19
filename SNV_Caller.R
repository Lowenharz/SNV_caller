#! /illumina/scratch/K2-I/Users/qliu2/JIRA/Code/R-3.4.1/bin/Rscript --vanilla
## Load Libraries ## 
VERSION = "1.0.0.1"
if (!require("stringr")) install.packages("stringr",repos = "http://cran.us.r-project.org")
library(stringr)

arg=commandArgs()
vcf_folder = unlist(strsplit(arg[ pmatch("--Folder",arg)], "="))[2]
output = unlist(strsplit(arg[ pmatch("--output",arg)], "="))[2]
cellline = unlist(str_split(arg[ pmatch("--cell_lines",arg)],"="))[2]
agnostic = unlist(str_split(arg[ pmatch("--agnostic",arg)],"="))[2]
cellline_list = unlist(str_split(cellline,","))
cellline_list = cellline_list[ cellline_list != "" ]

### If nothing is passed in celline default as agnostic, even input differently ###
if (length(cellline_list)==0) {
  cellline_list = ""
  agnostic = "TRUE"
}

if (startsWith(toupper(agnostic),"T")) {
  agnostic=TRUE
  print("Agnostic mode is on")
} else {
  agnostic=FALSE
  print("Agnostic mode is off")
}



## Check if file exists
if (file.exists(output)) {
  stop("Output file exists. Please use a new name or a different path.")
}

## Load Libraries ## 
library(stringr)


# Get list of files in folder
folder_list=list.files(vcf_folder)

folder_index=grep(pattern=".vcf",folder_list)
if (length(folder_index)==0){
  stop("Folder does not contain vcf files. Check again")
}
vcf_files=folder_list[folder_index]

############ PRINT VERSION and script name ###########
print(paste("####### SNV Caller", VERSION,"#######"))

### Read in CCLE and COSMIC SNP files
CCLE_table=read.csv("/illumina/scratch/K2-I/software/SNV_auto_caller/resource/CCLE_SNV_indel.csv",sep=",",stringsAsFactors = FALSE)
COSMIC_table=read.csv("/illumina/scratch/K2-I/software/SNV_auto_caller/resource/COSMIC_SNV_Indel.csv",sep=",",stringsAsFactors = FALSE)



startsWith <- function(x,prefix) {
  return(substring(x, 1, nchar(prefix)) == prefix)
}

endsWith <- function(x,suffix) {
  return(substring(x, (nchar(x)-nchar(suffix)+1), nchar(x)) == suffix)
}

CCLE_empty=FALSE
COSMIC_empty=FALSE
title=cbind("Sample_name","Cell_Line","Gene","Chr","Start",	"End",	"Variant_Type","Ref","Alt",	"DP","Observed variant frequency (1=100%)","Data Source", "CLLE VF if exist")
write.table(title,file=output, sep=",", col.names = FALSE,row.names = FALSE,append = TRUE)

for(j in 1:length(vcf_files)){
  vcf=vcf_files[j]
  print(paste("Processing",vcf, "..."))
  
  ### Can't make cell line agnostic yet because some cell lines have variants that others don't
  CCLE_rows=c()
  COSMIC_rows=c()
  filename=basename(vcf)
  sample_name=unlist(strsplit(filename,".vcf",fixed=TRUE))[1]
  ### Take the part with cellline name in CCLE table
  for (k in 1:length(cellline_list)) {
    cell_line_name = toupper(cellline_list[k])
    # print(cell_line_name)
    CCLE_compatible_name=gsub("[^[:alnum:]]", "", cell_line_name)
    for(i in 1:nrow(CCLE_table)) {
      check_name = unlist(strsplit(CCLE_table[i,1],"_"))[1]
      # print(check_name)
      # print(check_name==CCLE_compatible_name)
      if (check_name==CCLE_compatible_name){
        CCLE_rows=cbind(CCLE_rows,i)
      }
    }
    ### Take the part with cellline name in COSMIC table
    for(i in 1:nrow(COSMIC_table)) {
      if (COSMIC_table[i,1]==cell_line_name){
        COSMIC_rows=cbind(COSMIC_rows,i)
      }
    }
  }
  
  ## Create cellline table for CCLE and 
  #### If agnositic CCLE_cellline_only_table = CCLE_table ####
  ####      COSMIC_cellline_only_table = COSMIC_table     ####
  if (agnostic) { 
    CCLE_cellline_only_table=CCLE_table
    COSMIC_cellline_only_table=COSMIC_table
  } else {
    CCLE_cellline_only_table=CCLE_table[CCLE_rows,]
    COSMIC_cellline_only_table=COSMIC_table[COSMIC_rows,]
  }
  
  ### Check to see if table is empty ###
  if(nrow(CCLE_cellline_only_table)==0) { 
    CCLE_empty=TRUE
  } else {
    CCLE_empty=FALSE
  }
  
  if(nrow(COSMIC_cellline_only_table)==0) { 
    COSMIC_empty=TRUE
  } else {
    COSMIC_empty=FALSE
  }
  
  ### Get the start position of variant according to CCLE and COSMIC
  ### could make the list shorter if ignore same variant from diffrent
  ### cell line
  CCLE_start = CCLE_cellline_only_table$Start
  COSMIC_start = COSMIC_cellline_only_table$Start
  
  ### Read vcf file
  vcf.data=read.table(paste(vcf_folder,vcf,sep=''),stringsAsFactors=FALSE)
  vcf.data.POS=vcf.data[,2]
  
  ### Function for checking if the returned integer is empty or not
  is.integer0 <- function(x)
  {
    is.integer(x) && length(x) == 0L
  }
  
  
  if (!CCLE_empty){
    ################ Look through the CCLE data base
    CCLE_matched_list = c()
    CCLE_af_list=c()
    CCLE_DP_list=c()
    CCLE_EAF_list=c()
    for(i in 1:nrow(CCLE_cellline_only_table)) {
      CCLE_pos=CCLE_start[i]
      vcf.index = which(vcf.data.POS==CCLE_pos)
      vcf.index= vcf.index[1] # make sure getting the first matching position 
      if (!is.na(vcf.index)){
        variant = vcf.data[vcf.index,]
        vcf.result = variant[,10]
        vcf.dp = as.numeric(unlist(strsplit(variant[,8],"="))[2])
        vcf.af = as.numeric(unlist(strsplit(vcf.result,":"))[5])
        EAF = as.numeric(CCLE_table$Variant_Frequency)
        ## test if AF is greater than zero in VCF
        if(vcf.af > 0) {
          ## test if Ref==Ref and Alt==Tumor SNP
          if(CCLE_cellline_only_table[i,8] == variant[,4] &
             CCLE_cellline_only_table[i,9] == variant[,5]
          ) {
            CCLE_matched_list=cbind(CCLE_matched_list,i)
            CCLE_af_list =rbind(CCLE_af_list,vcf.af)
            CCLE_DP_list = rbind(CCLE_DP_list,vcf.dp)
            CCLE_EAF_list = rbind(CCLE_EAF_list, EAF)
          }
          ## test for DEL
          else if(endsWith(vcf.data[vcf.index-1,4],CCLE_cellline_only_table[i,8])) {
            CCLE_matched_list=cbind(CCLE_matched_list,i)
            CCLE_af_list =rbind(CCLE_af_list,vcf.af)
            CCLE_DP_list = rbind(CCLE_DP_list,vcf.dp)
            CCLE_EAF_list = rbind(CCLE_EAF_list, EAF)
          }
          ## test for INS
          else if(CCLE_cellline_only_table[i,8] == '-' &
                  str_length(vcf.data[vcf.index,4]) < str_length(vcf.data[vcf.index,5])) {
            CCLE_matched_list=cbind(CCLE_matched_list,i)
            CCLE_af_list =rbind(CCLE_af_list,vcf.af)
            CCLE_DP_list = rbind(CCLE_DP_list,vcf.dp)
            CCLE_EAF_list = rbind(CCLE_EAF_list, EAF)
          }
          ## test for DNP
          else if(CCLE_cellline_only_table[i,8] == paste(variant[,4],vcf.data[vcf.index+1,4],sep='') &
                  CCLE_cellline_only_table[i,9] == paste(variant[,5],vcf.data[vcf.index+1,5],sep='')) {
            CCLE_matched_list=cbind(CCLE_matched_list,i)
            CCLE_af_list =rbind(CCLE_af_list,vcf.af)
            CCLE_DP_list = rbind(CCLE_DP_list,vcf.dp)
            CCLE_EAF_list = rbind(CCLE_EAF_list, EAF)
          }
        }
      }
    }
    
    ################## Making CCLE Report Table ####################
    CCLE_columns=c(1:6,8,9,16)
    report_table_CCLE=CCLE_cellline_only_table[CCLE_matched_list,CCLE_columns]
    if (nrow(report_table_CCLE)==0) {
      CCLE_empty = TRUE
    } else {
      report_table_CCLE = cbind(report_table_CCLE,CCLE_af_list, CCLE_DP_list,CCLE_EAF_list)
      report_table_CCLE = report_table_CCLE[,c(1:8,11,10,9,12)] # Reorganize
    }
  }
  
  if(!COSMIC_empty) { 
    ############# Look through the COSMIC data base
    COSMIC_matched_list=c()
    COSMIC_af_list=c()
    COSMIC_type_list=c()
    COSMIC_DP_list = c()
    for(i in 1:nrow(COSMIC_cellline_only_table)) {
      COSMIC_pos=COSMIC_start[i]
      vcf.index = which(vcf.data.POS==COSMIC_pos)
      vcf.index=vcf.index[1]
      if (!is.na(vcf.index)){
        variant = vcf.data[vcf.index,]
        vcf.result = variant[,10]
        vcf.af = as.numeric(unlist(strsplit(vcf.result,":"))[5])
        vcf.dp = as.numeric(unlist(strsplit(variant[,8],"="))[2])
        ## test if AF is greater than zero in VCF
        if(vcf.af > 0) {
          ## test if Ref== Ref and Alt==Tumor
          if(COSMIC_cellline_only_table[i,6] == variant[,4] &
             COSMIC_cellline_only_table[i,7] == variant[,5]) {
            COSMIC_matched_list=cbind(COSMIC_matched_list,i)
            COSMIC_af_list=rbind(COSMIC_af_list,vcf.af)
            COSMIC_DP_list = rbind(COSMIC_DP_list,vcf.dp)
            if(str_length(COSMIC_cellline_only_table[i,6]) > str_length(COSMIC_cellline_only_table[i,7])) {
              COSMIC_type_list = rbind(COSMIC_type_list,"DEL")
            } else if (str_length(COSMIC_cellline_only_table[i,6]) < str_length(COSMIC_cellline_only_table[i,7])) {
              COSMIC_type_list = rbind(COSMIC_type_list,"INS")
            } else if (str_length(COSMIC_cellline_only_table[i,6])>1){
              COSMIC_type_list = rbind(COSMIC_type_list,"DNS")
            } else {
              COSMIC_type_list = rbind(COSMIC_type_list,"SNP")
            }
          }
        }
      }
    }
    
    
    ################## Making COSMIC Report Table ####################
    COSMIC_columns=c(1:8)
    report_table_COSMIC=COSMIC_cellline_only_table[COSMIC_matched_list,COSMIC_columns]
    report_table_COSMIC=cbind(report_table_COSMIC, COSMIC_af_list, COSMIC_DP_list)
    report_table_COSMIC=cbind(report_table_COSMIC, COSMIC_type_list, NA)
    if (nrow(report_table_COSMIC)==0) {
      COSMIC_empty = TRUE
    } else {
      report_table_COSMIC=report_table_COSMIC[,c(1:5,11,6:7,10,9,8,12)]
      colnames(report_table_COSMIC)=colnames(report_table_CCLE)
    }
  }
  
  ################# Eliminate repeats but keep the differences ###################
  
  if (!COSMIC_empty && !CCLE_empty) {
    index_delete_list=c()
    source_change_list=c()
    for(i in 1:nrow(report_table_CCLE)){
      gene=report_table_CCLE[i,2]
      index_delete=which(report_table_COSMIC[,2]==gene)
      if(!is.integer0(index_delete)){
        if(report_table_COSMIC[index_delete,4] == report_table_CCLE[i,4]&
           report_table_COSMIC[index_delete,5] == report_table_CCLE[i,5]){
          index_delete_list= cbind(index_delete_list,index_delete)
          source_change_list= cbind(source_change_list,i)
        }
      }
    }
    report_table_CCLE[source_change_list, 11] = "CCLE & COSMIC"
    if ( !is.null(index_delete_list)) {
      report_table_COSMIC= report_table_COSMIC[-index_delete_list,]  
    }
    combined_list=rbind(report_table_CCLE,report_table_COSMIC)
    combined_list$Sample_name= sample_name
    combined_list=combined_list[c(length(combined_list),1:(length(combined_list)-1))]
    # print(combined_list)
    write.table(combined_list,file=paste(output),sep=",", col.names = FALSE,row.names = FALSE,append = TRUE)
  } else if (!CCLE_empty){
    # print(report_table_CCLE)
    report_table_CCLE$Sample_name= sample_name
    report_table_CCLE=report_table_CCLE[c((length(report_table_CCLE)-2),1:(length(report_table_CCLE)-3),length(report_table_CCLE))]
    write.table(report_table_CCLE,file=paste(output),sep=",", col.names = FALSE,row.names = FALSE,append = TRUE)
  } else if (!COSMIC_empty) {
    # print(report_table_COSMIC)
    report_table_COSMIC$Sample_name=sample_name
    report_table_COSMIC=report_table_COSMIC[c(length(report_table_COSMIC),1:(length(report_table_COSMIC)-1))]
    write.table(report_table_COSMIC,file=paste(output),sep=",", col.names = FALSE,row.names = FALSE,append = TRUE)
  } else {
    print(paste("Both CCLE and COSMIC does not contain information of ",sample_name))
  }
}