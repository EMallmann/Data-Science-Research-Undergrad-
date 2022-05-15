library(rentrez);
library(readr);
snp_names <- read.csv("C:/Users/ericm/NewFolder2/PLINK_GRU/fullforest2.csv", header = FALSE, na.strings = c("NA",9), fill = TRUE);
snp_names = as.vector(snp_names[36950:37322,1]); #top 1% = 3:375 &&&&& bottom 1% = 36950:37322


##################################   TITLE SECTION   ###################################
titles <-data.frame(SNP=character(),Title=character(),stringsAsFactors = FALSE);
for (q in snp_names){
  ids_for_snp = entrez_search(db = "pubmed", term= q, use_history = TRUE); #gives the unique record ids associated with this SNP. 
  num_publications = as.integer(ids_for_snp[2]);
  if (num_publications==0){
    next();
  }
  summary = entrez_summary(db = "pubmed", web_history = ids_for_snp$web_history) #gives us the summary of publications associated with this SNP in the Pubmed database
  i = 1;
  
  while(i<=length(summary)){
    if (num_publications==1){
      title = summary$title;
    }
    else{
      title = summary[[i]]$title;
    }
    titles[dim(titles)[1]+1,] = c(q ,title);
    if(num_publications ==1){
      break;
    }
    i = i+1;
  }
}



##################################### SNP INFO SECTION ##################################################
snp_info <-data.frame(SNP=character(),Fxn_Class=character(), Clin_Sign=character(), stringsAsFactors = FALSE);

i = 1;
for(q in snp_names){
  ids_for_snp = entrez_search(db = "snp", term= q , use_history = TRUE); #gives the unique record ids associated with this SNP. 
  num_ids = as.integer(ids_for_snp[2]);
  
  if (num_ids == 0){
    snp_info[i,] = list(q,"","");
    i = i+1;
    next();
  }
  
  summary = entrez_summary(db = "snp", web_history = ids_for_snp$web_history)

  if (num_ids == 1){
    fxn_class = summary$fxn_class;
    clinical_significance = summary$clinical_significance;
  } else{
  fxn_class = summary[[1]]$fxn_class;
  clinical_significance = summary[[1]]$clinical_significance;
  }
  snp_info[i,] = list(q, fxn_class, clinical_significance);
  i = i+1;
}
  

