##Matt Kaufman and Vanessa Garayburu-Caruso vanessa.garayburu-caruso@pnnl.gov Pacific Northwest National Laboratory
##This script loads the list of SSS sites from the published metadata and loops through them all,
##running the stream metabolizer template for each one
rm(list=ls(all=T))
metadata.path="C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/SSS_metabolism/Published_Data/v3_SSS_Data_Package/" #from this data package: https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566

output.path= "Outputs/"

metadata=read.csv(paste(metadata.path,'v2_SSS_Metadata_IGSN-Mapping.csv',sep=''),skip=1,header=T)
metadata <- metadata[grepl("Water", metadata$Sample_Name, ignore.case = TRUE),]

SITE_LIST<-data.frame(substring(metadata$Sample_Name,0,6),metadata$Locality)
colnames(SITE_LIST)<-c('PARENT_ID','SITE_ID')
for(i in 1:length(SITE_LIST[,1])) { 
PARENT_ID<-SITE_LIST[i,1]
SITE_ID<-SITE_LIST[i,2]
htmlfilename=paste0(getwd(),"/",output.path,'v2_',PARENT_ID,'_',SITE_ID,'_SM_output.HTML')
  
rmarkdown::render("Scripts/v2_SSS_SM_final_template.Rmd",
                  output_file = htmlfilename,
                  params = list(
                    PARENT_ID=PARENT_ID,
                    SITE_ID=SITE_ID
))

}
