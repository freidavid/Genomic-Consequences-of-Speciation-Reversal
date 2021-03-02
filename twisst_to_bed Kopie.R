#Summarize twisst results in bed file
#Functions takes a list of twisst output files and makes a bed file 
#with all intorgressed blocks highlighted as 1. 

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
############################################################################
####################################################################################
####################################################################################
twisst_to_bed<-function(weights_list,window_data_file,topology_threshold){
  ########## read data ##################
  #weights = read.table("/Users/david/Desktop/twisst/new/gf.no_downsampling.phyml.50000_DF123470.csv.gz", header = T)
  weights= read.table(weights_list[1], header = T)
  
  
  #normalise rows so weights sum to 1
  weights <- weights / apply(weights, 1, sum)
  #retrieve the names of the topologies
  topoNames = names(weights)
  
  #window_data = read.table("/Users/david/Desktop/twisst/new/phyml.50000_no_downsampling.data.tsv", header = T)
  window_data = read.table(window_data_file, header = T)
  
  
  #exclude any rows where data is missing
  good_rows = which(is.na(apply(weights,1,sum)) == F)
  weights <- weights[good_rows,]
  window_data = window_data[good_rows,]
  #take only the relevant colums from window file
  borders<-window_data[,1:5]
  #Bind togehter with weights file
  df<-cbind(borders,weights)
  #make a bed file
  df<-df[,1:3]
  colnames(df) <-c("chrom","chromStart","chromEnd")
  
  
  ####################################################################################
  ####################################################################################
  #Add individuals
  #########################################################################################
  
  
  for (i in 1:length(weights_list)){
    
    weights= read.table(weights_list[i], header = T)
    
    
    #normalise rows so weights sum to 1
    weights <- weights / apply(weights, 1, sum)
    #retrieve the names of the topologies
    topoNames = names(weights)
    #only keep the good rows
    good_rows = which(is.na(apply(weights,1,sum)) == F)
    weights <- weights[good_rows,]
    #add a column of zeroes
    df<-cbind(df,rep(0,length(df[,1])) )
    #rename that colum
    colnames(df)[i+3] <-paste0("ind_",i)
    #highlight introgressed windows
    df[which(weights$topo3>topology_threshold),3+i]<-1
  }
  
  #Rmove columns that have only zeroes
  #df<-df[ which(rowSums(df[,4:(3+length(weights_list))])>0),]
  
  
  
  #return data frame
  return(df) 
  
}   
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

#define parameters
topology_threshold=2/3
#link to window file from phyml
window_data_file="/Users/david/Dropbox/Final_Figures/TWISST/phyml.50000_final_master.data.tsv"


#list of twisst output files
weights_list<-c("/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF123446.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF123448.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF131.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF121.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF122.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF123.csv",
                
                "/Users/david/Dropbox/Final_Figures/TWISST/gf.final_downsampling.phyml.50000_DF123470.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/gf.final_downsampling.phyml.50000_DF123458.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/gf.final_downsampling.phyml.50000_DF132.csv",
                
                "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF123477.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF123440.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF126.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF127.csv",
                "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF128.csv")
                
                
                
weights_list_bf<-c("/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF123446.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF123448.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF131.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF121.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF122.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/bf.final_downsampling.phyml.50000_DF123.csv")
                   
weights_list_sf<-c( "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF123477.csv",
                    "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF123440.csv",
                    "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF126.csv",
                    "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF127.csv",
                    "/Users/david/Dropbox/Final_Figures/TWISST/sf.final_downsampling.phyml.50000_DF128.csv")

weights_list_gf<-c("/Users/david/Dropbox/Final_Figures/TWISST/gf.final_downsampling.phyml.50000_DF123470.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/gf.final_downsampling.phyml.50000_DF123458.csv",
                   "/Users/david/Dropbox/Final_Figures/TWISST/gf.final_downsampling.phyml.50000_DF132.csv")
                   
total<-twisst_to_bed(weights_list,window_data_file,topology_threshold)
total_bf<-twisst_to_bed(weights_list_bf,window_data_file,topology_threshold)
total_gf<-twisst_to_bed(weights_list_gf,window_data_file,topology_threshold)
total_sf<-twisst_to_bed(weights_list_sf,window_data_file,topology_threshold)


####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
frequency_of_introgression<-function(output_from_twisst_to_bed){
  
  df<-output_from_twisst_to_bed[,4:length(output_from_twisst_to_bed[1,])]
  freq<-rowSums(df[,1:(length(df[1,]))])
  df<-cbind(output_from_twisst_to_bed[,1:3],freq)
  return(df)
  
}
####################################################################################
total_freq<-frequency_of_introgression(total)
bf_freq<-frequency_of_introgression(total_bf)
gf_freq<-frequency_of_introgression(total_gf)
sf_freq<-frequency_of_introgression(total_sf)



####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
new_overlap<-function(a_file,b_file){
  
  library(dplyr)
  joint1<-semi_join(x=a_file[,1:4],y=b_file[,1:3])
  joint2<-semi_join(x=b_file[,1:4],y=a_file[,1:3])
  
  df<-cbind(joint1,joint2$freq,joint1$freq+joint2$freq)
  colnames(df)[4]<-"freq_pop1"
  colnames(df)[5]<-"freq_pop2"
  colnames(df)[6]<-"freq_total"
  return(df)
}
####################################################################################

bf_gf<-new_overlap(bf_freq,gf_freq)
bf_sf<-new_overlap(bf_freq,sf_freq)
gf_sf<-new_overlap(gf_freq,sf_freq)

################################################################################################################################################################################################################################################################################################################################################
####################################################################################
####################################################################################
####################################################################################
#get only introgressed windows
final<-total_freq[total_freq$freq>0,]


#Write Bed Files with only introgressed windows
write.table(as.matrix(final[,1:3]),file="/Users/david/Dropbox/Final_Figures/TWISST/bed_file_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)






############################
#bed file with windows under selection (needs merged final from Final_TWISST_November_2020.R)
chroms<-as.vector(unique(total_freq$chrom))
merged_final2<-merged_final
c<-rep(chroms[1],length(merged_final2[which(merged_final2$chr==1),1]))

for (i in 2:length(chroms)){
  c<-c(c,rep(chroms[i],length(merged_final2[which(merged_final2$chr==i),1])))
         }
length(c)==length(merged_final2$chr)
merged_final2<-cbind(c,merged_final2[,2:(length(merged_final2[1,])-1)])
colnames(merged_final2)<-c("chrom","chromStart","chromEnd","nSL","freq")
merged_final2<-format(merged_final2,scientific=F)

#get only selected sites
selected_bed<-merged_final2[which(merged_final2$nSL>=top),1:3]
intro_bed<-merged_final2[which(merged_final2$freq>0),1:3]



non_intro_bed<-merged_final2[which(merged_final2$freq==0),1:3]
non_selected_non_intro_bed<-merged_final2[which(merged_final2$freq==0&merged_final2$nSL<top),1:3]
selected_and_intro_bed<-merged_final2[which(merged_final2$freq>0&merged_final2$nSL>=top),1:3]
selected_and_non_intro_bed<-merged_final2[which(merged_final2$freq==0&merged_final2$nSL>=top),1:3]
non_selected_and_intro_bed<- merged_final2[which(merged_final2$freq>0&merged_final2$nSL<top),1:3]
non_selected<- merged_final2[which(merged_final2$nSL<top),1:3]









write.table(format(selected_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top1_bed_file_all_selected_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)
write.table(format(intro_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_all_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)
write.table(format(non_intro_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_non_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)
write.table(format(non_selected,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_non_selected_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)

write.table(format(non_selected_non_intro_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_non_selected_non_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)
write.table(format(selected_and_intro_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_selected_and_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)
write.table(format(selected_and_non_intro_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_selected_and_not_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)
write.table(format(non_selected_and_intro_bed,scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/top_1_bed_file_not_selected_and_introgressed_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)

write.table(format(merged_final2[,1:3],scientific=F),file="/Users/david/Dropbox/Final_Figures/nSL/bed_files/all_windows.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names = FALSE)

help(write.table)
which(duplicated(rbind(merged_[,1:3],non_selected_and_intro_bed)[,1:3]))
help(format)

length(intro_bed[,1])
length(selected_and_intro_bed[,1])
length(non_selected_and_intro_bed[,1])


length(selected_bed[,1])
length(selected_and_intro_bed[,1])+length(selected_and_non_intro_bed[,1])





