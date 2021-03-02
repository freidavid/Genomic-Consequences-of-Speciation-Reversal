
###############################################################################################
#Preparation of Twisst output:
###############################################################################################

#Define colors
library(RColorBrewer)
colors<-brewer.pal(12, name="Paired")
greys<-brewer.pal(9, name="Greys")
GF = colors[4]
SF = colors[8]
BF = colors[2]
K = greys[3]

#function to generate transparent colors
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#First: import the data, and plot manhattan plots for each individuals, with the introgressed windows highlighted in color
#Function twisst_manhattan:
twisst_man<-function(weights_file=weights_file,topology_threshold=0.66,only_points=F,window_data_file=window_data_file,min_sites_per_window=0,max_values=0,wind_correction=0,limits=c(0,1),smooth=F,sep_lines=T,text_size=1.3,span_factor=1,Cols="black",line_width=1,chromo_range=NULL,bp_output=T,data_frame_known=F,data_frame=NULL){
  
  
  #function to generate transparent colors
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  
  
  
  if(data_frame_known==F){
    ########## read data ##################
    weights = read.table(weights_file, header = T)
    #normalise rows so weights sum to 1
    weights <- weights / apply(weights, 1, sum)
    #retrieve the names of the topologies
    topoNames = names(weights)
    
    window_data = read.table(window_data_file, header = T)
    
    #exclude any rows where data is missing
    good_rows = which(is.na(apply(weights,1,sum)) == F)
    weights <- weights[good_rows,]
    window_data = window_data[good_rows,]
    #take only the relevant colums from window file
    borders<-window_data[,1:5]
    #Bind togehter with weights file
    df<-cbind(borders,weights)
    
    
    #rename dataframe because i've written this function for fst actually
    fst_file<-df
    #get the chromosome names
    temp  <- fst_file[ (fst_file$sites >= min_sites_per_window),]
    
    
    chr<-unique(temp[,1])
    chromos<-chr[1:31]
    chr<-chromos
    
    
    temp2<-NULL
    temp2<-subset(temp,temp[,1]==chromos[1])
    for( i in 2:length(chromos)){
      temp2<-rbind(temp2,subset(temp,temp[,1]==chromos[i]))
    }
    
    temp<-temp2
    
    b<-vector()
    for(i in 1:length(chromos)){
      a<-rep(i,length(temp[temp[,1]==chromos[i],1]))
      b<-c(b,a)
    }
    temp[,1]<-b 

    # construct dataframe with with relevant information (Chr,Pos,Cov)
    chr <- as.numeric(temp$scaffold)
    chr[is.na(chr)]<-0
    pos<-as.numeric(((temp$start-1)+temp$end)/2)
    #pos <- as.numeric(temp$mid)
    pos<-pos+wind_correction
    pos[is.na(pos)]<-0
    fst <- as.numeric(temp$topo3)
    coverage<-as.numeric(temp$sites)
    wind_fst  <- data.frame(chr, pos, fst,coverage)
    pch<-rep(19,length(pos))
    cov_fst  <- data.frame(chr, pos, fst,coverage,pch)
    wind_fst
    
    #take only the ones with defined chromosome
    df11<-subset(wind_fst,wind_fst[,1]!=0)
    colnames(df11)[1] <- "CHR"
    colnames(df11)[2] <- "BP"
    colnames(df11)[3] <- "P"
    colnames(df11)[4] <- "Coverage"
    
    
    
    if(bp_output==T){
      df11_result<-df11
    }
    #############################################################
  }
  
  if(length(chromo_range)==1){
    chromo_range<-rep(chromo_range[1],2)
  }
  
  if(data_frame_known==T){
    df11<-data_frame
  }
  
  
  if(!is.null(chromo_range)){
    start_chromo<-chromo_range[1]
    end_chromo<-chromo_range[2]
    df11<-subset(df11,df11$CHR<=end_chromo)
    df11<-subset(df11,df11$CHR>=start_chromo)
    chromos<-seq(start_chromo,end_chromo,1)
  }
  
  #get longest of each chromosomes for all individuals
  chr<-unique(df11$CHR)
  if(max_values[1]==0){
    #loop to get length of each chromosome in each file in list
    y<-vector()
    for(j in 1:length(chromos)){
      y[j]<- max(as.numeric(temp[temp[,1]==chr[j],3]))
    }
    
    
    #get longest of each chromosomes
    max_values<-vector()
    for (g in 1:length(chromos)){
      max_values[g]<-max(y[g])
    }
  }else(max_values<-max_values)
  
  max_values
  if(!is.null(chromo_range)){
    max_values<-max_values[chromo_range[1]:chromo_range[2]]
  }
  ################################################################
  ################################################################
  #Here starts my actual own plot
  ################################################################
  ################################################################
  
  max_values
  #Set up an empty manhattan plot
  total<-sum(as.numeric(max_values))
  
  #get percentage on x-axis for each chromosome
  percentages<-max_values/total
  
  
  #sum it up to get 100 percent in total
  endlines<-percentages[1]
  for(c in 2:length(percentages)){
    endlines[c]<-c((percentages[c]+endlines[c-1]))
  }
  #Change to percent
  
  endlines<-endlines*100
  
  #move the chromosome label to the center of the chromosome
  axis<-endlines
  for(i in 1:length(chromos)){
    axis[i]<-axis[i]-((percentages[i]*100)/2)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  #make an empty plot
  if(only_points==F){
    par(mar=c(6,6,2,17),mgp=c(3.5,1,0),family="Times",xpd=T)
    
    plot(df11[,2],df11[,3],ylim=limits,xlim=c(0,max(endlines)),col=df11$Colour,pch=16,axes=F,ylab=expression("Topology weighting"),xlab="Chromosome",type="n",yaxs="i",cex.lab=text_size,cex.axis=text_size)
    
    
    #Create the background shading according to chromosome (alternating)
    rect(0,limits[1],endlines[length(chr)],limits[2],col=add.alpha("grey92",0.99),border=NA)
    for(i in 1:length(chromos)){
      rect(endlines[(i*2)-1],limits[1],endlines[i*2],limits[2],col="white",border=NA )
    }
  }
  

  ############################################################################################
  if(data_frame_known==F){
    #For the first data file: Calculate the percentage of each chromosome on the plot like above
    start<-c(0,endlines)
    for (h in 1:length(chr)){
      df11[df11$CHR==h,2]<-(df11[df11$CHR==h,2]/total*100)+start[h]
    }
  }
  
  if(only_points==T){
    
  }else{
    if(!is.null(chromo_range)){
      axis(side=1,at=axis,labels=seq(chromo_range[1],chromo_range[2],1),cex.axis=text_size)
      axis(side=2,cex.axis=text_size)}else{
        axis(side=1,at=axis,labels=seq(1,length(chromos),1),cex.axis=text_size)
        axis(side=2,cex.axis=text_size)
      }
  }
  par(xpd=F)
  #get other symbol for topology weight over threshold
  df11_h<-subset(df11,df11$P>=topology_threshold)
  df11_l<-subset(df11,df11$P<=topology_threshold)
  
  par(mar=c(6,6,2,17),mgp=c(3.5,1,0),family="Times",xpd=T)
  #Now make the actual plot: if sep lines is true, plot every chromosome separately, if false draw only one line over whole genome
  if(sep_lines==TRUE){
    for(t in 1:length(chr)){
      points(x=df11_l[df11_l$CHR==t,]$BP,y=df11_l[df11_l$CHR==t,]$P,col=add.alpha("grey",0.5),pch=20)
      
    }
  }else{points(x=df11_l$BP,y=df11_l$P,col=Cols[1],lwd=line_width,pch=20)}
  #legend(x="right",cex=l_text_size,inset=inset,xpd=T,bty="n",title=expression(bold("Pairwise Comparison")),legend=comp,col=Cols,lwd=line_width,pch=16)
  #Now make the actual plot: if sep lines is true, plot every chromosome separately, if false draw only one line over whole genome
  
  if(sep_lines==TRUE){
    for(t in 1:length(chr)){
      points(x=df11_h[df11_h$CHR==t,]$BP,y=df11_h[df11_h$CHR==t,]$P,col=Cols[1],lwd=line_width,pch=19)
      
    }
  }else{points(x=df11_h$BP,y=df11_h$P,col=Cols[1],lwd=line_width,pch=19)}
  #legend(x="right",cex=l_text_size,inset=inset,xpd=T,bty="n",title=expression(bold("Pairwise Comparison")),legend=comp,col=Cols,lwd=line_width,pch=16)
  
  
  if(bp_output==F){
    df11_result<-df11
  }
  if(data_frame_known==F){
    list<-list()
    list[[1]]<-df11_result
    list[[2]]<-max_values
    list[[3]]<-endlines
    return(list)
  }
}


#Manhattan plots for all individuals
int_121<-twisst_man(weights_file = "./bf.final_downsampling.phyml.50000_DF121.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=BF,min_sites_per_window = 25)
int_122<-twisst_man(weights_file = "./bf.final_downsampling.phyml.50000_DF122.csv",window_data_file = "./phyml.50000_final_master.data.tsv",Cols=BF ,min_sites_per_window = 25)
int_123<-twisst_man(weights_file = "./bf.final_downsampling.phyml.50000_DF123.csv",window_data_file = "./phyml.50000_final_master.data.tsv",Cols=BF,min_sites_per_window = 25 )
int_131<-twisst_man(weights_file = "./bf.final_downsampling.phyml.50000_DF131.csv",window_data_file = "./phyml.50000_final_master.data.tsv",Cols=BF ,min_sites_per_window = 25)
int_123446<-twisst_man(weights_file = "./bf.final_downsampling.phyml.50000_DF123446.csv",window_data_file = "./phyml.50000_final_master.data.tsv",Cols=BF,min_sites_per_window = 25 )
int_123448<-twisst_man(weights_file = "./bf.final_downsampling.phyml.50000_DF123448.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=BF,min_sites_per_window = 25)
 


int_123458<-twisst_man(weights_file = "./gf.final_downsampling.phyml.50000_DF123458.csv",window_data_file = "./phyml.50000_final_master.data.tsv",Cols=GF ,min_sites_per_window = 25)
int_123470<-twisst_man(weights_file = "./gf.final_downsampling.phyml.50000_DF123470.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=GF,min_sites_per_window = 25)
int_132<-twisst_man(weights_file = "./gf.final_downsampling.phyml.50000_DF132.csv",window_data_file = "./phyml.50000_final_master.data.tsv",Cols=GF,min_sites_per_window = 25 )


int_123477<-twisst_man(weights_file = "./sf.final_downsampling.phyml.50000_DF123477.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=SF,min_sites_per_window = 25)
int_123440<-twisst_man(weights_file = "./sf.final_downsampling.phyml.50000_DF123440.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=SF,min_sites_per_window = 25)
int_126<-twisst_man(weights_file = "./sf.final_downsampling.phyml.50000_DF126.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=SF,min_sites_per_window = 25)
int_127<-twisst_man(weights_file = "./sf.final_downsampling.phyml.50000_DF127.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=SF,min_sites_per_window = 25)
int_128<-twisst_man(weights_file = "./sf.final_downsampling.phyml.50000_DF128.csv",window_data_file = "./phyml.50000_final_master.data.tsv" ,Cols=SF,min_sites_per_window = 25)





#Now create the circos plot
library(circlize)


#define parameteres:
topology_threshold=0.66
snps_per_window=25
wind_size=50000

#get window data
jh<-int_121[[1]]

#Replace chromosome names with numbers
names<-c(1:31)
replace<-rep(names[1],length(jh[which(jh[,1]==1),1]))
for(i in 2:length(names)){
    replace<-c(replace,rep(names[i],length(jh[which(jh[,1]==i),1])))
}
jh$CHR<-replace


#Get data from all individuals
DF121<-subset(int_121[[1]],int_121[[1]]$P>topology_threshold)
DF122<-subset(int_122[[1]],int_122[[1]]$P>topology_threshold)
DF123<-subset(int_123[[1]],int_123[[1]]$P>topology_threshold)
DF123446<-subset(int_123446[[1]],int_123446[[1]]$P>topology_threshold)
DF123448<-subset(int_123448[[1]],int_123448[[1]]$P>topology_threshold)
DF131<-subset(int_131[[1]],int_131[[1]]$P>topology_threshold)

DF132<-subset(int_132[[1]],int_132[[1]]$P>topology_threshold)
DF123458<-subset(int_123458[[1]],int_123458[[1]]$P>topology_threshold)
DF123470<-subset(int_123470[[1]],int_123470[[1]]$P>topology_threshold)


DF123440<-subset(int_123440[[1]],int_123440[[1]]$P>topology_threshold)
DF123477<-subset(int_123477[[1]],int_123477[[1]]$P>topology_threshold)
DF126<-subset(int_126[[1]],int_126[[1]]$P>topology_threshold)
DF127<-subset(int_127[[1]],int_127[[1]]$P>topology_threshold)
DF128<-subset(int_128[[1]],int_128[[1]]$P>topology_threshold)





#Replace chromo names with numbers

for(i in 1:length(names)){
  DF121[which(DF121[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF122[which(DF122[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123[which(DF123[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123446[which(DF123446[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123448[which(DF123448[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF131[which(DF131[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123458[which(DF123458[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123470[which(DF123470[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF132[which(DF132[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123440[which(DF123440[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF123477[which(DF123477[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF126[which(DF126[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF127[which(DF127[,1]==i),1]<-names[i]}

for(i in 1:length(names)){
  DF128[which(DF128[,1]==i),1]<-names[i]}

###############################################################################################
###############################################################################################
#Preparation of selscan output (nSL):
###############################################################################################
library(readr)
#read in and store in a list
chromos<-c(1:3,5:6,8:16,18:21,23:27,29:31,33:34,36,39:40)
selscan<-list()
for(i in c(1:31)){
  selscan[[i]]<- read_delim(paste0("./scaffold",(chromos[i]-1),".selscan.nsl.out.20bins.norm.50kb.windows"), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)}


#Add a chromosome number to the list
selscan2<-list()
chromos<-1:31
for(i in c(1:31)){
  selscan2[[i]]<-cbind(rep(as.character(names[i]),length(as.vector(selscan[[chromos[i]]]$X1))),(selscan[[chromos[i]]]))}



#make one dataframe with all chromosomes
sel<-as.data.frame(selscan2[[1]])
for(i in c(2:31)){
  sel<-rbind(sel,selscan2[[i]])}
sel<-as.data.frame(sel)

#Replace chromosome names with numbers:
unique(sel[,1])

#add column with how many snps per window are above threshold
sel<-cbind(sel,round(sel[,4]*sel[,5],digits=3))
sel[which(sel[,5]==max(na.omit(sel[,5]))),]
#plot against index
plot(sel[which(sel[,1]==40),5])
length(sel[which(sel[,1]==1),5])
names<-colnames(sel)




#change names of df
colnames(sel)<-c("chr","start","end","snps","nSL",names[6:7],"above_cutoff")

#correct different window counting
sel[,3]<-sel[,3]-1

#Replace -1 with 0
sel$nSL[sel$nSL==-1]<-0


###############################################################################################
#Merge Twisst and Selscan output:
#Important: Here, the output from twisst_to_bed.R is used instead of the output from above. This
#is because it needs a file with introgressed windows across all individuals, which I did in the
#mentioned script already. File is called total_freq

#Replace chromosome names with numbers in total_freq
chr<-unique(total_freq[,1])
replace_with<-c(1:31)
a<-rep(replace_with[1],length(total_freq[total_freq[,1]==chr[1],1]))
for(i in 2:length(chr)){
  a<-c(a,rep(replace_with[i],length(total_freq[total_freq[,1]==chr[i],1])))}
total_freq2<-cbind(a,total_freq[,2:length(total_freq[1,])])


#Merge selection scan and total_freq2
#chromosome 1
sel_chr1<-sel[which(sel[,1]==1),]
total_freq_chr1<-total_freq2[total_freq2[,1]==1,]
sel_freq_merged<-merge(sel_chr1,total_freq_chr1,by.x="start",by.y="chromStart")
chromos<-c(1:31)
#all the rest
for(i in 2:length(chr)){
  sel_int<-sel[which(sel[,1]==chromos[i]),]
  total_freq_int<-total_freq2[total_freq2[,1]==chromos[i],]
  merged_int<-merge(sel_int,total_freq_int,by.x="start",by.y="chromStart")
  sel_freq_merged<-rbind(sel_freq_merged,merged_int)
}

#extract and resort relevant columns
merged_final<-data.frame(sel_freq_merged$chr,sel_freq_merged$start,sel_freq_merged$end,sel_freq_merged$nSL,sel_freq_merged$freq)
colnames(merged_final)<-c("chr","start","end","nSL","freq")


farben<-rep(NA,length(merged_final[,1]))
for (i in 1:length(farben)){
  #if((merged_final[i,5]>0)&(merged_final[i,4]>=0.2))
  if((merged_final[i,5]>0))
    {farben[i]<-"darkred"}
  else{farben[i]<-"black"}
}

merged_final<-cbind(merged_final,farben)


###############################################################################################


dev.off()
pdf("./Figure3_circos.pdf")
#Create circos tracks
library(circlize)
circos.clear()
circos.par("track.height" = 0.08,start.degree=90-(25/2),gap.degree=c(rep(1,30),25),cell.padding=c(0.005,0.005,0.005,0.005))#
circos.initialize(factors = jh$CHR, x = jh$BP)
#Create the first track with chromosome labels
BF<-add_transparency(BF,0.5)
GF<-add_transparency(GF,0.5)
SF<-add_transparency(SF,0.5)

#selection scan track
circos.track(factors = jh$CHR, ylim = c(0,1),bg.col=K,bg.border=NA,track.height=(0.075),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index ,niceFacing = T,cex=0.3)
             })
circos.yaxis("left",at=c(0,1),sector.index = 4,labels.cex=0.5,labels = c(0,1))

#Wartmanni Track
circos.track(factors = jh$CHR,ylim = c(0,1),bg.col=BF,bg.border=NA,track.height=0.15)
circos.yaxis("left",at=c(0+(1/12),1/6+(1/12),2/6+(1/12),3/6+(1/12),4/6+(1/12),5/6+(1/12)),sector.index = 1,labels.cex=0.5,labels = c(1:6,""))

#Macrophtlamus track
circos.track(factors = jh$CHR, y = jh$P,bg.col=GF,bg.border=NA,track.height=0.075)
circos.yaxis("left",at=c(0+(1/6),1/3+(1/6),2/3+(1/6)),sector.index = 1,labels.cex=0.5,labels = c(1:3,""))

#Arenicolus track
circos.track(factors = jh$CHR, y = jh$P,bg.col=SF,bg.border=NA,track.height=0.125)
circos.yaxis("left",at=c(0+(1/10),1/5+(1/10),2/5+(1/10),3/5+(1/10),4/5+(1/10)),sector.index = 1,labels.cex=0.5,labels = c(1:5,""))

#Plot selection scan into outermost tradck
for(i in 1:length(merged_final[,1])){
  circos.points(x=((merged_final[i,2]-1)+merged_final[i,3])/2,y=merged_final[i,4],track.index=1,col=merged_final[i,6],sector.index = merged_final[i,1],pch=16,cex=0.1)
}

#Now plot introgressed windows into tracks
#wartmanni
for(i in 1:length(DF121[,1])){
  circos.rect(DF121[i,2]-(wind_size/2),0, DF121[i,2]+(wind_size/2),1/6 ,track.index = 2,col="black",sector.index=DF121[i,1],border=NA)
}


for(i in 1:length(DF122[,1])){
  circos.rect(DF122[i,2]-(wind_size/2),1/6, DF122[i,2]+(wind_size/2),2/6 ,track.index = 2,col="black",sector.index=DF122[i,1],border=NA)
}

for(i in 1:length(DF123[,1])){
  circos.rect(DF123[i,2]-(wind_size/2),2/6, DF123[i,2]+(wind_size/2),3/6 ,track.index = 2,col="black",sector.index=DF123[i,1],border=NA)
}


for(i in 1:length(DF123446[,1])){
  circos.rect(DF123446[i,2]-(wind_size/2),3/6, DF123446[i,2]+(wind_size/2),4/6 ,track.index = 2,col="black",sector.index=DF123446[i,1],border=NA)
}


for(i in 1:length(DF123448[,1])){
  circos.rect(DF123448[i,2]-(wind_size/2),4/6, DF123448[i,2]+(wind_size/2),5/6 ,track.index = 2,col="black",sector.index=DF123448[i,1],border=NA)
}

for(i in 1:length(DF131[,1])){
  circos.rect(DF131[i,2]-(wind_size/2),5/6, DF131[i,2]+(wind_size/2),6/6 ,track.index = 2,col="black",sector.index=DF131[i,1],border=NA)
}


#macrophtalmus
for(i in 1:length(DF132[,1])){
  circos.rect(DF132[i,2]-(wind_size/2),0, DF132[i,2]+(wind_size/2),1/3 ,track.index = 3,col="black",sector.index=DF132[i,1],border=NA)
}

for(i in 1:length(DF123458[,1])){
  circos.rect(DF123458[i,2]-(wind_size/2),1/3, DF123458[i,2]+(wind_size/2),2/3 ,track.index = 3,col="black",sector.index=DF123458[i,1],border=NA)
}

for(i in 1:length(DF123470[,1])){
  circos.rect(DF123470[i,2]-(wind_size/2),2/3, DF123470[i,2]+(wind_size/2),3/3 ,track.index = 3,col="black",sector.index=DF123470[i,1],border=NA)
}

#arenicolus
for(i in 1:length(DF126[,1])){
  circos.rect(DF126[i,2]-(wind_size/2),0, DF126[i,2]+(wind_size/2),1/5 ,track.index = 4,col="black",sector.index=DF126[i,1],border=NA)
}

for(i in 1:length(DF127[,1])){
  circos.rect(DF127[i,2]-(wind_size/2),1/5, DF127[i,2]+(wind_size/2),2/5 ,track.index = 4,col="black",sector.index=DF127[i,1],border=NA)
}

for(i in 1:length(DF128[,1])){
  circos.rect(DF128[i,2]-(wind_size/2),2/5, DF128[i,2]+(wind_size/2),3/5 ,track.index = 4,col="black",sector.index=DF128[i,1],border=NA)
}


for(i in 1:length(DF123440[,1])){
  circos.rect(DF123440[i,2]-(wind_size/2),3/5, DF123440[i,2]+(wind_size/2),4/5 ,track.index = 4,col="black",sector.index=DF123440[i,1],border=NA)
}

for(i in 1:length(DF123477[,1])){
  circos.rect(DF123477[i,2]-(wind_size/2),4/5, DF123477[i,2]+(wind_size/2),5/5 ,track.index = 4,col="black",sector.index=DF123477[i,1],border=NA)
}




dev.off()



###############################################################################################



#Looking at the sharing of introgressed windows for heatmap in the middle
overlapping_windows<-function(weights_file_1,weights_file_2,window_data_file,topology_threshold,permutation=FALSE,iterations=100){
  
  
  ########## read data ##################
  weights_1 = read.table(weights_file_1, header = T)
  weights_2 = read.table(weights_file_2, header = T)
  
  #normalise rows so weights sum to 1
  weights_1 <- weights_1 / apply(weights_1, 1, sum)
  weights_2 <- weights_2 / apply(weights_2, 1, sum)
  
  
  window_data = read.table(window_data_file, header = T)
  
  
  
  #get the chromosomes
  chromos<-unique(window_data$scaffold)
  chromos<-chromos[1:31]
  chromosome_data<-list()
  #split window files for every chromosome and store it to list
  for(i in 1:31){
    chromosome_data[[i]]<-subset(window_data,window_data$scaffold==chromos[i])
  }
  
  #same for weights
  chromosome_weights_1<-list()
  chromosome_weights_1[[1]]<-weights_1[1:(length(chromosome_data[[1]]$scaffold)),]
  x<-length(chromosome_data[[1]]$scaffold)
  for (i in 2:31){
    chromosome_weights_1[[i]]<-weights_1[(x+1):(x+length(chromosome_data[[i]]$scaffold)),]
    x<-x+length(chromosome_data[[i]]$scaffold)
  }
  
  
  chromosome_weights_2<-list()
  chromosome_weights_2[[1]]<-weights_2[1:(length(chromosome_data[[1]]$scaffold)),]
  x<-length(chromosome_data[[1]]$scaffold)
  for (i in 2:31){
    chromosome_weights_2[[i]]<-weights_2[(x+1):(x+length(chromosome_data[[i]]$scaffold)),]
    x<-x+length(chromosome_data[[i]]$scaffold)
  }
  
  
  
  
  
  #get the positions of wheight of topo 3 above 0.8 for individual 1
  y<-0
  x<-0
  z<-0
  w<-0
  for(k in 1:31){
    x<-(chromosome_data[[k]][which(chromosome_weights_1[[k]]$topo3>=topology_threshold),1:4])
    if(length(x[,1])>0){
      y<-rbind(y,x)
    }
  }
  #get rid of first rows for individual 1
  
  
  
  y <- y[2:length(y[,1]),1:4]
  #get the positions of wheight of topo above 0.8 for individual 2
  topology_threshold
  for(k in 1:31){
    w<-(chromosome_data[[k]][which(chromosome_weights_2[[k]]$topo3>=topology_threshold),1:4])
    if(length(w[,1])>0){
      z<-rbind(z,w)
    }
  }
  
  #get rid of first rows for individual 2
  z <- z[2:length(z[,1]),1:4]
  
  
  cutoff="NA"
  #permutation  testing:
  if(permutation==T){
    
    observed<-length(y[which(duplicated(rbind(y[2:length(y[,1]),],z[2:length(z[,1]),]),fromLast = T)),1])
    m<-vector()
    for(i in 1:iterations){
      s1<-sample(c(1:length(na.omit(weights_1[,1]))),length(y[,1]),replace=F)
      s2<-sample(c(1:length(na.omit(weights_2[,1]))),length(z[,1]),replace=F)
      t<-duplicated(c(s1,s2))
      m[i]<-length(which(t==TRUE))
      # m[i]<-length((which(s1==s2)))
      
    }
    ######
    
    #####
    
    p_value<-length(which(m>=observed))/iterations
    m<-sort(m)
    cutoff<-m[iterations*0.05]
    
    if(length(which(m>=observed))==0){
      p_value<-paste0("p < ",1/iterations)}
    
    
    
    
  }else{p_value<-"No permutation testing applied"}
  
  
  
  #z is windows above 0.8 for the second file, y is windwos over 0.8 for the first file
  results<-list()
  results[[1]]<-"No overlapping windows"
  results[[2]]<-y
  results[[3]]<-z
  results[[4]]<-c(length(na.omit(weights_1[,1])),length(na.omit(weights_2[,1])))
  results[[5]]<-c(length(y[,1]),length(z[,1]))
  results[[6]]<-length(y[which(duplicated(rbind(y[2:length(y[,1]),],z[2:length(z[,1]),]),fromLast = T)),1])
  results[[7]]<-p_value
  results[[8]]<-cutoff
  #results[[9]]<-max(m)
  
  #get overlapping windows:
  overlapping<-y[which(duplicated(rbind(y[2:length(y[,1]),],z[2:length(z[,1]),]),fromLast = T)),]
  #store those in list[[1]]
  results[[1]]<-overlapping
  
  #name list elements
  names(results)<-c("overlapping","weights_1","weights_2","Number of windows with data (sample 1, sample 2)","Number of windows above topology threshold (sample 1, sample 2)","Number of overlapping windows","p-value","number of overlapping windows for p<0.05")
  if(length(overlapping)==0){
    print("No overlapping windows detected")
    return(results)
  }else{return(results)}
  
  
  
  
  
}
over_overlapping<-function(z,t){
  overlapping<-z[which(duplicated(rbind(z[2:length(z[,1]),],t[2:length(t[,1]),]),fromLast = T)),]
  return(na.omit(overlapping))
}


#with overlapping_windows function, assess all pairwise combinations of individuals and
#create table with proportion of shared windows. Eg:
overlapping_windows( "./sf.final_downsampling.phyml.50000_DF128.csv","./gf.final_downsampling.phyml.50000_DF132.csv", "./phyml.50000_august.data.tsv" ,topology_threshold = 2/3,permutation = T)
overlapping_windows( "./bf.final_downsampling.phyml.50000_DF121.csv","./gf.final_downsampling.phyml.50000_DF132.csv", "./phyml.50000_august.data.tsv" ,topology_threshold = 2/3,permutation = T)
overlapping_windows( "./bf.final_downsampling.phyml.50000_DF121.csv","./gf.final_downsampling.phyml.50000_DF123470.csv", "./phyml.50000_august.data.tsv" ,topology_threshold = 2/3,permutation = T)
#and so on....

#Load manually created table of pairwise overlap form above
library(readr)
#plot heatmap
overlap <- read_delim("./overlap.csv",  ";", escape_double = FALSE, col_names = T,   trim_ws = TRUE)
h<-as.data.frame(overlap)
h[is.na(h)] <- 0
h<-as.numeric(h)
h<-as.matrix(h)
count<-diag(h)
diag(h)<-0
help(matrix)
d<-matrix(NA,nrow=14,ncol=14)
diag(d)<-count



# use Rcolor brewer to define a color scale
rbPal<-colorRampPalette(c("NA",brewer.pal(9,"YlOrRd")))
nclasses<-rbPal(100)
z<-brewer.pal(9,"YlGnBu")
rbPal2<-colorRampPalette(z[4:9])
nclasses2<-rbPal2(100)
nclasses<-add.alpha(nclasses,0.6)
nclasses2<-add.alpha(nclasses2,0.6)

# define the breakpoints
breaksplot <- seq(min(h, na.rm=T)*0.99, max(h, na.rm=T)*1.01,length.out = length(nclasses)+1)
breaksplot2 <- seq(min(d, na.rm=T)*0.99, max(d, na.rm=T)*1.01,length.out = length(nclasses2)+1)


#Function for to plot heatmap:
image.scale <- function(z, zlim, col = rainbow(12), breaks, horiz=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){ylim<-c(0,1); xlim<-range(breaks)}
  if(!horiz){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", xlab="",...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

# plot the image
dev.off()
pdf("./Figure_3_heatmap.pdf",height=9.53,width=11)


#adjust layout
mat<-(rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3)))
layout(mat)

#plot
image(1:nrow(h), 1:nrow(h), t((h)), breaks=breaksplot, xlab="", ylab="", col=nclasses, main="",
      xaxt="n", yaxt="n")

image(1:nrow(d), 1:nrow(d), t(t(d)), breaks=breaksplot2, xlab="", ylab="", col=nclasses2, main="",
      xaxt="n", yaxt="n",add=T)

image.scale(t(h), zlim=c(min(h),max(h)), col = nclasses,breaks=breaksplot, horiz=FALSE, ylab="", main=expression(bold("A)"))
            , cex.lab=1.5, cex.axis=1.5, cex.main=1.5)

image.scale(t(h), zlim=c(min(d),max(d)), col = nclasses2,breaks=breaksplot2, horiz=FALSE, ylab="", main=expression(bold("B)"))
            , cex.lab=1.5, cex.axis=1.5, cex.main=1.5)



dev.off()






#############################################################################################################################################################################################################

#Test if sharing within and between species are different

#get all comparisons within species
within<-h[c(2,3,4,5,6,17,18,19,20,32,33,34,47,48,62,92,93,107,137,138,139,140,152,153,154,167,168,182)]

#get all comparisons between species
between<-h[c(7,8,9,10,11,12,13,14,21,22,23,24,25,26,27,28,35,36,37,38,39,40,41,42,49,50,51,52,53,54,55,56,63,64,65,66,67,68,69,70,77,78,79,80,81,82,83,84,94,95,96,97,98,108,109,110,111,112,122,123,124,125,126)]

#t-test
t.test(within,between)

#############################################################################################################################################################################################################




#############################################################################################################################################################################################################
#permutation testing across all for 
#enrichment of selected windows in all contmeporary species
#############################################################################################################################################################################################################
permutations=10000
expected <- vector()
observed<-length(merged_final[which(merged_final$nSL>0.51&merged_final$freq>=1),1])

for (i in 1:permutations){
#first: sample random number of introgressed windows of whole genome:
intro_numb<-length(merged_final[which(merged_final$freq>0),5])
index<-sort(sample(length(merged_final[,1]),intro_numb))
sampled<-merged_final[index,]

#calculate number of selected windows for this sample
sel_sampled<-sampled[which(sampled$nSL>=0.51),]
#calulate the number of windows under selection in sampled introgressed windows
expected[i]<-length(sel_sampled[,1])}

#calculated p_value
length(which(observed<=expected))/permutations


#############################################################################################################################################################################################################








