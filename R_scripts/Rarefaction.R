blocks<-c(1004,1678,2937,3719,4608,5089)
plot(blocks)
inds<-seq(1,6,1)

xlim=c(0,50)
ylim=c(0,6000)
fit = lm(blocks ~ log(inds))
coefficients<-coef(fit)
plot(inds, blocks,type='p',col='blue',bty="n",xlim=xlim,ylim=ylim)
saturation_curve<-curve(coefficients[1]+coefficients[2] *log(x),add=TRUE,col="red")
text(max(xlim),coefficients[1]+coefficients[2]*log(max(xlim)),round(coefficients[1]+coefficients[2]*log(max(xlim)),digits=1),xpd=T,pos=3)

library(vegan)
rarefy()
help(rarefy)
test
bf_25$weight_1
rarefy(50000,c(400,400,500,500,600,600),se=T,MARGIN=1)
data(dune)



window_data = read.table("twisst_window_data.tsv", header = T)

#exclude any rows where data is missing
good_rows = which(is.na(apply(weights,1,sum)) == F)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]
#take only the relevant colums from window file
borders<-window_data[,1:5]
length(borders[,1])
i=2



########## read data ##################
weights_1 = read.table(weights_file_1, header = T)
weights_2 = read.table(weights_file_2, header = T)

#normalise rows so weights sum to 1
weights_1 <- weights_1 / apply(weights_1, 1, sum)
weights_2 <- weights_2 / apply(weights_2, 1, sum)


window_data = read.table(window_data_file, header = T)







wind_id<-vector()
#Create a vector with chr_position as "species list"
for(i in 1:length(borders$scaffold)){
  wind_id[i]<-paste0(borders[i,1],"_",borders[i,2])
}

#create a data frame
zero<-rep(0,length(borders$scaffold))
frame<-rbind(zero,zero,zero,zero,zero,zero)


intr_list<-list()
intr_list[[1]]<-bf_25$weights_1
intr_list[[2]]<-bf_25$weights_2
intr_list[[3]]<-sf_25$weights_1
intr_list[[4]]<-sf_25$weights_2
intr_list[[5]]<-gf_25$weights_1
intr_list[[6]]<-gf_25$weights_2

wind_id_intr_list<-list()

for(k in 1:length(intr_list)){

data_file<-intr_list[[k]]

ind_id<-vector()
#Create a vector with chr_position as "species list"
for(p in 1:length(data_file[,1])){
  ind_id[p]<-paste0(data_file[p,1],"_",data_file[p,2])
}
wind_id_intr_list[[k]]<-ind_id






}

#Loop through all individuals and mark introgressed windows with 1 in "frame"
for(j in 1:length(wind_id_intr_list)){
  vec<-vector()
  relevant_winds<-wind_id_intr_list[[j]]
  
  for(i in 1:length(relevant_winds)){
    
    for(g in 1:length(wind_id))
      if(wind_id[g]==relevant_winds[i]){
        frame[j,g]<-1}
  }
}



vec
length(vec)
legnth()
get_back[2]
rarefy(frame,10000,MARGIN=1)

inds<-seq(1,100,1)
inds
test<-((0.93536063^inds)*1255.333333)
test2<-vector()
test2[1]<-1255.333333
for(i in 2:length(test)){
  test2[i]<-test2[i-1]+test[i]
}

#Transform to percentage of windows
test3<-test2/76409

plot(test3,type="l",col="red",xlab="Individuals",ylab="Proportion of the genome")
text(-8,max(test3)*0.95,round(max(test3),digits=2),xpd=T,pos=3,font=2)
abline(h=max(test3),xpd=F,lty=3)
#points(blocks,pch=19)

18224/76409

