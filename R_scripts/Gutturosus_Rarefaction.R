#Needs input files from twisst_to_bed.R
library(iNEXT)

incfeq_list<-list(as.incfreq(total[4:length(total[1,])]),as.incfreq(total_bf[,4:length(total_bf[1,])]),as.incfreq(total_gf[,4:length(total_gf[1,])]),as.incfreq(total_sf[,4:length(total_sf[1,])]))
#all together (something not working)
out_all<-iNEXT(incfeq_list[[1]],datatype = "incidence_freq",endpoint = 75)
#bf
out_bf<-iNEXT(incfeq_list[[2]],datatype = "incidence_freq",endpoint = 30)
#gf
out_gf<-iNEXT(incfeq_list[[3]],datatype = "incidence_freq",endpoint = 30)
#sf
out_sf<-iNEXT(incfeq_list[[4]],datatype = "incidence_freq",endpoint = 30)


#Estimatation of kilch variation if more individuals would be sampled with Chao Estimator
rare_total<-ChaoRichness(incfeq_list[[1]],datatype = "incidence_freq")
rare_bf<-ChaoRichness(incfeq_list[[2]],datatype = "incidence_freq")
rare_gf<-ChaoRichness(incfeq_list[[3]],datatype = "incidence_freq")
rare_sf<-ChaoRichness(incfeq_list[[4]],datatype = "incidence_freq")

library(ggplot2)

dev.off()
par(mfrow=c(2,2))
#Plot rarefaction curve
total_plot<-ggiNEXT(x=out_all,type=1,grey=T) +
  xlab("Number of individuals sampled") + 
  ylab("Introgressed windows") +   
  theme(legend.position="none") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10, face="bold")) +
        ggtitle("Total")
bf_plot<-ggiNEXT(x=out_bf,type=1,grey=T) +
  xlab("Number of individuals sampled") + 
  ylab("Introgressed windows") +
  theme(legend.position="none") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10, face="bold.italic")) +
        ggtitle("C. wartmanni") 
gf_plot<-ggiNEXT(x=out_gf,type=1,grey=T) +
  xlab("Number of individuals sampled") + 
  ylab("Introgressed windows") +
  theme(legend.position="none") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10, face="bold.italic")) +
        ggtitle("C. macrophthalmus") 
sf_plot<-ggiNEXT(x=out_sf,type=1,grey=T) +
  xlab("Number of individuals sampled") + 
  ylab("Introgressed windows") +
  theme(legend.position="none") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=10, face="bold.italic")) +
        ggtitle("C. arenicolus") 

library(magrittr)
library(multipanelfigure)
figure1 <- multi_panel_figure(columns = 2, rows = 2, panel_label_type = "none")
figure1 %<>%
  fill_panel(total_plot, column = 1, row = 1) %<>%
  fill_panel(bf_plot, column = 2, row = 1) %<>%
  fill_panel(gf_plot, column = 1, row = 2) %<>%
  fill_panel(sf_plot, column = 2, row = 2)
figure1


#Calculation of percentatge for total 
#Estimated
rare_total$Estimator/length(total_bf[,1])
#Observed
rare_total$Observed/length(total_bf[,1])



#Bf
#Estimated
rare_bf$Estimator/length(total_bf[,1])
#Observed
rare_bf$Observed/length(total_bf[,1])


#Gf
#Estimated
rare_gf$Estimator/length(total_gf[,1])
#Observed
rare_gf$Observed/length(total_gf[,1])

#Sf
#Estimated
rare_sf$Estimator/length(total_sf[,1])
#Observed
rare_sf$Observed/length(total_sf[,1])

