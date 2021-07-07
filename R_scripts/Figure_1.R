#Define color for each species
dev.off()
library(RColorBrewer)
colors<-brewer.pal(12, name="Paired")
greys<-brewer.pal(9, name="Greys")


#Figure 1
pdf("Figure_1.pdf",height=7.5,width=10.8)
#layout(rbind(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,4),c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)))
layout(rbind(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3),c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)))

#Panel a - phosphorus curve
library(readr)
library(circlize)
p <- read_delim("P_data.csv",  ";", escape_double = FALSE, col_names = F,   trim_ws = TRUE)
p<-as.data.frame(p)
x=c(1930,p[,1],2021)
y=c(NA,p[,2],NA)
polyCurve <- function(x, y, from, to, n = 50, miny,col = "red", border = col) {
  drawPoly <- function(fun, from, to, n = 50, miny, col, border) {
    Sq <- seq(from = from, to = to, length = n)
    polygon(x = c(Sq[1], Sq, Sq[n]),
            y = c(miny, fun(Sq), miny),
            col = col, border = border)
  }
  lf <- length(from)
  stopifnot(identical(lf, length(to)))
  if(length(col) != lf)
    col <- rep(col, length.out = lf)
  if(length(border) != lf)
    border <- rep(border, length.out = lf)
  if(missing(miny))
    miny <- min(y)
  interp <- approxfun(x = x, y = y)
  mapply(drawPoly, from = from, to = to, col = col, border = border,
         MoreArgs = list(fun = interp, n = n, miny = miny))
  invisible()
}
plot(x,y,ylim=c(0,130),xlab="Year",ylab=" ",type="l",col="cornflowerblue",xlim=c(1930,2018),yaxt="n",xaxt="n")
polyCurve(c(1949,p[1:length(p[,1]),1],2020),c(0,p[1:length(p[,1]),2],0),n=1000,from=1950,to=2019,col=add_transparency("cornflowerblue",0.8),border=NA)
mtext(side = 2, text = expression(paste("P"[tot]," [mg/",m^3,"]" )), line = 2.75,cex=0.7)
axis(side=1,labels=c(seq(1930,2010,10),2019),at=c(seq(1930,2010,10),2019))
axis(side=2,labels=seq(0,100,10),at=seq(0,100,10))
text(labels="Eutrophication",x=1980,y=20,col="cornflowerblue")



lines(rbind(c(1935,-10),c(1935,120)),type="l",lty=3)
lines(rbind(c(1936,-10),c(1936,125)),type="l",lty=3)
lines(rbind(c(1937,-10),c(1937,130)),type="l",lty=3)
lines(rbind(c(1946,-10),c(1946,120)),type="l",lty=3)
lines(rbind(c(2015,-10),c(2015,120)),type="l",lty=3)

points(1946,120,xpd=T,pch=3,col=colors[1],lwd=3,cex=1.5)
points(1935,120,xpd=T,pch=3,col=colors[3],lwd=3,cex=1.5)
points(1936,125,xpd=T,pch=3,col=colors[7],lwd=3,cex=1.5)
points(1937,130,xpd=T,pch=3,col=greys[3],lwd=3,cex=1.5)
points(2015,131,xpd=T,pch=19,col=colors[8],lwd=3,cex=1.5)
points(2015,125.5,xpd=T,pch=19,col=colors[2],lwd=3,cex=1.5)
points(2015,120,xpd=T,pch=19,col=colors[4],lwd=3,cex=1.5)


#Panel b - PCA
library(readr)
PCA_file_covMat <- read_delim("PCAngsd_output.cov",  ",", escape_double = FALSE, col_names = FALSE,   trim_ws = TRUE)
m<-as.matrix(PCA_file_covMat)
e<-eigen(m)
cols<-c(rep(greys[3],11),rep(colors[7],3),rep(colors[3],2),rep(colors[1],2),rep(colors[8],5),rep(colors[4],3),rep(colors[2],6))
plot(main=" ",rbind(e$vectors[1:31,c(1,2)],jitter(e$vectors[32,c(1,2)],factor=0.1)),ylab=paste0("PC2 (",round(e$values[2]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC1 (",round(e$values[1]/sum(e$values)*100,digits=2),"%)"),col=cols,pch=c(rep(3,18),rep(19,14)),lwd=3,cex=1.5)


#Panel c - Legend
plot.new()
#plot.new()

l_text<-c(expression(paste(italic("C. gutturosus")," (pre)")),expression(paste(italic("C. arenicolus")," (pre)")),expression(paste(italic("C. arenicolus")," (post)")),expression(paste(italic("C. wartmanni")," (pre)")),expression(paste(italic("C. wartmanni")," (post)")),expression(paste(italic("C. macrophthalmus")," (pre)")),expression(paste(italic("C. macrophthalmus")," (post)")))
l_text_2<-rep(" ",length(l_text))
fake<-rep(" ",length(l_text))
#legend(x="topleft",inset=c(-0.425,-0.225),xpd=T,bty="n",title=expression(bold("")),legend=l_text,col=c(greys[3],colors[7],colors[8],colors[1],colors[2],colors[3],colors[4]),pch=c(3,3,NA,3,NA,3,NA),y.intersp=3,pt.lwd=3,pt.cex=1.5)
#legend(x="topleft",inset=c(-0.425,-0.225),xpd=T,bty="n",title=expression(bold("")),legend=l_text_2,col=c(greys[3],colors[7],colors[8],colors[1],colors[2],colors[3],colors[4]),pch=c(NA,NA,19,NA, 19,NA,19),y.intersp=3,pt.lwd=3,pt.cex=1.5)

legend(x="topleft",inset=c(-0.1,-0.225),xpd=T,bty="n",title=expression(bold("")),legend=l_text,col=c(greys[3],colors[7],colors[8],colors[1],colors[2],colors[3],colors[4]),pch=c(3,3,NA,3,NA,3,NA),y.intersp=3,pt.lwd=3,pt.cex=1.5)
legend(x="topleft",inset=c(-0.1,-0.225),xpd=T,bty="n",title=expression(bold("")),legend=l_text_2,col=c(greys[3],colors[7],colors[8],colors[1],colors[2],colors[3],colors[4]),pch=c(NA,NA,19,NA, 19,NA,19),y.intersp=3,pt.lwd=3,pt.cex=1.5)


#Panel d - admixture plot
#plot.new()
admix = read.table( "PCAngsd_output.qopt",sep=",")
sorted<-rbind(admix[1:14,],admix[19:23,],admix[15:16,],admix[24:26,],admix[17:18,],admix[27:28,],admix[c(29:32),])
barplot(height=t(sorted),col=c(greys[3],colors[4],colors[2],colors[8]),ylab="Admixture proportion",las=2,space=c(rep(0.1,11),1,rep(0.1,2),0.3,rep(0.1,4),1,0.1,0.3,0.1,0.1,1,0.1,0.3,rep(0.1,5)))
help(barplot)
dev.off()


