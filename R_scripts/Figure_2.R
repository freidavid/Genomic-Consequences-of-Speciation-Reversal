pdf("./d_stat.pdf",height=5,width=12)
par(mar=c(8.5,5,5,5))
number<-c(1:3,7,11:13,17:18)
#gutturosus into macropthalmus, gutturosus into arenicolus, macrophthalmus into arenicolus
#=deep into shallow
d_t_s <- c(0.018622,0.026814,0.020142)
#arenicolus into macrophthalmus
#=shallow into deep
s_t_d <- c(0.00327)
#gutturosus into wartmanni, macrophtalmus into wartmanni, areniclus into wartmanni
#= bentic into pelagic
b_t_p <- c(0.033183,0.033174,0.028265)
#wartmanni into macrophthalmus, wartmanni into arenicolus
#=pelagic into benthic
p_t_b <- c(0.005573,0.003391)
#make one vectore out of data to plot
data<-c(d_t_s,s_t_d,b_t_p,p_t_b)
#create an empty plot
plot(xlab="",ylab="D",cex=1.5,ylim=c(-0.0025,1.1*max(data)),number,data,col=NA,pch=1,axes=F)
#make axes
axis(side=1,at=number,labels=c("a","b","c","a","a","b","c","a","b"))
axis(side=2)
#abline at 0
abline(h=0,lty=2)
#define colors
col=c("black","black","black","darkgrey","black","black","black","darkgrey","darkgrey")
#draw squares
points(y=data,x=number,col=col,pch=15)
#text above different categories
mtext(side=3,line=1,text="1 deeper into shallower",at=2)
mtext(side=3,line=1,text="2 shallower into deeper",at=7)
mtext(side=3,line=1,text="3 benthic into pelagic",at=12)
mtext(side=3,line=1,text="4 pelagic into benthic",at=17.5)
box()
#error bars
arrows(x0=number[1], y0=data[1]-0.0014, x1=number[1], y1=data[1]+0.0014, code=3, angle=90, length=0.1,lwd=2)
arrows(x0=number[2], y0=data[2]-0.0014, x1=number[2], y1=data[2]+0.0014, code=3, angle=90, length=0.1,lwd=2)
arrows(x0=number[3], y0=data[3]-0.0014, x1=number[3], y1=data[3]+0.0014, code=3, angle=90, length=0.1,lwd=2)
arrows(x0=number[4], y0=data[4]-0.0014, x1=number[4], y1=data[4]+0.0014, code=3, angle=90, length=0.1,col="darkgrey",lwd=2)
arrows(x0=number[5], y0=data[5]-0.0014, x1=number[5], y1=data[5]+0.0014, code=3, angle=90, length=0.1,lwd=2)
arrows(x0=number[6], y0=data[6]-0.0014, x1=number[6], y1=data[6]+0.0014, code=3, angle=90, length=0.1,lwd=2)
arrows(x0=number[7], y0=data[7]-0.0014, x1=number[7], y1=data[7]+0.0014, code=3, angle=90, length=0.1,lwd=2)
arrows(x0=number[8], y0=data[8]-0.0014, x1=number[8], y1=data[8]+0.0014, code=3, angle=90, length=0.1,col="darkgrey",lwd=2)
arrows(x0=number[9], y0=data[9]-0.0014, x1=number[9], y1=data[9]+0.0014, code=3, angle=90, length=0.1,col="darkgrey",lwd=2)


dev.off()

