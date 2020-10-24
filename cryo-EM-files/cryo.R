library(latex2exp)


k=read.table("cryo-em-data.csv",sep=",",header=T)
cos_dist<-NULL
ang_dist<-NULL

for(i in 1:nrow(k)){
xc<-k[i,1]
yc<-k[i,2]
xi<-k[i,3]
yi<-k[i,4]
xo<-k[i,5]
yo<-k[i,6]

x_phc<-(k[i,3]+k[i,5])/2
y_phc<-(k[i,4]+k[i,6])/2


mem_vector<-c(x_phc-xc,y_phc-yc)
mem_vector<-mem_vector/sqrt((x_phc-xc)**2+(y_phc-yc)**2)

phd_vector<-c(xo-xi,yo-yi)
phd_vector<-phd_vector/sqrt((xo-xi)**2+(yo-yi)**2)
cos_dist<-c(cos_dist,(mem_vector[1]*phd_vector[1])+(mem_vector[2]*phd_vector[2]))
angle<-acos((mem_vector[1]*phd_vector[1])+(mem_vector[2]*phd_vector[2]))
ang_dist<-c(ang_dist,angle)
}

p=hist(ang_dist*(180/3.14),plot=F)



hist(ang_dist*(180/3.14),col="darkgreen",xlab="Angle (deg)",ylab="Frequency", main="",probability=T)

k=cbind(k,ang_dist*(180/3.14))
colnames(k)[8]="Angle(deg)"

tiff("cryo-em-freq-dist.tiff", width = 5, height = 4, units = 'in', res = 300)
hist(ang_dist*(180/3.14),col="darkgreen",xlab="Angle (deg)",ylab="Frequency", main="")
dev.off()


tiff("cryo-em-pdist.tiff", width = 5, height = 4, units = 'in', res = 300)
plot(0,type="n",xlim=c(0,60),ylim=c(0,0.3),xlab=TeX('$\\theta$ (Deg)'),ylab="probability")

for(i in 1:6){
rect(xleft=p$breaks[i],xright=p$breaks[i+1],ybottom=0,ytop=p$counts[i]/sum(p$counts),col="darkblue")
}
dev.off()


