library(ggplot2)
library(plyr)

result<-read.csv("par.csv",sep="\t", header = FALSE)
result<-result[1:56092,]
t<-seq(1:56092)
rs<-data.frame(t,result)
colnames(rs)<-c("runs","piv","beta","dEs","dv","E0s","v0","R0")

p<-ggplot() + geom_line(data=rs,aes(x=runs,y=R0)) + 
  xlab("Accepted chain") +
  ylab("Values")+ggtitle("R_0")+xlim(0,10000)+
  theme(plot.title = element_text(hjust = 0.5)) 
p

result<-read.csv("par.csv",sep="\t", header = FALSE)
result<-result[10000:56092,]
colnames(result)<-c("piv","beta","dEs","dv","E0s","v0","R0")
library("HDInterval")
a<-hdi(result$R0)
R0min<-a[[1]]
R0max<-a[[2]]
R0mean<-mean(result$R0)
R0min
R0mean
R0max

