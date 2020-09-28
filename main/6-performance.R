dir = 'E:\\pluseTD\\pulseTD-2'
load(file.path(dir,'rpkmTRUE.RData'))
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
labexon = rpkmRep$labexon
totexon = rpkmRep$totexon
totintr = rpkmRep$totintr
##
load(file.path(dir,'rpkmSim.RData'))
##
########################
gnum = 1000
sim_TL_sample = rpkmSim$labexon[1:gnum, ]
sim_PT_sample = rpkmSim$totintr[1:gnum, ]
sim_TT_sample = rpkmSim$totexon[1:gnum, ]

####TrueData--SimData--cor######
genename = rownames(sim_TL_sample)
par(mfrow = c(1, 1))
library(LSD)
heatscatter(log(as.vector(sim_TL_sample)),log(as.vector(as.matrix(labexon[genename,]))),
            cor=TRUE,method='pearson',main='LabelExons', ylab='log(TrueData)', xlab='log(SimData)',
            ylim=c(-5,10), xlim=c(-5,10))
abline(0,1)

heatscatter(log(as.vector(sim_TT_sample)), log(as.vector(as.matrix(totexon[genename,]))),
            cor=TRUE,method='pearson',main='TotalExons', ylab='log(TrueData)', xlab='log(SimData)',ylim=c(-5,15), xlim=c(-5,15))
abline(0,1)

heatscatter(log(as.vector(sim_PT_sample)), log(as.vector(as.matrix(totintr[genename,]))),
            cor=TRUE,method='pearson',main='TotalIntrons', ylab='log(TrueData)', xlab='log(SimData)',ylim=c(-5,6), xlim=c(-5,6))
abline(0,1)

####self--accuracy##############
#自身准确性
load(file.path(dir,'pulseRates.RData'))
load(file.path(dir,'ratesTrue.RData'))

at = ratesTrue$a[genename,]
bt = ratesTrue$b[genename,]
ct = ratesTrue$c[genename,]

genename=pulseRates@genenames
sim_a = as.matrix(getRates(pulseRates, 'transcription'))[genename,]
sim_b = as.matrix(getRates(pulseRates, 'degradation')[genename,]/(sim_TT_sample[genename,]-sim_PT_sample[genename,]))
sim_c = as.matrix(getRates(pulseRates, 'processing')[genename,]/sim_PT_sample[genename,])
same_a = as.matrix(at[genename,])
same_b = as.matrix(bt[genename,])
same_c = as.matrix(ct[genename,])
### 过滤
filter.outers <- function(data, fold=1.5){
  filter_na = na.omit(data)
  eval_var=apply(filter_na, 1, var)
  eval_var = na.omit(eval_var)
  fenwei = quantile(eval_var)
  filter_var <- filter_na[(eval_var < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) &
                             eval_var > (fenwei[2] - (fenwei[4]-fenwei[2])*fold)),]
  return(filter_var)
}
par(mfrow = c(1, 1))
scatterdata = data.frame(pluse_a= as.vector(sim_a),pluse_b= as.vector(sim_b),pluse_c= as.vector(sim_c),
                         same_a = as.vector(same_a),same_b = as.vector(same_b),same_c = as.vector(same_c))

scatterdata = log(filter.outers(scatterdata))

heatscatter((as.vector(scatterdata[,1])), (as.vector(as.matrix(scatterdata[,4]))),
            cor=TRUE,method='pearson',main='Transcription-rates',
            ylab='log-trueRates', xlab='log-modelRates',
            ylim=c(-5,0.5), xlim=c(-5,0.5))
abline(0,1)
scatterdata = filter.outers(scatterdata, 2)
heatscatter((as.vector(scatterdata[,2])), (as.vector(as.matrix(scatterdata[,5]))),
            cor=TRUE,method='pearson',main='degradation',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-5,-1), xlim=c(-5,-1))
abline(0,1)
heatscatter((as.vector(scatterdata[,3])), (as.vector(as.matrix(scatterdata[,6]))),
            cor=TRUE,method='pearson',main='processing',
            ylab='log-trueRates', xlab='log-modelRates'
            ,ylim=c(-5,1), xlim=c(-5,1))
abline(0,1)

load(file.path(dir,'solver.RData'))
solver_a = solver$a[rownames(same_a),]
solver_b = solver$b[rownames(same_a),]
solver_c = solver$c[rownames(same_a),]

scatterdata = data.frame(solver_a= as.vector(solver_a),solver_b= as.vector(solver_b),pluse_c= as.vector(solver_c),
                         same_a = as.vector(same_a),same_b = as.vector(same_b),same_c = as.vector(same_c))
scatterdata = log(filter.outers(scatterdata))
heatscatter((as.vector(scatterdata[,1])), (as.vector(as.matrix(scatterdata[,4]))),
            cor=TRUE,method='pearson',main='Transcription-rates',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-5,0.5), xlim=c(-5,0.5))
abline(0,1)
heatscatter((as.vector(scatterdata[,2])), (as.vector(as.matrix(scatterdata[,5]))),
            cor=TRUE,method='pearson',main='degradation',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-10,0.5), xlim=c(-10,0.5))
abline(0,1)
heatscatter((as.vector(scatterdata[,3])), (as.vector(as.matrix(scatterdata[,6]))),
            cor=TRUE,method='pearson',main='processing',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-7,0.5), xlim=c(-7,0.5))
abline(0,1)


library(ggplot2)
library(reshape2)
library(grid)
scatterdata = data.frame(pluse_a= as.vector(sim_a),pluse_b= as.vector(sim_b),pluse_c= as.vector(sim_c),
                         same_a = as.vector(same_a),same_b = as.vector(same_b),same_c = as.vector(same_c))
scatterdata = log(filter.outers(scatterdata))

g1 = ggplot(log(scatterdata[,c(1,4)]), aes(x=pluse_a, y=same_a))+geom_point()+geom_abline(intercept=0,slope=1)+ggtitle('A')
g2 = ggplot(log(scatterdata[,c(2,5)]), aes(x=pluse_b, y=same_b))+geom_point()+geom_abline(intercept=0,slope=1)+ggtitle('B')
g3 = ggplot(log(scatterdata[,c(3,6)]), aes(x=pluse_c, y=same_c))+geom_point()+geom_abline(intercept=0,slope=1)+ggtitle('C')

grid.newpage() # 生成新画布
pushViewport(viewport(layout = grid.layout(1,3))) # 指定一个1乘2的layout
print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # 画
print(g2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) # 画
print(g3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3)) # 画
popViewport() # 退出当前位置

####Compared--accuracy##############
load(file.path(dir,'reg2.RData'))
load(file.path(dir,'reg3.RData'))
load(file.path(dir,'reg4.RData'))
load(file.path(dir,'solver.RData'))


pluse_error_a =  abs(as.vector(filter.outers(sim_a - same_a)))
pluse_error_b =  abs(as.vector(filter.outers(sim_b - same_b)))
pluse_error_c =  abs(as.vector(filter.outers(sim_c - same_c)))
reg_error_a2 =   abs(as.vector(filter.outers(reg2$a - at)))
reg_error_b2 =   abs(as.vector(filter.outers(reg2$b - bt)))
reg_error_c2 =   abs(as.vector(filter.outers(reg2$c - ct)))
reg_error_a3 =   abs(as.vector(filter.outers(reg3$a - at)))
reg_error_b3 =   abs(as.vector(filter.outers(reg3$b - bt)))
reg_error_c3 =   abs(as.vector(filter.outers(reg3$c - ct)))
reg_error_a4 =   abs(as.vector(filter.outers(reg4$a - at)))
reg_error_b4 =   abs(as.vector(filter.outers(reg4$b - bt)))
reg_error_c4 =   abs(as.vector(filter.outers(reg4$c - ct)))
solver_error_a = abs(as.vector(filter.outers(solver$a - at)))
solver_error_b = abs(as.vector(filter.outers(solver$b - bt)))
solver_error_c = abs(as.vector(filter.outers(solver$c - ct)))

bardata = data.frame(pluse=c(mean(pluse_error_a), mean(pluse_error_b), mean(pluse_error_c)),
                     solve=c(mean(solver_error_a),mean(solver_error_b), mean(solver_error_c)),
                     ploy2=c(mean(reg_error_a2),  mean(reg_error_b2), mean(reg_error_c2)),
                     ploy3=c(mean(reg_error_a3),  mean(reg_error_b3), mean(reg_error_c3)),
                     ploy4=c(mean(reg_error_a4),  mean(reg_error_b4), mean(reg_error_c4)),
                     id = c('a', 'b', 'c'))
bardata=melt(bardata, id='id', value=c())

ggplot(bardata, aes(x=variable, y=value, fill=id))+
  #, position = 'dodge'
  geom_bar(stat="identity", position = 'dodge')+
  scale_fill_manual(values=c("#E87E75","#C2909A","#96ADB9"))+
  ylab("error")+xlab("")+
  theme_bw()+labs(fill='')+
  # theme(legend.position='bottom') 标签位置
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))


library(magrittr)
library(ggpubr)
data <- data.frame(x = c("Alpha","Bravo","Charlie","Delta"),y=c(200,20,10,15))
#画下面
p1 <- ggplot(bardata, aes(x=variable, y=value, fill=id)) + geom_bar(stat='identity',position=position_dodge()) +
  labs(x=NULL,y=NULL,fill=NULL)+    #可自定义标签名字
  coord_cartesian(ylim = c(0,500))   #设置下面一半的值域
p2 <- ggplot(bardata, aes(x=variable, y=value, fill=id)) + geom_bar(stat='identity',position=position_dodge()) +
  labs(x=NULL,y=NULL,fill=NULL) +   #不要标签
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +     #去掉X轴和X轴的文字
  coord_cartesian(ylim = c(1000,50000)) +  #设置上面一半的值域
  scale_y_continuous(breaks = c(1000,50000,1000)) #以5为单位划分Y轴
ggarrange(p2,p1,heights=c(2/5, 3/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
#######预测对比########################################################
load(file.path(dir,'error_5.RData'))
load(file.path(dir,'error_7.RData'))
load(file.path(dir,'error_9.RData'))
load(file.path(dir,'error_11.RData'))

meltdf = rbind(error_5,error_7,error_9,error_11)

library(reshape)
library(ggplot2)
meltdf$id=c(rep(5,dim(error_5)[1]),rep(7,dim(error_5)[1]),rep(9,dim(error_5)[1]),rep(11,dim(error_5)[1]))
meltdf=na.omit(meltdf)
head(melt(meltdf[,c(2,12,10,4,6,8,13)], id='id'))
ggplot(melt(log(filter.outers(meltdf[,c(2,4,6,8,13)])), id='id'),
       aes(x = as.factor(variable), y = value, fill=as.factor(id))) +
  geom_boxplot()
###########################################
#######速率误差对比##########
load(file.path(dir,'rates_error_5.RData'))
load(file.path(dir,'rates_error_7.RData'))
load(file.path(dir,'rates_error_9.RData'))
load(file.path(dir,'rates_error_11.RData'))

##a
melta = as.data.frame(rbind(as.matrix(rates_error_5[,c(1,2,3,4)]),
                            as.matrix(rates_error_7[,c(1,2,3,4)]),
                            as.matrix(rates_error_9[,c(1,2,3,4)]),
                            as.matrix(rates_error_11[,c(1,2,3,4)])))
melta$id=c(rep(5,dim(rates_error_5)[1]),
           rep(7,dim(rates_error_7)[1]),
           rep(9,dim(rates_error_9)[1]),
           rep(11,dim(rates_error_11)[1]))
melta = na.omit(melta)

head(melt(melta, id='id'))

ggplot(melt(filter.outers(melta), id='id'),
       aes(x = as.factor(variable), y = log(value), fill=as.factor(id))) +
  geom_boxplot()+
  #scale_fill_manual(values=c("#C6B093", "#EDAF67", "#DDBC62","#C3CA68"))+
  labs(title="Transcription",y = "log-error")+
  xlab("")+
  #theme_bw()+labs(fill='')+
  # theme(legend.position='bottom') 标签位置
  scale_fill_manual(values=c("#EDAF67", "#EE7541", "#C2909A","#96A2B8"))+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

##b
meltb = as.data.frame(rbind(as.matrix(rates_error_5[,c(6,7,8,9)]),
                            as.matrix(rates_error_7[,c(6,7,8,9)]),
                            as.matrix(rates_error_9[,c(6,7,8,9)]),
                            as.matrix(rates_error_11[,c(6,7,8,9)])))
meltb$id=c(rep(5,dim(rates_error_5)[1]),
           rep(7,dim(rates_error_7)[1]),
           rep(9,dim(rates_error_9)[1]),
           rep(11,dim(rates_error_11)[1]))
meltb = na.omit(meltb)
head(melt(meltb, id='id'))

ggplot(melt(filter.outers(meltb), id='id'),
       aes(x = as.factor(variable), y = log(value), fill=as.factor(id))) +
  geom_boxplot()+
  labs(title="Degradation",y = "log-error")+
  xlab("")+
  # theme(legend.position='bottom') 标签位置
  scale_fill_manual(values=c("#EDAF67", "#EE7541", "#C2909A","#96A2B8"))+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
##c
meltc = as.data.frame(rbind(as.matrix(rates_error_5[, c(11,12:14)]),
                            as.matrix(rates_error_7[, c(11,12:14)]),
                            as.matrix(rates_error_9[, c(11,12:14)]),
                            as.matrix(rates_error_11[,c(11,12:14)])))
meltc$id=c(rep(5,dim(rates_error_5)[1]),
           rep(7,dim(rates_error_7)[1]),
           rep(9,dim(rates_error_9)[1]),
           rep(11,dim(rates_error_11)[1]))
meltc = na.omit(meltc)
head(melt(meltc), id='id')
ggplot(melt(filter.outers(log(meltc)), id='id'),
       aes(x = as.factor(variable), y = (value), fill=as.factor(id))) +
  geom_boxplot()+
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","#79A72C"))+
  #theme(legend.position='none',plot.title = element_text(size=25,hjust = 0.5))+
  labs(title="procressing",y = "log-error")+
  xlab("")+
  # theme(legend.position='bottom') 标签位置
  scale_fill_manual(values=c("#EDAF67", "#EE7541", "#C2909A","#96A2B8"))+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
##c


######Predicted expression value##############
library(grid)
singleGene = 'NM_001002011'
load(file.path(system.file(package="pluseTD"),'data','pluseRates.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse5.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse7.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse9.RData'))
load(file.path(system.file(package="pluseTD"),'data','pluse11.RData'))

plotRates(pluseRates, singleGene, predict=c(0,180,15))
preExp = predictExpression(pluse7, 180,1)
df = data.frame(preExp[[singleGene]])
df$time = seq(0,180,1)
idx=1:7
ggplot(df)+
  geom_ribbon(aes(x=time,ymin=downTT,ymax=upTT), fill="grey",alpha=0.4)+
  geom_line(aes(x=time, y=TT, color='lm'), size=1)+
  scale_color_manual(values=c("#0000FF"))+
  annotate("point",x=TimeGrid, y=as.vector(as.matrix(sim_TT_sample[singleGene,])), color='red', size=2)+
  annotate("point",x=TimeGrid[idx], y=as.vector(as.matrix(sim_TT_sample[singleGene,idx])), size=2)+
  labs(title="predict mRNA expression",y = "Expression", x='Time(min)')+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
        axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))


g1 = ggplot(df)+
  geom_ribbon(aes(x=time,ymin=downPT,ymax=upPT), fill="grey",alpha=0.4)+
  geom_line(aes(x=time, y=PT))+
  annotate("point",x=TimeGrid, y=as.vector(as.matrix(sim_PT_sample[singleGene,])), color='red')+
  annotate("point",x=TimeGrid[idx], y=as.vector(as.matrix(sim_PT_sample[singleGene,idx])))
g2=ggplot(df)+
  geom_ribbon(aes(x=time,ymin=downTT,ymax=upTT), fill="grey",alpha=0.4)+
  geom_line(aes(x=time, y=TT))+
  annotate("point",x=TimeGrid, y=as.vector(as.matrix(sim_TT_sample[singleGene,])), color='red')+
  annotate("point",x=TimeGrid[idx], y=as.vector(as.matrix(sim_TT_sample[singleGene,idx])))+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))


grid.newpage() # 生成新画布
pushViewport(viewport(layout = grid.layout(1,2))) # 指定一个1乘2的layout
print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # 画
print(g2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) # 画
popViewport() # 退出当前位置
#########################
######ROC#######
library(pROC)
###
corr_pluse5 = correctionParams(pluse5)
genename5=pluse5@genenames[-c(corr_pluse5@fitfailure)]
roc_5 = ROCdata(corr_pluse5, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                same_a,same_b,same_c,
                genename5, TimeGrid, tL, idx=1:5, idx2 = 6:7)
p5 = apply(pData(pluse5@score.degradation)[genename5,], 1, mean)
#
genename=pluse7@genenames[-c(pluse7@fitfailure)]
roc_7 = ROCdata(pluse7, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                same_a,same_b,same_c,
                genename, TimeGrid, tL, idx=1:7, idx2 = 8:13)
#
load(file.path(system.file(package="pluseTD"),'data','pluse9.RData'))
corr_pluse9 = correctionParams(pluse9)
genename9=pluse9@genenames[-c(corr_pluse9@fitfailure)]
roc_9 = ROCdata(corr_pluse9, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                same_a,same_b,same_c,
                genename9, TimeGrid, tL, idx=1:9, idx2 = 10:13)
#
genename=pluse11@genenames[-c(pluse11@fitfailure)]
roc_11 = ROCdata(pluse11, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                same_a,same_b,same_c,
                genename, TimeGrid, tL, idx=1:11, idx2 = 12:13)


load(file.path(system.file(package="pluseTD"),'data','pluse9.RData'))
corr_pluse9 = correctionParams(pluse9)
genename=corr_pluse9@genenames[-c(corr_pluse9@fitfailure)]
p9 = apply(pData(pluse9@score.degradation)[genename,], 1, sum)
###

par(mfrow=c(1,1))
r5=plot.roc(roc_5$label, cv, col='red', lwd=4)
r7=plot.roc(roc_7$label,roc_7$socre, col='deepskyblue', lwd=4, add=TRUE)
r9=plot.roc(roc_9$label,pData(pluse9@score.degradation)[genename,2], col='navy', lwd=4, add=TRUE)
r11=plot.roc(roc_11$label,roc_11$socre, col='green', lwd=4, add=TRUE)

legend('bottomright',
       legend= paste(c('pluse5', 'pluse7', 'pluse9', 'pluse11'),
                     ' - AUC=',signif(c(as.numeric(r5$auc), as.numeric(r7$auc),
                                        as.numeric(r9$auc),as.numeric(r11$auc)), 4),
                     sep=''),
       col=c('red', 'deepskyblue', 'navy', 'green'), lty=1, lwd=4)

library(pheatmap)
a = pData(corr_pluse5@score.transcription)[genename5,2:6]
a = t(apply(a,1,function(x)(x-min(x))/(max(x)-min(x))))
tmpdf=cbind(apply(a[,1:2], 1, min),
            apply(a[,3:4], 1, min),
            #a[,5],
            roc_5$label)
tmpdf[which(tmpdf>0.5)]=1
tmpdf[which(tmpdf<=0.5)]=0

tmpdf = cbind(tmpdf, cv=apply(tmpdf[,1:3], 1, sum))
tmpdf = tmpdf[order(tmpdf[,7]),]
pheatmap(tmpdf,
         cluster_row=TRUE,
         cluster_col=FALSE,scale="none",legend=TRUE)
tmpdf = as.matrix(tmpdf)
lm(tmpdf[,7]~tmpdf[,2:5])
