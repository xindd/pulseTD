library(pulseTD)
# rpkmRep1
load('E:/pluseTD/pulseTD-2/rpkmRep1.RData')
# rpkmRep2
load('E:/pluseTD/pulseTD-2/rpkmRep2.RData')
##############################################

low.exp <- function(data){
  filter_na = na.omit(data)
  ord = order(apply(filter_na, 1, mean))
  return(c(ord[1:200], ord[500:1000]))
}

labexon1 = rpkmRep1$labexon
totexon1 = rpkmRep1$totexon
totintr1 = rpkmRep1$totintr
labexon2 = rpkmRep2$labexon
totexon2 = rpkmRep2$totexon
totintr2 = rpkmRep2$totintr

lowgene = intersect(low.exp(labexon1),low.exp(totexon1))
lowgene = intersect(lowgene, low.exp(totintr1))
lowgene = intersect(lowgene, low.exp(labexon2))
lowgene = intersect(lowgene, low.exp(totexon2))
lowgene = intersect(lowgene, low.exp(totintr2))

##
###############################################################################
t_time <- c(0, 2, 4, 6, 8)
tL <- 1
gnum = 1000
# ratesRep1
load('E:/pluseTD/pulseTD-2/ratesRep1.RData')
# ratesRep2
load('E:/pluseTD/pulseTD-2/ratesRep2.RData')



plotRates(rep1, 614)
plotRates(rep2, '91452')

a1 = getRates(rep1, 'transcription', timevector = t_time)
b1 = getRates(rep1, 'degradation', timevector = t_time)
c1 = getRates(rep1, 'processing', timevector = t_time)

a2 = getRates(rep2, 'transcription', timevector = t_time)
b2 = getRates(rep2, 'degradation', timevector = t_time)
c2 = getRates(rep2, 'processing', timevector = t_time)
samegene = intersect(rep1@genenames, rep2@genenames)


filter.outers <- function(data, fold=1.5){
  filter_na = na.omit(data)
  eval_var=apply(filter_na, 1, var)
  fenwei = quantile(eval_var)
  filter_var <- filter_na[(eval_var < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) &
                             eval_var > (fenwei[2] - (fenwei[4]-fenwei[2])*fold)),]
  return(filter_var)
}

## rep cor ####
library(LSD)
par(mfrow=c(2,3))
ag = intersect(rownames(filter.outers(as.matrix(a1[samegene,]))),rownames(filter.outers(as.matrix(a2[samegene,]))))
heatscatter(as.vector(a1[ag,]),as.vector(a2[ag,]),
            cor=TRUE,method='pearson',main='pulseTD-alpha',
            ylab='rep2', xlab='rep1')
abline(0,1)

bg = intersect(rownames(filter.outers(as.matrix(b1[samegene,]))),rownames(filter.outers(as.matrix(b2[samegene,]))))
heatscatter(as.vector(b1[bg,]),as.vector(b2[bg,]),
            cor=TRUE,method='pearson',main='pulseTD-beta',
            ylab='rep2', xlab='rep1')
abline(0,1)
cg = intersect(rownames(filter.outers(as.matrix(c1[samegene,]))),rownames(filter.outers(as.matrix(c2[samegene,]))))
heatscatter(as.vector(c1[cg,]),as.vector(c2[cg,]),
            cor=TRUE,method='pearson',main='pulseTD-gamma',
            ylab='rep2', xlab='rep1')
abline(0,1)

### 不过滤
filter.outers <- function(data, fold=0){
  filter_na = na.omit(data)
  eval_var=apply(filter_na, 1, var)
  fenwei = quantile(eval_var)
  filter_var <- filter_na[(eval_var < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) &
                             eval_var > (fenwei[2] - (fenwei[4]-fenwei[2])*fold)),]
  return(data)
}
par(mfrow=c(2,3))
ag = intersect(rownames(filter.outers(as.matrix(a1[samegene,]))),rownames(filter.outers(as.matrix(a2[samegene,]))))
heatscatter(as.vector(a1[ag,]),as.vector(a2[ag,]),
            cor=TRUE,method='pearson',main='pulseTD-alpha',
            ylab='rep2', xlab='rep1',
            xlim = c(0,600), ylim = c(0,600))
abline(0,1)

cg = intersect(rownames(filter.outers(as.matrix(c1[samegene,]))),rownames(filter.outers(as.matrix(c2[samegene,]))))
heatscatter(as.vector(c1[cg,]),as.vector(c2[cg,]),
            cor=TRUE,method='pearson',main='pulseTD-gamma',
            ylab='rep2', xlab='rep1',
            xlim = c(0,400), ylim = c(0,400))
abline(0,1)
bg = intersect(rownames(filter.outers(as.matrix(b1[samegene,]))),rownames(filter.outers(as.matrix(b2[samegene,]))))
heatscatter(as.vector(b1[bg,]),as.vector(b2[bg,]),
            cor=TRUE,method='pearson',main='pulseTD-beta',
            ylab='rep2', xlab='rep1',
            xlim = c(0,200), ylim = c(0,200))
abline(0,1)
###
allgene = intersect(ag, bg)
allgene = intersect(allgene, cg)
heatscatter(c(as.vector(a1[allgene,]), as.vector(c1[allgene,]), as.vector(b1[allgene,])),
            c(as.vector(a2[allgene,]), as.vector(c2[allgene,]), as.vector(b2[allgene,])),
            cor=TRUE,method='pearson',main='Transcriptional dynamics', ylab='replication2', xlab='replication1',
            xlim = c(0,30), ylim = c(0,30))
abline(0,1)
## rep inspect ###
# rpkmRep1
#load('E:/pluseTD/pluseTD/data/rpkmRep1.RData')
# rpkmRep2
#load('E:/pluseTD/pluseTD/data/rpkmRep2.RData')
library(INSPEcT)
t_time <- c(0, 2, 4, 6, 8)
tL <- 1
gnum = 970
names_g = rownames(rpkmRep1$labexon[1:gnum,])
mycerIds1 <- newINSPEcT(t_time, tL,
                        as.data.frame(rpkmRep1$labexon[names_g,]),
                        as.data.frame(rpkmRep1$totexon[names_g,]),
                        as.data.frame(rpkmRep1$labintr[names_g,]),
                        as.data.frame(rpkmRep1$totintr[names_g,]), BPPARAM=SerialParam())

myat1 = as.matrix(ratesFirstGuess(mycerIds1, 'synthesis'))
mybt1 = as.matrix(ratesFirstGuess(mycerIds1, 'degradation'))
myct1 = as.matrix(ratesFirstGuess(mycerIds1, 'processing'))

mycerIds2 <- newINSPEcT(t_time, tL,
                        as.data.frame(rpkmRep2$labexon[names_g,]),
                        as.data.frame(rpkmRep2$totexon[names_g,]),
                        as.data.frame(rpkmRep2$labintr[names_g,]),
                        as.data.frame(rpkmRep2$totintr[names_g,]), BPPARAM=SerialParam())

myat2 = as.matrix(ratesFirstGuess(mycerIds2, 'synthesis'))
mybt2 = as.matrix(ratesFirstGuess(mycerIds2, 'degradation'))
myct2 = as.matrix(ratesFirstGuess(mycerIds2, 'processing'))

names_same = intersect(rownames(myat1), rownames(myat2))

myat1 = myat1[names_same,]
mybt1 = mybt1[names_same,]
myct1 = myct1[names_same,]
myat2 = myat2[names_same,]
mybt2 = mybt2[names_same,]
myct2 = myct2[names_same,]

library(LSD)

ag = intersect(rownames(filter.outers(as.matrix(myat1))),rownames(filter.outers(as.matrix(myat2))))
heatscatter(as.vector(myat1[ag,]),as.vector(myat2[ag,]),
            cor=TRUE,method='pearson',main='INSPEcT-alpha',
            ylab='rep2', xlab='rep1')
abline(0,1)

cg = intersect(rownames(filter.outers(as.matrix(myct1))),rownames(filter.outers(as.matrix(myct2))))
heatscatter(as.vector(myct1[cg,]),as.vector(myct2[cg,]),
            cor=TRUE,method='pearson',main='INSPEcT-gamma',
            ylab='rep2', xlab='rep1')
abline(0,1)
bg = intersect(rownames(filter.outers(as.matrix(mybt1))),rownames(filter.outers(as.matrix(mybt2))))
heatscatter(as.vector(mybt1[bg,]),as.vector(mybt2[bg,]),
            cor=TRUE,method='pearson',main='INSPEcT-beta',
            ylab='rep2', xlab='rep1')
abline(0,1)

### 不过滤
par(mfrow=c(1,1))
ag = intersect(rownames(filter.outers(as.matrix(myat1))),rownames(filter.outers(as.matrix(myat2))))
heatscatter(as.vector(myat1[ag,]),as.vector(myat2[ag,]),
            cor=TRUE,method='pearson',main='INSPEcT-alpha',
            ylab='rep2', xlab='rep1',
            xlim = c(0,600), ylim = c(0,600))
abline(0,1)

bg = intersect(rownames(filter.outers(as.matrix(mybt1))),rownames(filter.outers(as.matrix(mybt2))))
heatscatter(as.vector(mybt1[bg,]),as.vector(mybt2[bg,]),
            cor=TRUE,method='pearson',main='INSPEcT-beta',
            ylab='rep2', xlab='rep1',
            xlim = c(0,200), ylim = c(0,200))
abline(0,1)

cg = intersect(rownames(filter.outers(as.matrix(myct1))),rownames(filter.outers(as.matrix(myct2))))
heatscatter(as.vector(myct1[cg,]),as.vector(myct2[cg,]),
            cor=TRUE,method='pearson',main='INSPEcT-gamma',
            ylab='rep2', xlab='rep1',
            xlim = c(0,1000), ylim = c(0,1000))
abline(0,1)
###
allgene = intersect(ag, bg)
allgene = intersect(allgene, cg)
heatscatter(c(as.vector(myat1[allgene,]), as.vector(myct1[allgene,]), as.vector(mybt1[allgene,])),
            c(as.vector(myat2[allgene,]), as.vector(myct2[allgene,]), as.vector(mybt2[allgene,])),
            cor=TRUE,method='pearson',main='Transcriptional dynamics', ylab='replication2', xlab='replication1'
            )
abline(0,1)




#########################
library(LSD)
par(mfrow=c(2,3))
lg = intersect(lowgene, samegene)

ag = intersect(rownames(a1[lg,]),rownames(a2[lg,]))
heatscatter(as.vector(a1[ag,]),as.vector(a2[ag,]),
            cor=TRUE,method='pearson',main='pulseTD-alpha',
            ylab='rep2', xlab='rep1')
abline(0,1)

cg = intersect(rownames(c1[lg,]),rownames(c2[lg,]))
heatscatter(as.vector(c1[cg,]),as.vector(c2[cg,]),
            cor=TRUE,method='pearson',main='pulseTD-gamma',
            ylab='rep2', xlab='rep1')
abline(0,1)

bg = intersect(rownames(b1[lg,]),rownames(b2[lg,]))
heatscatter(as.vector(b1[bg,]),as.vector(b2[bg,]),
            cor=TRUE,method='pearson',main='pulseTD-beta',
            ylab='rep2', xlab='rep1')
abline(0,1)


##
lg = intersect(lowgene, rownames(myat1))

ag = intersect(rownames(myat1[lg,]),rownames(myat2[lg,]))
heatscatter(as.vector(myat1[ag,]),as.vector(myat2[ag,]),
            cor=TRUE,method='pearson',main='INSPEcT-alpha',
            ylab='rep2', xlab='rep1')
abline(0,1)

cg = intersect(rownames(myct1[lg,]),rownames(myct2[lg,]))
heatscatter(as.vector(myct1[cg,]),as.vector(myct2[cg,]),
            cor=TRUE,method='pearson',main='INSPEcT-gamma',
            ylab='rep2', xlab='rep1')
abline(0,1)
bg = intersect(rownames(mybt1[lg,]),rownames(mybt2[lg,]))
heatscatter(as.vector(mybt1[bg,]),as.vector(mybt2[bg,]),
            cor=TRUE,method='pearson',main='INSPEcT-beta',
            ylab='rep2', xlab='rep1')
abline(0,1)

#########################
### rep distribution####
ff = data.frame(a1= apply(a1[ag,], 1, mean),
                a2= apply(a2[ag,], 1, mean))
bff = data.frame(b1=apply(b1[bg,], 1, mean),
                 b2=apply(b2[bg,], 1, mean))
cff = data.frame(c1=apply(c1[cg,], 1, mean),
                 c2=apply(c2[cg,], 1, mean))

library(ggplot2)
library(reshape2)
ggplot(melt(ff), aes(x = value, color = variable)) +
  geom_density(size=1.5)+
  scale_color_manual(values=c("#555180","#EE7541"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(title="rep1-rep2 transcription distribution", x='rates', fill="rates")+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
        axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))


ggplot(melt(bff), aes(x = value, color = variable)) +
  geom_density(size=1.5)+
  scale_color_manual(values=c("#555180","#EE7541"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(title="rep1-rep2 degradation distributionn", x='rates', fill="rates")+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
        axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))



ggplot(melt(cff), aes(x = value, color = variable)) +
  geom_density(size=1.5)+
  geom_density(size=1.5)+
  scale_color_manual(values=c("#555180","#EE7541"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  labs(title="rep1-rep2 processing distributionn", x='rates', fill="rates")+
  theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+
  theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
        axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
  theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
        axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

