library(pulseTD)

## 1 expression
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
test_path <- file.path(system.file(package="pulseTD"),'data','test1.sorted.bam')
test_path2 <- file.path(system.file(package="pulseTD"),'data','test2.sorted.bam')
res <- estimateExpression(txdb,c(test_path,test_path2), by='gene')
head(res$total_exp)
head(res$pre_exp)

## 2 estimateParams
load(file.path(system.file(package="pulseTD"),'data','rpkmSim.RData'))
sim_TL_sample <- rpkmSim$labexon[1:20,]
sim_PT_sample <- rpkmSim$totintr[1:20,]
sim_TT_sample <- rpkmSim$totexon[1:20,]
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
pulseRates <- estimateParams(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, loopnumber=50)

## 3 correction params
load(file.path(system.file(package="pluseTD"),'data','pluseRates.RData'))
pulseRates_correct <- correctionParams(pulseRates)

## 4 getParams
load(file.path(system.file(package="pulseTD"),"data","pulseRates.RData"))
pulseRates_correct <- correctionParams(pulseRates)
transcription_params = getParams(pulseRates,'transcription')
degradation_params = getParams(pulseRates, 'degradation')
processing_params = getParams(pulseRates, 'processing')
head(transcription_params)
head(degradation_params)
head(processing_params)
transcription_params = getParams(pulseRates, 'transcription', genename=c(1,2,3))
head(transcription_params)
###
transcription_pulse = pulseModel(as.matrix(transcription_params[1,]), TimeGrid)
degradation_pulse = pulseModel(as.matrix(degradation_params[1,]), TimeGrid)
processing_pulse = pulseModel(as.matrix(processing_params[1,]), TimeGrid)

## 5 get Rates
load(file.path(system.file(package="pulseTD"),"data","pulseRates.RData"))
pulseRates_correct <- correctionParams(pulseRates)
transcription = getRates(pulseRates_correct,'transcription')
degradation = getRates(pulseRates_correct, 'degradation')
processing = getRates(pulseRates_correct, 'processing')
head(transcription)
head(degradation)
head(processing)
#### Scale Factor
if(length(pulseRates_correct@fitfailure)==0){
  genename=pulseRates_correct@genenames
}else{
  genename=pulseRates_correct@genenames[-pulseRates_correct@fitfailure]
}

transcription_factorr = getRates(pulseRates_correct,'transcription')
degradation_factor = getRates(pulseRates_correct, 'degradation')/(sim_TT_sample[genename,]-sim_PT_sample[genename,])
processing_factor = getRates(pulseRates_correct, 'processing')/(sim_PT_sample[genename,])
head(transcription_factor)
head(degradation_factor)
head(processing_factor)

## 6 predict Rates
transcription_pre = getRates(pulseRates_correct, 'transcription', timevector=seq(0,360, 15))
degradation_pre = getRates(pulseRates_correct, 'degradation', timevector=seq(0,360, 15))
processing_pre = getRates(pulseRates_correct, 'processing', timevector=seq(0,360, 15))
head(transcription_pre)
head(degradation_pre)
head(processing_pre)

## 7 Select specific genes calculation rate
transcription_genes = getRates(pulseRates_correct, 'transcription', genename=c(1,2,3))
degradation_genes = getRates(pulseRates_correct, 'degradation', genename=c(1,2,3))
processing_genes = getRates(pulseRates_correct, 'processing', genename=c(1,2,3))
head(transcription_genes)
head(degradation_genes)
head(processing_genes)

## 8 plot a gene
load(file.path(system.file(package="pulseTD"),"data","pulseRates.RData"))
pulseRates_correct <- correctionParams(pulseRates)
plotRates(pulseRates_correct, 15)

## 9 plot a gene with precict
plotRates(pulseRates, 15, predict=c(0,180,20))

## 10 predict expression
load(file.path(system.file(package="pulseTD"),'data','pulseRates.RData'))
pulseRates_correct = correctionParams(pulseRates)
TimeGrid = seq(0,180,15)
preExp = predictExpression(pulseRates_correct, tg=TimeGrid)
df = data.frame(preExp[['NM_001002011']])
head(df)

###
load(file.path(system.file(package="pulseTD"),'data','pulse9.RData'))
pulse9 <- correctionParams(pulse9)
preExp <- predictExpression(pulse9, end=180,interval=15)

TimeGrid <- seq(0,180,15)
df <- data.frame(preExp[['NM_001002011']])
head(df)
df$time <- seq(0,180,15)
idx <- 1:9
library(ggplot2)
library(grid)
g1 <- ggplot(df)+
  geom_ribbon(aes(x=time,ymin=downPT,ymax=upPT), fill="grey",alpha=0.4)+
  geom_line(aes(x=time, y=PT))+
  annotate("point",x=TimeGrid, y=as.vector(as.matrix(sim_PT_sample['NM_001002011',])), color='red')+
  annotate("point",x=TimeGrid[idx], y=as.vector(as.matrix(sim_PT_sample['NM_001002011',idx])))
g2 <- ggplot(df)+
  geom_ribbon(aes(x=time,ymin=downTT,ymax=upTT), fill="grey",alpha=0.4)+
  geom_line(aes(x=time, y=TT))+
  annotate("point",x=TimeGrid, y=as.vector(as.matrix(sim_TT_sample['NM_001002011',])))
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
popViewport()











