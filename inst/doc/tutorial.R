## ----import packages, message=FALSE-------------------------------------------
library(pulseTD)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## ----estimateExpression-------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
test_path <- file.path(system.file(package="pulseTD"),'extdata/test1.sorted.bam')
test_path2 <- file.path(system.file(package="pulseTD"),'extdata/test2.sorted.bam')

## ---- eval=FALSE--------------------------------------------------------------
#  rpkmres <- estimateExpression(txdb,c(test_path,test_path2), by='gene')

## -----------------------------------------------------------------------------
data('rpkmres', package='pulseTD')
head(rpkmres$total_exp)
head(rpkmres$pre_exp)

## ----estimateParams, eval=FALSE-----------------------------------------------
#  data('rpkmSim', package='pulseTD')
#  rpkm_TL <- rpkmSim$labexon[1:2,]
#  rpkm_PT <- rpkmSim$totintr[1:2,]
#  rpkm_TT <- rpkmSim$totexon[1:2,]
#  TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
#  tL <- 10
#  pulseRates<- estimateParams(rpkm_TL, rpkm_TT, rpkm_PT,TimeGrid, tL, clusterNumber=1,loopnumber=10)

## ----correctionParams, eval=FALSE---------------------------------------------
#  data('pulseRates', package='pulseTD')
#  pulseRates_correct = correctionParams(pulseRates)
#  pulseRates_correct@fitfailure

## ----getParams----------------------------------------------------------------
data('pulseRates', package='pulseTD')
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
transcription_params = getParams(pulseRates,'transcription')
degradation_params = getParams(pulseRates, 'degradation')
processing_params = getParams(pulseRates, 'processing')
head(transcription_params)
head(degradation_params)
head(processing_params)
transcription_params = getParams(pulseRates, 'transcription', genename=c(1,2,3))
head(transcription_params)
# get pulse Model value
transcription_pulse = pulseModel(as.matrix(transcription_params[1,]), TimeGrid)
degradation_pulse = pulseModel(as.matrix(degradation_params[1,]), TimeGrid)
processing_pulse = pulseModel(as.matrix(processing_params[1,]), TimeGrid)

## ----getRates-----------------------------------------------------------------
data('pulseRates', package='pulseTD')
pulseRates_correct <- correctionParams(pulseRates)
transcription = getRates(pulseRates_correct,'transcription')
degradation = getRates(pulseRates_correct, 'degradation')
processing = getRates(pulseRates_correct, 'processing')
head(transcription)

## ----getRates by solver-------------------------------------------------------
if(length(pulseRates_correct@fitfailure)==0){
  genename=pulseRates_correct@genenames
}else{
  genename=pulseRates_correct@genenames[-pulseRates_correct@fitfailure]
}
data('rpkmSim', package='pulseTD')
simTL <- rpkmSim$labexon[c(1,2,3),]
simPT <- rpkmSim$totintr[c(1,2,3),]
simTT <- rpkmSim$totexon[c(1,2,3),]
trans_factor <- getRates(pulseRates_correct,'transcription', genename=c(1,2,3))
degr_factor <- getRates(pulseRates_correct, 'degradation', genename=c(1,2,3)) /(simTT-simPT)
proc_factor <- getRates(pulseRates_correct,'processing',genename=c(1,2,3))/simPT
head(degr_factor)

## ----getRates, predict--------------------------------------------------------
transcription_pre <- getRates(pulseRates_correct, 'transcription', timevector=seq(0,360, 15))
degradation_pre <- getRates(pulseRates_correct, 'degradation', timevector=seq(0,360, 15))
processing_pre <- getRates(pulseRates_correct, 'processing', timevector <- seq(0,360, 15))
head(degradation_pre)

## ----plotRatest---------------------------------------------------------------
data('pulseRates', package='pulseTD')
plotRates(pulseRates, 15)

## ----plotRatest,any time------------------------------------------------------
plotRates(pulseRates, 15, predict=c(0,360,15))

## ----predictExpression--------------------------------------------------------
data('pulseRates', package='pulseTD')
pulseRates_correct = correctionParams(pulseRates)
TimeGrid = seq(0,180,15)

## ---- eval=FALSE--------------------------------------------------------------
#  preExp = predictExpression(pulseRates_correct, tg=TimeGrid)

## -----------------------------------------------------------------------------
data('preExp', package='pulseTD')
df = data.frame(preExp[['NM_001002011']])
head(df)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

