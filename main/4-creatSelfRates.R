library(pulseTD)
# rpkmTRUE

load('E:\\pluseTD\\pulseTD-2\\rpkmTRUE.RData')
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
labexon = rpkmRep$labexon
totexon = rpkmRep$totexon
totintr = rpkmRep$totintr
##
load('E:\\pluseTD\\pulseTD-2\\rpkmSim.RData')
##
########################
gnum = 1000
sim_TL_sample = rpkmSim$labexon[1:gnum, ]
sim_PT_sample = rpkmSim$totintr[1:gnum, ]
sim_TT_sample = rpkmSim$totexon[1:gnum, ]

pulseRates = estimateParams(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, loopnumber=4)
# = correctionParams(pluseRates)
save(pulseRates, file = "pulseRates.RData")
# load(file.path(system.file(package="pluseTD"),'data','pluseRates.RData'))
plotRates(pulseRates, 10, predict=c(0,300,20))
predictExpression(pulseRates, TimeGrid)

#######################
a1 = getRates(pulseRates, 'transcription', timevector = TimeGrid)
b1 = getRates(pulseRates, 'degradation', timevector = TimeGrid)/as.matrix(sim_TT_sample[pulseRates@genenames,]-sim_PT_sample[pulseRates@genenames,])
c1 = getRates(pulseRates, 'processing', timevector = TimeGrid)/as.matrix(sim_PT_sample[pulseRates@genenames,])

load('E:/pluseTD/pulseTD-2/ratesTRUE.RData')
a = ratesTrue$a
b = ratesTrue$b
c = ratesTrue$c
par(mfrow = c(1, 3))
scatterdata = data.frame(pluse_a= as.vector(a1),pluse_b= as.vector(b1),pluse_c= as.vector(c1),
                         same_a = as.vector(a),same_b = as.vector(b),same_c = as.vector(c))

scatterdata = log(filter.outers(scatterdata))
heatscatter((as.vector(scatterdata[,1])), (as.vector(as.matrix(scatterdata[,4]))),
            cor=TRUE,method='pearson',main='Transcription-rates',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-5,1), xlim=c(-10,3))
abline(0,1)
heatscatter((as.vector(scatterdata[,2])), (as.vector(as.matrix(scatterdata[,5]))),
            cor=TRUE,method='pearson',main='degradation',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-5,1), xlim=c(-10,3))
abline(0,1)
heatscatter((as.vector(scatterdata[,3])), (as.vector(as.matrix(scatterdata[,6]))),
            cor=TRUE,method='pearson',main='processing',
            ylab='log-trueRates', xlab='log-modelRates',ylim=c(-5,1), xlim=c(-10,3))
abline(0,1)
#######################

simgene = pulseRates@genenames
###########
pluse5 = estimateParams(sim_TL_sample[simgene,1:5],
                        sim_TT_sample[simgene,1:5],
                        sim_PT_sample[simgene,1:5], TimeGrid[1:5], tL, loopnumber=4)
#pluse5 = correctionParams(pluse5)
save(pluse5, file = "pluse5.RData")
pluse7 =  estimateParams(sim_TL_sample[simgene,1:7],
                         sim_TT_sample[simgene,1:7],
                         sim_PT_sample[simgene,1:7], TimeGrid[1:7], tL, loopnumber=4)
#pluse7 = correctionParams(pluse7)
save(pluse7, file = "pluse7.RData")
pluse9 =  estimateParams(sim_TL_sample[simgene,1:9],
                         sim_TT_sample[simgene,1:9],
                         sim_PT_sample[simgene,1:9], TimeGrid[1:9], tL, loopnumber=4)
#pluse9 = correctionParams(pluse9)
save(pluse9, file = "pluse9.RData")
pluse11= estimateParams(sim_TL_sample[simgene,1:11],
                        sim_TT_sample[simgene,1:11],
                        sim_PT_sample[simgene,1:11], TimeGrid[1:11], tL, loopnumber=4)
#pluse11 = correctionParams(pluse11)
save(pluse11, file = "pluse11.RData")

