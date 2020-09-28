#######加载数据###############
dirpath = 'E:\\pluseTD\\RASG\\data'
rpkms_rep1=list()
rpkms_rep1$foursu_exons = read.table(paste(dirpath, '\\expression_gene.4sU.M.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$foursu_introns = read.table(paste(dirpath, '\\expression_gene.4sU.P.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$total_exons = read.table(paste(dirpath, '\\expression_gene.total.M.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$total_introns = read.table(paste(dirpath, '\\expression_gene.total.P.txt',sep=''),sep='\t', header = TRUE)

rownames(rpkms_rep1$foursu_exons) = rpkms_rep1$foursu_exons[,1]
rownames(rpkms_rep1$foursu_introns) = rpkms_rep1$foursu_introns[,1]
rownames(rpkms_rep1$total_exons) = rpkms_rep1$total_exons[,1]
rownames(rpkms_rep1$total_introns) = rpkms_rep1$total_introns[,1]

rpkms_rep1$foursu_exons = rpkms_rep1$foursu_exons[,-1]
rpkms_rep1$foursu_introns = rpkms_rep1$foursu_introns[,-1]
rpkms_rep1$total_exons = rpkms_rep1$total_exons[,-1]
rpkms_rep1$total_introns = rpkms_rep1$total_introns[,-1]

labexon = rpkms_rep1$foursu_exons
labintr = rpkms_rep1$foursu_introns
totexon = rpkms_rep1$total_exons
totintr = rpkms_rep1$total_introns

genelist = intersect(rownames(labexon),rownames(totexon))
genelist = intersect(genelist,rownames(totintr))
genelist = intersect(genelist,rownames(labintr))

labexon = labexon[genelist,]
labintr = labintr[genelist,]
totexon = totexon[genelist,]
totintr = totintr[genelist,]

rpkmRep=list()
rpkmRep$labexon=labexon[1:1000,]
rpkmRep$labintr=labintr[1:1000,]
rpkmRep$totexon=totexon[1:1000,]
rpkmRep$totintr=totintr[1:1000,]

save(rpkmRep, file = "rpkmTRUE.RData")
#####################################################
library(INSPEcT)
.find_tt_par <- function(tpts){
  cvLogTpts <- function(a , tpts) {
    newtime <- log2(tpts + a )
    sd(diff(newtime)) / mean(diff(newtime))}
  optimize(f=cvLogTpts, interval=c(0,5), tpts=tpts )$minimum
}
.makeEmptyModel <- function(tpts) {
  model <- matrix(NA, nrow=length(tpts), 5)
  colnames(model) <- c('alpha','beta','gamma','preMRNA','total')
  as.data.frame(model)
}
.time_transf <- function(t, log_shift) {
  newtime <- log2(t+log_shift)
  return(newtime)
}
.rxnrate <- function(t,c,parms){

  # rate constant passed through a list called parms
  alpha <- parms$alpha
  beta  <- parms$beta
  gamma <- parms$gamma

  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1] <- alpha(t) - gamma(t) * c["p"]
  r[2] <- alpha(t) - beta(t) * (c["t"] - c["p"] )

  # c is the concentration of species

  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))

}
.makeModel <- function(tpts, hyp, log_shift){
  params <- list()
  params$alpha <- function(x)
    hyp$alpha$fun$value(.time_transf(x, log_shift), hyp$alpha$par)
  params$beta  <- function(x)
    hyp$beta$fun$value(.time_transf(x, log_shift), hyp$beta$par)
  params$gamma <- function(x)
    hyp$gamma$fun$value(.time_transf(x, log_shift), hyp$gamma$par)
  cinit <- c(params$alpha(tpts[1]) / params$gamma(tpts[1]),
             params$alpha(tpts[1]) / params$beta(tpts[1]) +
               params$alpha(tpts[1]) / params$gamma(tpts[1]))
  names(cinit) <- c('p', 't')
  model <- as.data.frame(
    deSolve::ode(y=cinit, times=tpts, func=.rxnrate, parms=params))
  model$alpha <- params$alpha(tpts)
  model$beta  <- params$beta(tpts)
  model$gamma <- params$gamma(tpts)
  colnames(model)[2:3] <- c('preMRNA','total')
  return(model)
}
makeSimDataset <- function(object, tpts, nRep) {
  ## create the clean concentrations and rates for each gene
  ratesSpecs <- object@ratesSpecs
  nGenes <- length(ratesSpecs)
  log_shift <- .find_tt_par(tpts)
  cleanRates <- lapply(1:nGenes, function(i) {
    tryCatch(
      .makeModel(tpts, ratesSpecs[[i]][[1]], log_shift)
      , error=function(e)
        .makeEmptyModel(tpts)
    )
  })
  ## store total, preMRNA and alpha
  totalSim <- t(sapply(cleanRates, function(x) x$total))
  preMRNASim <- t(sapply(cleanRates, function(x) x$preMRNA))
  alphaSim <- t(sapply(cleanRates, function(x) x$alpha))
  ## get noise variance form the object
  totalSim_noisevar <- object@params$sim$noiseVar$total
  preMRNASim_noisevar <- object@params$sim$noiseVar$pre
  alphaSim_noisevar <- object@params$sim$noiseVar$alpha
  ## simulate the noise
  addNoise <- function(signal, noiseVar) {
    nConditions <- ncol(signal)
    noise <- t(sapply(sqrt(noiseVar),
                      function(sd) rnorm(nConditions,mean=0,sd=sd)))
    out <- signal + noise
    out[out < 0] = 0
    return(out)
  }
  totalSimReplicates <- do.call('cbind',
                                lapply(1:nRep, function(i) addNoise(totalSim,totalSim_noisevar)))
  rownames(totalSimReplicates) <- 1:nGenes
  preMRNASimReplicates <- do.call('cbind',
                                  lapply(1:nRep, function(i) addNoise(preMRNASim,preMRNASim_noisevar)))
  rownames(preMRNASimReplicates) <- 1:nGenes
  alphaSimReplicates <- do.call('cbind',
                                lapply(1:nRep, function(i) addNoise(alphaSim,alphaSim_noisevar)))
  rownames(alphaSimReplicates) <- 1:nGenes
  tpts <- rep(tpts, nRep)
  res = list(tt = totalSimReplicates, tp = preMRNASimReplicates, tl=alphaSimReplicates)
  return(res)
}
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
sim_inspect <- newINSPEcT(TimeGrid, tL,
                        labexon[1:1000,], totexon[1:1000,],
                        labintr[1:1000,], totintr[1:1000,], BPPARAM=SerialParam())
simRates_inspect <- makeSimModel(sim_inspect, 1000, newTpts=NULL)
#TimeGrid  =  c(0, 15, 30, 45, 60)
simData <- makeSimDataset(simRates_inspect, TimeGrid, 2)
simData1 = list(tt=simData$tt[,1:13], tp=simData$tp[,1:13], tl=simData$tl[,1:13])
simData2 = list(tt=simData$tt[,14:26], tp=simData$tp[,14:26], tl=simData$tl[,14:26])

#####################################################
gg=100
testinspect1 <- newINSPEcT(TimeGrid, tL,
                          rpkms_4su_exons=simData1$tl[1:gg,], rpkms_total_exons=simData1$tt[1:gg,],
                          rpkms_total_introns=simData1$tp[1:gg,], BPPARAM=SerialParam())

testmodel1 = estimateParams(simData1$tl[1:gg,],
                            simData1$tt[1:gg,],
                            simData1$tp[1:gg,], TimeGrid, tL, loopnumber=10)

testinspect2 <- newINSPEcT(TimeGrid, tL,
                           rpkms_4su_exons=simData2$tl[1:gg,], rpkms_total_exons=simData2$tt[1:gg,],
                           rpkms_total_introns=simData2$tp[1:gg,], BPPARAM=SerialParam())

testmodel2 = estimateParams(simData2$tl[1:gg,],
                            simData2$tt[1:gg,],
                            simData2$tp[1:gg,], TimeGrid, tL, loopnumber=10)

plotRates(testmodel1,3)

a_in_1 = as.matrix(ratesFirstGuess(testinspect1, 'synthesis'))[testmodel1@genenames,]
b_in_1 = as.matrix(ratesFirstGuess(testinspect1, 'degradation'))[testmodel1@genenames,]

a_in_2 = as.matrix(ratesFirstGuess(testinspect2, 'synthesis'))[testmodel1@genenames,]#[as.numeric(rownames(a_in_1)),]
b_in_2 = as.matrix(ratesFirstGuess(testinspect2, 'degradation'))[testmodel1@genenames,]#[as.numeric(rownames(a_in_1)),]


a1 = getRates(testmodel1, 'transcription', timevector = TimeGrid)
b1 = getRates(testmodel1, 'degradation', timevector = TimeGrid)/as.matrix(simData1$tt[testmodel1@genenames,]-simData1$tp[testmodel1@genenames,])
c1 = getRates(testmodel1, 'processing', timevector = TimeGrid)/as.matrix(simData1$tp[testmodel1@genenames,])

a2 = getRates(testmodel2, 'transcription', timevector = TimeGrid)[testmodel1@genenames,]
b2 = getRates(testmodel2, 'degradation', timevector = TimeGrid)[testmodel1@genenames,]/as.matrix(simData1$tt[testmodel1@genenames,]-simData1$tp[testmodel1@genenames,])
c2 = getRates(testmodel2, 'processing', timevector = TimeGrid)[testmodel1@genenames,]/as.matrix(simData1$tp[testmodel1@genenames,])

par(mfrow = c(2, 3))
scatterdata = data.frame(pluse_a= as.vector(a1),pluse_b= as.vector(b1),pluse_c= as.vector(c1),
                         same_a = as.vector(a2),same_b = as.vector(b2),same_c = as.vector(c2))
scatterdata = log(filter.outers(scatterdata))
heatscatter((as.vector(scatterdata[,1])), (as.vector(as.matrix(scatterdata[,4]))),
            cor=TRUE,method='pearson',main='Transcription-rates',
            ylab='log-trueRates', xlab='log-modelRates')
abline(0,1)
heatscatter((as.vector(scatterdata[,2])), (as.vector(as.matrix(scatterdata[,5]))),
            cor=TRUE,method='pearson',main='degradation',
            ylab='log-trueRates', xlab='log-modelRates')
abline(0,1)
heatscatter((as.vector(scatterdata[,3])), (as.vector(as.matrix(scatterdata[,6]))),
            cor=TRUE,method='pearson',main='processing',
            ylab='log-trueRates', xlab='log-modelRates')
abline(0,1)
scatterdata = data.frame(solver_a= as.vector(a_in_1),solver_b=as.vector(b_in_1),
                         same_a = as.vector(a_in_2),same_b = as.vector(b_in_2))
scatterdata = log(filter.outers(scatterdata))
heatscatter((as.vector(scatterdata[,1])), (as.vector(as.matrix(scatterdata[,3]))),
            cor=TRUE,method='pearson',main='Transcription-rates',
            ylab='log-trueRates', xlab='log-modelRates')
abline(0,1)
heatscatter((as.vector(scatterdata[,2])), (as.vector(as.matrix(scatterdata[,4]))),
            cor=TRUE,method='pearson',main='degradation',
            ylab='log-trueRates', xlab='log-modelRates')
abline(0,1)
heatscatter((as.vector(scatterdata[,2])), (as.vector(as.matrix(scatterdata[,4]))),
            cor=TRUE,method='pearson',main='degradation',
            ylab='log-trueRates', xlab='log-modelRates')
abline(0,1)

######生成仿真#######################################
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
gnum = 1000
#计算速率
resmodel20 = estimateParams(labexon[1:gnum,],
                            totexon[1:gnum,],
                            totintr[1:gnum,], TimeGrid, tL, loopnumber=50)
#normFactors(labexon[1:gnum,], totexon[1:gnum,],TimeGrid, tL,batch = 2, samumber=5, iternumber=5)

save(resmodel20, file = "ratesTRUE.RData")
################################
genename = resmodel20@genenames
a20 = getRates(resmodel20, 'transcription', timevector = TimeGrid)
b20 = getRates(resmodel20, 'degradation', timevector = TimeGrid)/as.matrix(totexon[genename,]-totintr[genename,])
c20 = getRates(resmodel20, 'processing', timevector = TimeGrid)/as.matrix(totintr[genename,])
###############################
### 过滤
filter.outers <- function(data, fold=1.5){
  filter_na = na.omit(data)
  fenwei = quantile(filter_na)
  filter_var <- filter_na[(filter_na < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) &
                             filter_na > (fenwei[2] - (fenwei[4]-fenwei[2])*fold))]
  return(filter_var)
}
##随机抽样
# a
sim_gene = 1000
sim_time = TimeGrid#c(1,3,5,7,9,11)
sim_tl = tL# 1
mean_a = mean(filter.outers(as.vector(a20)))
var_a = var(filter.outers(as.vector(a20)))
sim_a = c()
for(i in 1:sim_gene){
  sim_a = rbind(sim_a, abs(rnorm(length(sim_time), mean_a, var_a)))
}
# b
k = diag(cor(a20, b20))
sim_b = c()
for(i in 1:sim_gene){
  tmp_b = c()
  for(it in 1:length(sim_time)){
    tmp_b = c(tmp_b, abs(rnorm(1, mean_a*k[it], (sqrt(var_a)*k[it])^2)))
  }
  sim_b = rbind(sim_b, tmp_b)
}
# c
kc = diag(cor(a20, c20))
sim_c = c()
for(i in 1:sim_gene){
  tmp_c = c()
  for(it in 1:length(sim_time)){
    tmp_c = c(tmp_c, abs(rnorm(1, mean_a*kc[it], (sqrt(var_a)*kc[it])^2)))
  }
  sim_c = rbind(sim_c, tmp_c)
}
#龙哥库塔计算表达至
rownames(sim_a) = 1:sim_gene
rownames(sim_b) = 1:sim_gene
rownames(sim_c) = 1:sim_gene
ratesSim = list(a=sim_a, b=sim_b, c=sim_c)
save(ratesSim, file = "ratesSim.RData")

at = sim_a
bt = sim_b
ct = sim_c
# 计算T，P，TL初始值
mean_p0 = mean(filter.outers(as.vector(totintr[,1])))
var_p0 = var(filter.outers(as.vector(totintr[,1])))
sim_P0 = abs(rnorm(sim_gene, mean_p0, var_p0))

mean_t0 = mean(filter.outers(as.vector(totexon[,1])))
var_t0 = var(filter.outers(as.vector(totexon[,1])))
sim_T0 = abs(rnorm(sim_gene, mean_t0, var_t0))

sim_PT = c()
sim_TT = c()
sim_TL = at*sim_tl
t = sim_time
RungFunction2 <- function(t, x, parms){
  times = parms[,4]
  tmp = abs(t-times)
  pos = which(tmp==min(tmp))
  if(length(pos)!=1){
    pos = pos[1]
  }
  a = parms[pos,1]
  b = parms[pos,2]
  c = parms[pos,3]
  P = x[1]
  C = x[2]
  dP  <- a - c * P
  dC  <- a - b * (C-P)
  res = c(dP,dC)
  return(list(res))
}
for(i in 1:dim(at)[1]){
  parms  <- cbind(spline(sim_time,at[i,], length(t))$y,
                  spline(sim_time,bt[i,], length(t))$y,
                  spline(sim_time,ct[i,], length(t))$y,
                  t)
  P_T = deSolve::rk(c(P=sim_P0[i], C=sim_T0[i]), t, RungFunction2, parms)
  sim_PT = rbind(sim_PT, P_T[,2])
  sim_TT = rbind(sim_TT, P_T[,3])
}
sim_PL = at / ct * (1 - exp(-ct * sim_tl))

rownames(sim_PT) = 1:sim_gene
rownames(sim_TT) = 1:sim_gene
rownames(sim_TL) = 1:sim_gene
rownames(sim_PL) = 1:sim_gene
colnames(sim_PT) = t
colnames(sim_TT) = t
colnames(sim_TL) = colnames(at)
colnames(sim_PL) = colnames(at)

sim_TL_sample = sim_TL[1:gnum,]
sim_PT_sample = sim_PT[1:gnum,]
sim_TT_sample = sim_TT[1:gnum,]
sim_PL_sample = sim_PL[1:gnum,]

rpkmSim = list()
rpkmSim$labexon = sim_TL_sample
rpkmSim$totintr = sim_PT_sample
rpkmSim$totexon = sim_TT_sample
rpkmSim$labintr = sim_PL_sample
save(rpkmSim, file = "rpkmSim.RData")





