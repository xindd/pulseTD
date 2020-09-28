RatesErrorMat <- function(predictRates,sim_TL_sample, sim_TT_sample, sim_PT_sample,
                          same_a,same_b,same_c,
                          filtergene, TimeGrid, tL, idx=1:5, idx2 = 6:13){
  trainTL = sim_TL_sample[filtergene,idx]
  trainTT = sim_TT_sample[filtergene,idx]
  trainPT = sim_PT_sample[filtergene,idx]
  trainTime = TimeGrid[idx]
  message('计算pluse')
  # pluse
  # predictRates = estimateParams(trainTL, trainTT, trainPT, trainTime, tL, loopnumber=60)
  pluse_a = getRates(predictRates, 'transcription', timevector = TimeGrid)[filtergene,]
  pluse_b = getRates(predictRates, 'degradation', timevector = TimeGrid)[filtergene,]/as.matrix(sim_TT_sample[filtergene,]-sim_PT_sample[filtergene,])
  pluse_c = getRates(predictRates, 'transcription', timevector = TimeGrid)[filtergene,]/as.matrix(sim_PT_sample[filtergene,])
  #poly1
  message('计算poly2')
  reg2 = regress(trainTL, trainTT, trainPT, trainTime, tL, n=2)
  message('计算poly3')
  reg3 = regress(trainTL, trainTT, trainPT, trainTime, tL, n=3)
  message('计算poly4')
  reg4 = regress(trainTL, trainTT, trainPT, trainTime, tL, n=4)
  #solver
  message('计算solver')
  mycerIds <- newINSPEcT(trainTime, tL, trainTL, trainTT,
                         labintr[rownames(trainTL),idx], trainPT, BPPARAM=SerialParam())
  ##
  sat = na.omit(as.matrix(ratesFirstGuess(mycerIds, 'synthesis')))
  sbt = na.omit(as.matrix(ratesFirstGuess(mycerIds, 'degradation')))
  sct = na.omit(as.matrix(ratesFirstGuess(mycerIds, 'processing')))
  myg = intersect(rownames(sat), rownames(sbt))
  myg = intersect(myg, rownames(sct))
  message('NA值数量',length(filtergene)-length(myg))
  dira=c()
  dirb=c()
  dirc=c()
  for(g in myg){
    fitexps = c()
    fitexps2 = c()
    fitexps3 = c()
    for(i in 1:50){
      fitexp = fitRates(sat[g,], trainTime)
      fitexps = rbind(fitexps,fitexp)

      fitexp2 = fitRates(sct[g,], trainTime)
      fitexps2 = rbind(fitexps2,fitexp2)

      fitexp3 = fitRates(sbt[g,], trainTime)
      fitexps3 = rbind(fitexps3,fitexp3)
    }

    bestres_a = fitexps[which(fitexps[,7] == min(fitexps[,7])), 1:6]
    bestres_b = fitexps3[which(fitexps3[,7] == min(fitexps3[,7])), 1:6]
    bestres_c = fitexps2[which(fitexps2[,7] == min(fitexps2[,7])), 1:6]
    tg = TimeGrid
    fit_a = pulseModel(bestres_a,tg)
    fit_c = pulseModel(bestres_c,tg)
    fit_b = pulseModel(bestres_b,tg)
    dira=rbind(dira, fit_a)
    dirb=rbind(dirb, fit_b)
    dirc=rbind(dirc, fit_c)
  }
  tmp = matrix(NA, nrow = length(filtergene), ncol = length(TimeGrid))
  rownames(tmp) = filtergene
  tmp[myg,]=dira
  dira = tmp

  tmp = matrix(NA, nrow = length(filtergene), ncol = length(TimeGrid))
  rownames(tmp) = filtergene
  tmp[myg,]=dirb
  dirb = tmp

  tmp = matrix(NA, nrow = length(filtergene), ncol = length(TimeGrid))
  rownames(tmp) = filtergene
  tmp[myg,]=dirc
  dirc = tmp
  ###############################
  message('计算error')
  errorDF = data.frame()
  for(i in 1:length(filtergene)){
    gn = filtergene[i]
    pluse_a_er = pluse_a[gn,] #
    pluse_b_er = pluse_b[gn,] #
    pluse_c_er = pluse_c[gn,] #

    pre_a2 = regressionModel(reg2$para[i,],TimeGrid)
    pre_b2 = regressionModel(reg2$parb[i,],TimeGrid)
    pre_c2 = regressionModel(reg2$parc[i,],TimeGrid)

    pre_a3 = regressionModel(reg3$para[i,],TimeGrid)
    pre_b3 = regressionModel(reg3$parb[i,],TimeGrid)
    pre_c3 = regressionModel(reg3$parc[i,],TimeGrid)

    pre_a4 = regressionModel(reg4$para[i,],TimeGrid)
    pre_b4 = regressionModel(reg4$parb[i,],TimeGrid)
    pre_c4 = regressionModel(reg4$parc[i,],TimeGrid)

    direct_a = dira[gn,]
    direct_b = dirb[gn,]
    direct_c = dirc[gn,]

    errorDF = rbind(errorDF, c(sum(( pluse_a_er[idx2]  -same_a[gn,idx2])^2),
                               sum((     pre_a2[idx2]  -same_a[gn,idx2])^2),
                               sum((     pre_a3[idx2]  -same_a[gn,idx2])^2),
                               sum((     pre_a4[idx2]  -same_a[gn,idx2])^2),
                               sum((   direct_a[idx2]  -same_a[gn,idx2])^2),
                               sum(( pluse_b_er[idx2]  -same_b[gn,idx2])^2),
                               sum((     pre_b2[idx2]  -same_b[gn,idx2])^2),
                               sum((     pre_b3[idx2]  -same_b[gn,idx2])^2),
                               sum((     pre_b4[idx2]  -same_b[gn,idx2])^2),
                               sum((   direct_b[idx2]  -same_b[gn,idx2])^2),
                               sum(( pluse_c_er[idx2]  -same_c[gn,idx2])^2),
                               sum((     pre_c2[idx2]  -same_c[gn,idx2])^2),
                               sum((     pre_c3[idx2]  -same_c[gn,idx2])^2),
                               sum((     pre_c4[idx2]  -same_c[gn,idx2])^2),
                               sum((   direct_c[idx2]  -same_c[gn,idx2])^2)
    ))
  }
  colnames(errorDF)=c('pluse_a', 'poly2a','poly3a','poly4a', 'direct_a',
                      'pluse_b', 'poly2b','poly3b','poly4b', 'direct_b',
                      'pluse_c', 'poly2c','poly3c','poly4c', 'direct_c')
  return(errorDF)
}
rgFunc<-function(a,b,c, initValue, TimeGrid, tL){
  # 4 计算T，P，TL初始值
  sim_P0 = initValue[,1]#PT[genename,1]#at[,1]/ct[,1]
  sim_T0 = initValue[,2]#TT[genename,1]#at[,1]/bt[,1]+sim_P0
  # 5 龙格库塔计算T，P，TL
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
  sim_PT = c()
  sim_TT = c()
  sim_TL = a*tL
  for(i in 1:dim(a)[1]){
    if(is.na(a[i,1])||is.na(b[i,1])||is.na(c[i,1])){
      sim_PT = rbind(sim_PT, rep(NA,length(TimeGrid)))
      sim_TT = rbind(sim_TT, rep(NA,length(TimeGrid)))

    }else{
      parms  <- cbind(a[i,],
                      b[i,],
                      c[i,],
                      TimeGrid)
      P_T = deSolve::rk(c(P=sim_P0[i], C=sim_T0[i]), TimeGrid, RungFunction2, parms)
      sim_PT = rbind(sim_PT, P_T[,2])
      sim_TT = rbind(sim_TT, P_T[,3])
    }

  }

  rownames(sim_PT) = rownames(a)
  rownames(sim_TT) = rownames(a)
  rownames(sim_TL) = rownames(a)
  return(list(sim_TL=sim_TL, sim_PT=sim_PT, sim_TT=sim_TT))

}
fitRates<-function(rates, TimeGrid){
  errorFun<-function(pars,abc){
    res = pulseModel(pars,abc$TimeGrid)
    error = sum((res - abc$y)^2)
    return(error)
  }
  ###################
  lower = c(0,0.01,0,0,0,0)
  upper = rep(Inf,6)
  info_pt = list(TimeGrid=TimeGrid, y=rates)
  ini = abs(rnorm(6,0,1))
  res_pt = nlminb(start=ini,
                  objective=function(params)errorFun(params,info_pt),
                  gradient = NULL,
                  lower = lower, upper = upper)
  return(c(res_pt$par,res_pt$objective))
}
fitExpression <- function(PT, TT, TimeGrid){
  pulseModel<-function(xita,x){
    h0 = xita[1]
    h1 = xita[2]
    h2 = xita[3]
    t1 = xita[4]
    t2 = xita[5] + t1
    beta=xita[6]
    ( h0 + (h1-h0)/(1+exp(-beta*(x-t1))) ) * ( h2 + (h1-h2)/(1+exp(beta*(x-t2))) ) / h1
  }

  errorFun<-function(pars,abc){
    res = pulseModel(pars,abc$TimeGrid)
    error = sum((res - abc$y)^2)
    return(error)
  }

  ###################
  lower = c(0,0.01,0,0,0,0)
  upper = rep(Inf,6)
  info_pt = list(TimeGrid=TimeGrid, y=PT)
  ini = abs(rnorm(6,0,1))
  res_pt = nlminb(start=ini,
                  objective=function(params)errorFun(params,info_pt),
                  gradient = NULL,
                  lower = lower, upper = upper)

  info_tt = list(TimeGrid=TimeGrid, y=TT)
  res_tt = nlminb(start=ini,
                  objective=function(params)errorFun(params,info_tt),
                  gradient = NULL,
                  lower = lower, upper = upper)

  return(c(res_pt$par,res_pt$objective, res_tt$par, res_tt$objective))
}
errorMat <- function(predictRates,sim_TL_sample, sim_TT_sample, sim_PT_sample, filtergene, TimeGrid, tL, idx=1:5, idx2 = 6:13){
  trainTL = sim_TL_sample[filtergene,idx]
  trainTT = sim_TT_sample[filtergene,idx]
  trainPT = sim_PT_sample[filtergene,idx]
  trainTime = TimeGrid[idx]
  message('计算pluse')
  #pluse
  # predictRates = estimateParams(trainTL, trainTT, trainPT, trainTime, tL, loopnumber=60)
  preExp = predictExpression(predictRates, TimeGrid)
  #poly1
  message('计算poly2')
  reg2 = regress(trainTL, trainTT, trainPT, trainTime, tL, n=2)
  message('计算poly3')
  reg3 = regress(trainTL, trainTT, trainPT, trainTime, tL, n=3)
  message('计算poly4')
  reg4 = regress(trainTL, trainTT, trainPT, trainTime, tL, n=4)
  #solver
  message('计算solver')
  mycerIds <- newINSPEcT(trainTime, tL, trainTL, trainTT,
                         labintr[rownames(trainTL),idx], trainPT, BPPARAM=SerialParam())
  ##
  sat = na.omit(as.matrix(ratesFirstGuess(mycerIds, 'synthesis')))
  sbt = na.omit(as.matrix(ratesFirstGuess(mycerIds, 'degradation')))
  sct = na.omit(as.matrix(ratesFirstGuess(mycerIds, 'processing')))
  myg = intersect(rownames(sat), rownames(sbt))
  myg = intersect(myg, rownames(sct))
  message('NA值数量',length(filtergene)-length(myg))
  dira=c()
  dirb=c()
  dirc=c()
  for(g in myg){
    fitexps = c()
    fitexps2 = c()
    fitexps3 = c()
    for(i in 1:50){
      fitexp = fitRates(sat[g,], trainTime)
      fitexps = rbind(fitexps,fitexp)

      fitexp2 = fitRates(sct[g,], trainTime)
      fitexps2 = rbind(fitexps2,fitexp2)

      fitexp3 = fitRates(sbt[g,], trainTime)
      fitexps3 = rbind(fitexps3,fitexp3)
    }

    bestres_a = fitexps[which(fitexps[,7] == min(fitexps[,7])), 1:6]
    bestres_b = fitexps3[which(fitexps3[,7] == min(fitexps3[,7])), 1:6]
    bestres_c = fitexps2[which(fitexps2[,7] == min(fitexps2[,7])), 1:6]
    tg = TimeGrid
    fit_a = pulseModel(bestres_a,tg)
    fit_c = pulseModel(bestres_c,tg)
    fit_b = pulseModel(bestres_b,tg)
    dira=rbind(dira, fit_a)
    dirb=rbind(dirb, fit_b)
    dirc=rbind(dirc, fit_c)
  }
  tmp = matrix(NA, nrow = length(filtergene), ncol = length(TimeGrid))
  rownames(tmp) = filtergene
  tmp[myg,]=dira
  dira = tmp

  tmp = matrix(NA, nrow = length(filtergene), ncol = length(TimeGrid))
  rownames(tmp) = filtergene
  tmp[myg,]=dirb
  dirb = tmp

  tmp = matrix(NA, nrow = length(filtergene), ncol = length(TimeGrid))
  rownames(tmp) = filtergene
  tmp[myg,]=dirc
  dirc = tmp
  message('direct Pluse')
  directPluse = rgFunc(dira,dirb,dirc, cbind(trainPT[,1], trainTT[,1]), tg, tL)
  #exp
  message('计算exp')
  exp_PP=c()
  exp_TT=c()
  for(g in filtergene){
    expRes = c()
    for(i in 1:50){
      trainExp = fitExpression(trainPT[g,], trainTT[g,], trainTime)
      expRes = rbind(expRes, trainExp)
    }
    exp_PP = rbind(exp_PP, pulseModel(expRes[order(expRes[,7])[1], 1:6], TimeGrid))
    exp_TT = rbind(exp_TT, pulseModel(expRes[order(expRes[,14])[1], 8:13], TimeGrid))
  }
  ###############################
  message('计算error')
  errorDF = data.frame()
  for(i in 1:length(filtergene)){
    gn = filtergene[i]
    preE = data.frame(preExp[[gn]])
    plusePT = sum((preE$PT[idx2]-sim_PT_sample[gn,idx2])^2)
    pluseTT = sum((preE$TT[idx2]-sim_TT_sample[gn,idx2])^2)
    pre_p2 = integrateFunc(reg2$para[i,], TimeGrid)-integrateFunc(reg2$parc[i,], TimeGrid) + sim_PT_sample[gn,1]
    pre_t2 = integrateFunc(reg2$para[i,], TimeGrid)-integrateFunc(reg2$parb[i,], TimeGrid) + sim_TT_sample[gn,1]
    pre_p3 = integrateFunc(reg3$para[i,], TimeGrid)-integrateFunc(reg3$parc[i,], TimeGrid) + sim_PT_sample[gn,1]
    pre_t3 = integrateFunc(reg3$para[i,], TimeGrid)-integrateFunc(reg3$parb[i,], TimeGrid) + sim_TT_sample[gn,1]
    pre_p4 = integrateFunc(reg4$para[i,], TimeGrid)-integrateFunc(reg4$parc[i,], TimeGrid) + sim_PT_sample[gn,1]
    pre_t4 = integrateFunc(reg4$para[i,], TimeGrid)-integrateFunc(reg4$parb[i,], TimeGrid) + sim_TT_sample[gn,1]
    direct_p = sum((directPluse$sim_PT[gn,idx2]-sim_PT_sample[gn,idx2])^2)
    direct_t = sum((directPluse$sim_TT[gn,idx2]-sim_TT_sample[gn,idx2])^2)
    exp_p = sum((exp_PP[i,idx2]-sim_PT_sample[gn,idx2])^2)
    exp_t = sum((exp_TT[i,idx2]-sim_TT_sample[gn,idx2])^2)

    errorDF = rbind(errorDF, c(plusePT,
                               pluseTT,
                               sum((pre_p2[idx2]-sim_PT_sample[gn,idx2])^2),
                               sum((pre_t2[idx2]-sim_TT_sample[gn,idx2])^2),
                               sum((pre_p3[idx2]-sim_PT_sample[gn,idx2])^2),
                               sum((pre_t3[idx2]-sim_TT_sample[gn,idx2])^2),
                               sum((pre_p4[idx2]-sim_PT_sample[gn,idx2])^2),
                               sum((pre_t4[idx2]-sim_TT_sample[gn,idx2])^2),
                               direct_p,direct_t,
                               exp_p, exp_t))
  }
  colnames(errorDF)=c('plusePT','pluseTT','pre_p2','pre_t2',
                      'pre_p3','pre_t3','pre_p4','pre_t4',
                      'direct_p', 'direct_t',
                      'exp_p', 'exp_t')
  return(errorDF)
}
regressionModel <- function(xita, t){
  tmp = as.function(polynom::polynomial(xita))(t)
  #tmp = exp(tmp)
  return(tmp)
}
calLabel <- function(xita, timegrid, tl){
  res=sapply(timegrid,function(xt)(regressionModel(xita,xt-tl) + regressionModel(xita,xt))*tl/2)
  #res=sapply(timegrid,function(xt)(integrate(function(tt)exp(as.function(polynom::polynomial(xita))(tt)), lower=xt-tl, upper=xt)$value))
  return(res)
}
integrateFunc <- function(xita, timegrid){
  res=sapply(timegrid,function(xt)(integrate(function(tt)regressionModel(xita, tt), lower=0, upper=xt)$value))
  return(res)
}
regerrorAFunc = function(params, abc){
  na = abc[[1]]+1
  errorLabel = sum((calLabel(params, abc$TimeGrid,abc$tL)-abc$l)^2)
  error = errorLabel/2
  return(error)
}
regerrorCFunc = function(params, abc){
  error = sum((integrateFunc(as.matrix(params), abc$TimeGrid) - abc$p)^2)
  error = error/2
  return(error)
}
regressFun = function(PT,TT,TL,TimeGrid,tL, na=3,nb=3,nc=3){
  error=c()
  ############### aaaaaaa #####################
  info_a = list(na, TimeGrid=TimeGrid, l=TL,tL=tL)
  ini_a = rnorm(na+1,0,0.1)
  tmpres = optim(par = ini_a,
                 fn = function(params)regerrorAFunc(params,info_a),
                 #gr = function(params)gradAFunc(params,info_a),
                 method = "BFGS")
  pre_a = as.matrix(tmpres$par)
  error = c(error, tmpres$value)
  ############### ccccccc #####################
  PTC = integrateFunc(pre_a, TimeGrid)-(PT - PT[1])
  info_c = list(na=na, nc=nc, parm_a=pre_a,TimeGrid=TimeGrid, p=PTC,tL=tL,lamda=0.1)
  ini_c = rnorm(nc+1,0,0.1)
  integrateFunc(ini_c, TimeGrid)

  tmpres = optim(par = ini_c,
                 fn = function(params)regerrorCFunc(params,info_c),
                 #gr = function(params)gradCFunc(params,info_c),
                 method = "BFGS")
  error = c(error, tmpres$value)
  pre_c = as.matrix(tmpres$par)
  ############### bbbbbbb #####################
  TTB = integrateFunc(pre_a, TimeGrid)-(TT - TT[1])
  info_b = list(na=na,nb=nb, parm_a=pre_a,TimeGrid=TimeGrid, p = TTB,tL=tL)
  ini_b = rnorm(nb+1,0,0.000001)
  tmpres = optim(par = ini_b,
                 fn = function(params)regerrorCFunc(params,info_b),
                 #gr = function(params)gradAFunc(params,info_c),
                 method = "BFGS")
  error = c(error, tmpres$value)
  pre_b = as.matrix(tmpres$par)
  ######################

  par(mfrow = c(2, 3))
  tg=TimeGrid#seq(0,200,1)
  pre_TL = calLabel(pre_a,TimeGrid, tL)
  pre_p = integrateFunc(pre_a, tg)-integrateFunc(pre_c, tg) + PT[1]
  pre_t = integrateFunc(pre_a, tg)-integrateFunc(pre_b, tg) + TT[1]

  predict_a = regressionModel(pre_a, tg)
  predict_b = regressionModel(pre_b, tg)
  predict_c = regressionModel(pre_c, tg)
  #plot(tg, predict_a, type='o')
  #plot(tg, predict_c, type='o')
  #plot(tg, predict_b, type='o')
  #return(list(pa=pre_a, pb=pre_b, pc=pre_c, error=error))
  return(list(a=predict_a, b=predict_b, c=predict_c, para=pre_a, parb=pre_b, parc=pre_c))
}
regress <- function(TL, TT, PT, TimeGrid, tL, n=3){
  reg_a2= c()
  reg_b2= c()
  reg_c2= c()
  reg_para=c()
  reg_prab=c()
  reg_prac=c()
  for(i in rownames(TL)){
    regres = list(a=rep(NA,length(TimeGrid)),b=rep(NA,length(TimeGrid)),c=rep(NA,length(TimeGrid)),
                  para=as.matrix(rep(0,n+1)),parb=as.matrix(rep(0,n+1)),parc=as.matrix(rep(0,n+1)))

    tryCatch({
      regres = regressFun(PT[i,],TT[i,],TL[i,],TimeGrid,tL, na=n,nb=n,nc=n)
    },
    error = function(e){}
    )
    reg_a2= rbind(reg_a2, regres$a)
    reg_b2= rbind(reg_b2, regres$b)
    reg_c2= rbind(reg_c2, regres$c)
    reg_para=rbind(reg_para,t(regres$para))
    reg_prab=rbind(reg_prab,t(regres$parb))
    reg_prac=rbind(reg_prac,t(regres$parc))
  }
  return(list(a=reg_a2,
              b=as.matrix(reg_b2/(TT-PT)),
              c=as.matrix(reg_c2/PT),
              para=reg_para,
              parb=reg_prab,
              parc=reg_prac
  ))
}

dropNaModel<-function(object){
  dropvector = c()
  for(g in 1:length(object@genenames)){
    if(sum(is.na(object@ratesPar.transcription@data[g,]))!=0 ||
       sum(is.na(object@ratesPar.processing@data[g,]))!=0 ||
       sum(is.na(object@ratesPar.degradation@data[g,]))!=0 ||
       sum(is.na(object@filterExpression.total_introns@data[g,]))!=0 ||
       sum(is.na(object@filterExpression.total_exons@data[g,]))!=0 ||
       sum(is.na(object@filterExpression.foursu_exons@data[g,]))!=0)
      dropvector = c(dropvector, g)
  }
  fg = object@genenames[-dropvector]

  pulsemodel = new('pulseTDmodel')
  pulsemodel@filterExpression.total_introns = new('AnnotatedDataFrame',data=object@filterExpression.total_introns@data[fg,])
  pulsemodel@filterExpression.total_exons = new('AnnotatedDataFrame',data=object@filterExpression.total_exons@data[fg,])
  pulsemodel@filterExpression.foursu_exons = new('AnnotatedDataFrame',data=object@filterExpression.foursu_exons@data[fg,])

  pulsemodel@ratesPar.transcription =new('AnnotatedDataFrame', data=object@ratesPar.transcription@data[fg,])
  pulsemodel@ratesPar.processing = new('AnnotatedDataFrame', data=object@ratesPar.processing@data[fg,])
  pulsemodel@ratesPar.degradation = new('AnnotatedDataFrame', data=object@ratesPar.degradation@data[fg,])

  pulsemodel@w1 = object@w1
  pulsemodel@w2 = object@w2
  pulsemodel@chisq = object@chisq[-dropvector]

  pulsemodel@genenames = fg
  pulsemodel@t_time = object@t_time
  pulsemodel@tL = object@tL
  return(pulsemodel)
}


library(pulseTD)
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10

dir = 'E:\\pluseTD\\pulseTD-2'
##
load(file.path(dir,'rpkmTRUE.RData'))
labexon = rpkmRep$labexon
labintr = rpkmRep$labintr
totexon = rpkmRep$totexon
totintr = rpkmRep$totintr

load(file.path(dir,'rpkmSim.RData'))
sim_TL_sample = rpkmSim$labexon
sim_PT_sample = rpkmSim$totintr
sim_TT_sample = rpkmSim$totexon
########
#poly
reg2 = regress(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, n=2)
reg3 = regress(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, n=3)
reg4 = regress(sim_TL_sample, sim_TT_sample, sim_PT_sample, TimeGrid, tL, n=4)
save(reg2, file = "reg2.RData")
save(reg3, file = "reg3.RData")
save(reg4, file = "reg4.RData")
#solver
library(INSPEcT)
mycerIds1 <- newINSPEcT(TimeGrid, tL, labexon[rownames(sim_TL_sample),], totexon[rownames(sim_TL_sample),],
                        labintr[rownames(sim_TL_sample),], totintr[rownames(sim_TL_sample),], BPPARAM=SerialParam())

myat = as.matrix(ratesFirstGuess(mycerIds1, 'synthesis'))
mybt = as.matrix(ratesFirstGuess(mycerIds1, 'degradation'))
myct = as.matrix(ratesFirstGuess(mycerIds1, 'processing'))
solver = list(a = myat, b=mybt, c=myct)
save(solver, file = "solver.RData")

######################
load(file.path(dir,'pluse5.RData'))
load(file.path(dir,'pluse7.RData'))
load(file.path(dir,'pluse9.RData'))
load(file.path(dir,'pluse11.RData'))
##
pluse5 = dropNaModel(pluse5)
pluse7 = dropNaModel(pluse7)
pluse9 = dropNaModel(pluse9)
pluse11 = dropNaModel(pluse11)

genename = pluse5@genenames
####################
error_5 = errorMat(pluse5, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:5, idx2 = 6:13)
error_7 = errorMat(pluse7, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:7, idx2 = 8:13)
error_9 = errorMat(pluse9, sim_TL_sample, sim_TT_sample, sim_PT_sample, genename, TimeGrid, tL, idx=1:9, idx2 = 10:13)
error_11= errorMat(pluse11, sim_TL_sample, sim_TT_sample, sim_PT_sample,genename, TimeGrid, tL, idx=1:11, idx2 = 12:13)
save(error_5, file = "error_5.RData")
save(error_7, file = "error_7.RData")
save(error_9, file = "error_9.RData")
save(error_11, file = "error_11.RData")
#####rates error###########
load(file.path(dir,'ratesTrue.RData'))
# resmodel20 = correctionParams(resmodel20)
#genename = resmodel20@genenames
same_a = ratesTrue$a[genename,]
same_b = ratesTrue$b[genename,]
same_c = ratesTrue$c[genename,]
##
genename=pluse5@genenames
rates_error_5 = RatesErrorMat(pluse5, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                              same_a,same_b,same_c,
                              genename, TimeGrid, tL, idx=1:5, idx2 = 6:13)
save(rates_error_5, file = "rates_error_5.RData")
##
genename=pluse7@genenames
rates_error_7 = RatesErrorMat(pluse7, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                              same_a,same_b,same_c,
                              genename, TimeGrid, tL, idx=1:7, idx2 = 8:13)
save(rates_error_7, file = "rates_error_7.RData")
##
genename=pluse9@genenames
rates_error_9 = RatesErrorMat(pluse9, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                              same_a,same_b,same_c,
                              genename, TimeGrid, tL, idx=1:9, idx2 = 10:13)
save(rates_error_9, file = "rates_error_9.RData")
##
genename=pluse11@genenames
rates_error_11=RatesErrorMat(pluse11, sim_TL_sample, sim_TT_sample, sim_PT_sample,
                             same_a,same_b,same_c,
                             genename, TimeGrid, tL, idx=1:11, idx2 = 12:13)
save(rates_error_11, file = "rates_error_11.RData")


