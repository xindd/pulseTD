#' Estimated pulse model parameters!
#'
#' @description It is used to estimate pulse model parameters.
#' Transcription rates, processing rates, and degradation correspond to different pulse model parameters, respectively.
#' The transcription rate has 6 parameters, the processing rate has 5 parameters, and the degradation rate has 4 parameters.
#'
#' @param TimeGrid A vector of points in time at which experimental data is collected, not allowed to be repeated
#' @param tL 4sU labeled time during the experiment
#' @param labexon A matrix containing expression levels of 4su exons.For the calculation of the expression value, see \code{\link{estimateExpression}}
#' @param totexon A matrix containing expression levels of total exons.For the calculation of the expression value, see \code{\link{estimateExpression}}
#' @param totintr A matrix containing expression levels of total introns.For the calculation of the expression value, see \code{\link{estimateExpression}}
#' @param clusterNumber Given the number of cluster cores, the default is the maximum number of available cores
#' @param loopnumber The number of iterations of the gradient descent when solving the parameter. The default is 50.
#' @param message Whether to print the log, the default is TRUE
#' @return A 'pulseTDmodel' containing the expression values of the filtered genes, a list of solved parameters, and some basic parameter information
#' @import parallel
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom graphics par
#' @importFrom stats integrate nlminb qt quantile rnorm na.omit var
#' @export
#' @examples
#' data('rpkmSim', package='pulseTD')
#' rpkm_TL <- rpkmSim$labexon[1:2,]
#' rpkm_PT <- rpkmSim$totintr[1:2,]
#' rpkm_TT <- rpkmSim$totexon[1:2,]
#' TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
#' tL <- 10
#' \donttest{
#' pulseRates<- estimateParams(rpkm_TL,rpkm_TT,rpkm_PT,TimeGrid,tL,clusterNumber=1,loopnumber=10)
#' }
#'
estimateParams<-function(labexon, totexon, totintr, TimeGrid, tL,
                         clusterNumber=NULL,loopnumber=50,message=TRUE, normal='None'){
  ## get Platform
  plat = c(sessionInfo()$running[1],sessionInfo()$platform)
  if(length(grep(pattern='ubutun', x=plat, ignore.case = TRUE))!=0 ||
     length(grep(pattern='linux', x=plat, ignore.case = TRUE))!=0 ||
     length(grep(pattern='red hat', x=plat, ignore.case = TRUE))!=0){
    platform = 'linux'
  }else{
    platform = 'windows'
  }
  if(message){
    message('Operating platform ', platform)
  }
  ############################################
  labexon_filter = filterzeros(labexon)
  totexon_filter = filterzeros(totexon)
  totintr_filter = filterzeros(totintr)

  genenamelist=intersect(rownames(labexon_filter)
                         ,rownames(totexon_filter))
  genenamelist=intersect(genenamelist,rownames(totintr_filter))

  totintr = rbind(totintr[genenamelist,])
  totexon = rbind(totexon[genenamelist,])
  labexon = rbind(labexon[genenamelist,])
  rownames(totintr) = genenamelist
  rownames(totexon) = genenamelist
  rownames(labexon) = genenamelist
  normFunc <- function(data){
    return((data-min(data)) / (max(data)-min(data)))
  }
  if(normal == 'log'){
    totintr = log2(totintr+1)
    totexon = log2(totexon+1)
    labexon = log2(labexon+1)
  }else if(normal == 'min-max-normalization'){
    totintr = apply(totintr, 2, normFunc)
    totexon = apply(totexon, 2, normFunc)
    labexon = apply(labexon, 2, normFunc)

  }else if(normal =='None'){

  }else{
    message('Error: Please enter the correct normalization method.')
  }

  if(message){
    message('Partially zero expressed genes have been filtered.')
    message('The number of genes used to calculate the rate of transcriptional dynamics is ',length(genenamelist))
  }
  #############################################
  genenames = genenamelist
  geneNumbers = length(genenames)
  rates_paramslist=list()

  if(is.null(clusterNumber)){
    cl.cores = detectCores()
  }else{
    cl.cores = clusterNumber
  }
  if(message){
    message('The number of threads used is ', cl.cores)
  }
  cl <- makeCluster(cl.cores)
  clusterExport(cl, c("splitA6Params","splitB4Params", "splitC5Params",
                      "loopnumber",'tL',"TimeGrid"),envir = environment())

  # normalization factor
  if(message){
    message('Estimating global normalization factor...')
  }
  clusterExport(cl,c("labexon", "totexon", "totintr"),envir = environment())
  factors_w <- parLapply(cl,1:(loopnumber*2), function(x){normFactors(labexon, totexon,
                                                            TimeGrid, tL,
                                                            batch = 2, samumber=5, iternumber=5)})
  factors_w <- do.call('rbind',factors_w)
  bestres_w = as.vector(round(apply(factors_w, 2, median), 4))
  w1 = bestres_w[1:length(TimeGrid)]
  w2 = bestres_w[(length(TimeGrid)+1):(2*length(TimeGrid))]
  # w1 = rep(0,length(TimeGrid))
  # w2 = rep(1,length(TimeGrid))
  #########################
  if(message){
    message('Estimating pulseModel parameters')
  }
  if(message){
    pro = 0
    sp = 10
    pre = 0
    p = geneNumbers/sp
    message('Completed 0%')
    t1=Sys.time()
  }
  for(gene in genenames){
    TL = as.matrix(labexon[gene,])
    TT = as.matrix(totexon[gene,])
    PT = as.matrix(totintr[gene,])

    clusterExport(cl,c("PT", "TT", "TL", 'w1', 'w2'),envir = environment())
    ##############
    dta <- parLapply(cl,1:(loopnumber*2), function(x){splitA6Params(TT, TL, tL, TimeGrid, w1, w2)})
    res_A <- do.call('rbind',dta)
    bestres_a = res_A[which(res_A[,7]==min(res_A[,7]))[1], 1:6]
    clusterExport(cl, "bestres_a",envir = environment())
    ##############
    dtc <- parLapply(cl,1:loopnumber, function(x){splitB4Params(PT, bestres_a, TimeGrid)})
    res_C <- do.call('rbind',dtc)
    bestres_c = res_C[which(res_C[,5]==min(res_C[,5]))[1], 1:4]
    ###############
    dtb <- parLapply(cl,1:loopnumber, function(x){splitB4Params(TT, bestres_a, TimeGrid)})
    res_B <- do.call('rbind',dtb)
    bestres_b = res_B[which(res_B[,5]==min(res_B[,5]))[1], 1:4]
    ###############
    rates_paramslist$at=rbind(rates_paramslist$at, bestres_a)
    rates_paramslist$ct=rbind(rates_paramslist$ct, c(bestres_a[1], bestres_c[1],bestres_a[3], bestres_c[2:4]))
    rates_paramslist$bt=rbind(rates_paramslist$bt, c(bestres_a[1], bestres_b[1],bestres_a[3], bestres_b[2:4]))
    if(message){
      pro = pro+1
      if(pre != pro%/%p){
        t2=Sys.time()

        message('Completed ',round(pro%/%p/sp,2)*100,'%',
                ', cost: ', round(as.numeric(t2-t1),2),units((t2-t1)),
                ', eta:', round(as.numeric(t2-t1)/(pro%/%p)*10-as.numeric(t2-t1),2),units((t2-t1)))
        pre = pro%/%p
      }
    }
  }
  stopCluster(cl)


  rownames(rates_paramslist$at) = genenames
  rownames(rates_paramslist$ct) = genenames
  rownames(rates_paramslist$bt) = genenames
  ###################
  pre_a = rates_paramslist$at
  pre_c = rates_paramslist$ct
  pre_b = rates_paramslist$bt

  chisq_res = c()
  for(i in genenames){
    tmp = NaN
    tryCatch({
      pre_TL = w1 * totexon[i,] + w2 * integrateA(pre_a[i,],TimeGrid, tL)
      pre_PT = integrateAC(pre_a[i,],TimeGrid) - integrateAC(pre_c[i,],TimeGrid) + totintr[i,1]
      pre_TT = integrateAC(pre_a[i,],TimeGrid) - integrateAC(pre_b[i,],TimeGrid) + totexon[i,1]
      Ox = c(pre_TL, pre_PT, pre_TT)
      Tx = c(labexon[i,], totintr[i,], totexon[i,])
      options(warn=-1)
      tmp = NaN
      tab = rbind(Ox, Tx)
      tab[tab<0]=0
      tryCatch({
        tmp = chisq.test(tab)$p.value
      },
      error = function(e) {}
      )
    },error = function(e) {}
    )
    chisq_res = c(chisq_res, tmp)
    # R2 = 1-(sum((Ox-Tx)^2)/sum((Tx-mean(Tx)^2)))
  }
  ###################################
  pulsemodel = new('pulseTDmodel')
  pulsemodel@filterExpression.total_introns = new('AnnotatedDataFrame',data=as.data.frame(totintr))
  pulsemodel@filterExpression.total_exons = new('AnnotatedDataFrame',data=as.data.frame(totexon))
  pulsemodel@filterExpression.foursu_exons = new('AnnotatedDataFrame',data=as.data.frame(labexon))

  pulsemodel@ratesPar.transcription =new('AnnotatedDataFrame', data=as.data.frame(rates_paramslist$at))
  pulsemodel@ratesPar.processing = new('AnnotatedDataFrame', data=as.data.frame(rates_paramslist$ct))
  pulsemodel@ratesPar.degradation = new('AnnotatedDataFrame', data=as.data.frame(rates_paramslist$bt))

  pulsemodel@w1 = w1
  pulsemodel@w2 = w2
  pulsemodel@chisq = chisq_res

  pulsemodel@genenames = genenames
  pulsemodel@t_time = TimeGrid
  pulsemodel@tL = tL
  return(pulsemodel)
}




