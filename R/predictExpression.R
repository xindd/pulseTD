#' predictExpression-method
#' @description
#' It is used to predict the expression of all gene at a given time,
#' including the expression of pre-mRNA (Precursor RNA) and the expression of total mRNA (Mature RNA).
#' End time and time interval can be arbitrarily defined
#'
#' @param object a 'pulseTDmodel' that has been calculated with \code{\link{estimateParams}}
#' @param tg A vector of points in time at which experimental data is collected, not allowed to be repeated
#' @return Returns a list containing predicted values for each gene and a 0.95 confidence interval
#' @examples
#' data('pulseRates', package='pulseTD')
#' pulseRates_correct = correctionParams(pulseRates)
#' TimeGrid = seq(0,180,15)
#' \donttest{preExp = predictExpression(pulseRates_correct, tg=TimeGrid)}
#' data('preExp', packages='pulseTD')
#' df = data.frame(preExp[['NM_001002011']])
#' head(df)
#'
setMethod('predictExpression', 'pulseTDmodel', function(object, tg){
  a = as.matrix(pData(object@ratesPar.transcription))
  b = as.matrix(pData(object@ratesPar.degradation))
  c = as.matrix(pData(object@ratesPar.processing))
  TimeGrid = object@t_time
  tL = object@tL
  ptdata = as.matrix(pData(object@filterExpression.total_introns))
  ttdata = as.matrix(pData(object@filterExpression.total_exons))
  pldata = as.matrix(pData(object@filterExpression.foursu_exons))
  n=length(tg)
  res = list()
  genenamelist = rownames(ptdata)
  for(gene in genenamelist){
    pre_a = a[gene,]
    pre_c = c[gene,]
    pre_b = b[gene,]
    PT = ptdata[gene,]
    TT = ttdata[gene,]

    pre_PT = integrateAC(pre_a,tg) - integrateAC(pre_c,tg) + PT[1]
    pre_TT = integrateAC(pre_a,tg) - integrateAC(pre_b,tg) + TT[1]

    up_PT = c()
    down_PT = c()
    up_TT = c()
    down_TT =c()
    for(i in 1:n){
      X = tg[i]
      Yj = pre_PT[i]
      tmp = Yj + qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_PT-Yj)^2) / (n-2)) * sqrt(1/n + (X -mean(tg))^2 / ((n-1)*var(tg))) #sum((tg - mean(tg))^2))
      tmp2 = Yj - qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_PT-Yj)^2) / (n-2)) * sqrt(1/n + (X -mean(tg))^2 / ((n-1)*var(tg))) #sum((tg - mean(tg))^2))
      tmp3 = pre_TT[i] + qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_TT-pre_TT[i])^2) / (n-2)) * sqrt(1/n + (X -mean(tg))^2 / ((n-1)*var(tg))) #sum((tg - mean(tg))^2))
      tmp4 = pre_TT[i] - qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_TT-pre_TT[i])^2) / (n-2)) * sqrt(1/n + (X -mean(tg))^2 / ((n-1)*var(tg))) #sum((tg - mean(tg))^2))

      up_PT = c(up_PT, tmp)
      down_PT = c(down_PT, tmp2)
      up_TT = c(up_TT, tmp3)
      down_TT = c(down_TT, tmp4)
    }
    res[[gene]] = list(PT=pre_PT,TT=pre_TT,
                           upPT=up_PT,downPT=down_PT,
                           upTT=up_TT,downTT=down_TT)
  }

  return(res)

})
