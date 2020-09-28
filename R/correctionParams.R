#' Correct the parameters of the pulse model
#' @description It is used to correct pulse model parameters.
#' During the estimation of parameters in \code{\link{estimateParams}},
#' some NA will occur, which may be caused by random initial values.
#' Therefore, some NA values can be estimated from the new ones.
#' @param object a 'pulseTDmodel' that has been calculated with \code{\link{estimateParams}}
#' @return a 'pulseTDmodel' that has been modified
#' @export
#' @importFrom Biobase pData
#' @examples
#' data('pulseRates', package='pulseTD')
#' pulseRates_correct = correctionParams(pulseRates)
#' pulseRates_correct@fitfailure

correctionParams<-function(object){
  pa = pData(object@ratesPar.transcription)
  pb = pData(object@ratesPar.degradation)
  pc = pData(object@ratesPar.processing)
  labexon = pData(object@filterExpression.foursu_exons)
  totexon = pData(object@filterExpression.total_exons)
  totintr = pData(object@filterExpression.total_introns)
  nanum = c()
  for(i in 1:dim(pa)[1]){
    if('NA' %in% pa[i,] || 'NA' %in% pb[i,] || 'NA' %in% pc[i,]){
      nanum=c(nanum, i)
    }
  }
  filter = c()
  if(length(nanum)!=0){
    message('The number of genes that are incorrectly fitted is ',length(nanum))
    message('The parameters are being corrected...')
    #for(j in nanum){
    #  test = estimateParams(labexon[j,],
    #                        totexon[j,],
    #                        totintr[j,], object@t_time, object@tL, loopnumber=50, message=FALSE)
    #  pa[j,]=pData(test@ratesPar.transcription)
    #  pb[j,]=pData(test@ratesPar.degradation)
    #  pc[j,]=pData(test@ratesPar.processing)
    #  #message(pData(test@ratesPar.transcription), '\n', pData(test@ratesPar.degradation), '\n', pData(test@ratesPar.processing))
    #  if('NA' %in% pa[j,] || 'NA' %in% pb[j,] || 'NA' %in% pc[j,]){
    #    filter=c(filter, j)
    #  }
    #}
  }

  pulsemodel = new('pulseTDmodel')
  if(length(filter)!=0){
    message('The number of genes failed to fit is ', length(filter))
    totintr = totintr[-filter, ]
    totexon = totexon[-filter, ]
    labexon = labexon[-filter, ]
    pa = pa[-filter, ]
    pc = pc[-filter, ]
    pb = pb[-filter, ]
    pulsemodel@fitfailure = filter
  }else{
    message('There are no parameters that need to be corrected')
  }
  pulsemodel@filterExpression.total_introns = new('AnnotatedDataFrame',data=as.data.frame(totintr))
  pulsemodel@filterExpression.total_exons = new('AnnotatedDataFrame',data=as.data.frame(totexon))
  pulsemodel@filterExpression.foursu_exons = new('AnnotatedDataFrame',data=as.data.frame(labexon))

  pulsemodel@ratesPar.transcription =new('AnnotatedDataFrame', data=as.data.frame(pa))
  pulsemodel@ratesPar.processing = new('AnnotatedDataFrame', data=as.data.frame(pc))
  pulsemodel@ratesPar.degradation = new('AnnotatedDataFrame', data=as.data.frame(pb))

  pulsemodel@w1 = object@w1
  pulsemodel@w2 = object@w2
  pulsemodel@genenames = object@genenames
  pulsemodel@t_time = object@t_time
  pulsemodel@tL = object@tL
  return(pulsemodel)
}
