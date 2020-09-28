#' plotRate-Method
#' @description
#' Draw the image of the fitting result.
#' @param object An object of class 'pulseTDmodel'
#' @param genename The name of the gene you want to draw.
#' @param predict The time interval you want to draw.
#' @return No return value
#'
setMethod('plotRates', 'pulseTDmodel', function(object, genename,predict=FALSE) {
  if(length(genename)!=1){
    stop('please input correct gene name, but only one')
  }else if(!genename %in% object@genenames && !genename  %in% 1:length(object@genenames)){
    stop('Gene name does not existe')
  }else{
    plotSingleGene(parlist=list(pre_a = as.matrix(pData(object@ratesPar.transcription)[genename,]),
                                pre_c = as.matrix(pData(object@ratesPar.processing)[genename,]),
                                pre_b = as.matrix(pData(object@ratesPar.degradation)[genename,])),
                   PT = as.matrix(pData(object@filterExpression.total_introns)[genename,]),
                   TT = as.matrix(pData(object@filterExpression.total_exons)[genename,]),
                   TL = as.matrix(pData(object@filterExpression.foursu_exons)[genename,]),
                   w1 = object@w1,
                   w2 = object@w2,
                   tL = object@tL,
                   TimeGrid = object@t_time,
                   predict=predict
                   )
  }
})

