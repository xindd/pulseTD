#' getParams-method
#' @description
#' It is used to obtain transcription, processing and degradation rates.
#' @param object An object of class 'pulseTDmodel'
#' @param stage A character such as transcription, processing or degradation for rates.
#' @param genename A vector of gene names, default is NULL.
#' @return A numeric matrix containing the rates

setMethod('getParams', 'pulseTDmodel', function(object, stage, genename=NULL) {
  timevector=object@t_time
  if(is.null(genename)){
    if(stage=='transcription'){
      tmp=pData(object@ratesPar.transcription)
      colnames(tmp) = c('h_0', 'h_1', 'h_2', 't_1', 't_2', 'beta')
      return(tmp)
    }else if(stage=='processing'){
      tmp=pData(object@ratesPar.processing)
      colnames(tmp) = c('h_0', 'h_1', 'h_2', 't_1', 't_2', 'beta')
      return(tmp)
    }else if(stage=='degradation'){
      tmp=pData(object@ratesPar.degradation)
      colnames(tmp) = c('h_0', 'h_1', 'h_2', 't_1', 't_2', 'beta')
      return(tmp)
    }else message('Please enter the correct stage parameters, you can try to use transcription, processing or degradation')
  }else{
    if(stage=='transcription'){
      tmp=pData(object@ratesPar.transcription)[genename,]
      colnames(tmp) = c('h_0', 'h_1', 'h_2', 't_1', 't_2', 'beta')
      return(tmp)
    }else if(stage=='processing'){
      tmp = pData(object@ratesPar.processing)[genename,]
      colnames(tmp) = c('h_0', 'h_1', 'h_2', 't_1', 't_2', 'beta')
      return(tmp)
    }else if(stage=='degradation'){
      tmp=pData(object@ratesPar.degradation)[genename,]
      colnames(tmp) = c('h_0', 'h_1', 'h_2', 't_1', 't_2', 'beta')
      return(tmp)
    }else message('Please enter the correct stage parameters, you can try to use transcription, processing or degradation')
  }

})
