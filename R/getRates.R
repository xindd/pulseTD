#' getRates-method
#' @description
#' It is used to obtain transcription, processing and degradation rates.
#' @param object An object of class 'pulseTDmodel'
#' @param stage A character such as transcription, processing or degradation for rates.
#' @param timevector A vector of times, default is NULL.
#' @param genename A vector of gene names, default is NULL.
#' @return A numeric matrix containing the rates


setMethod('getRates', 'pulseTDmodel', function(object, stage, timevector=NULL,genename=NULL) {
  if(is.null(timevector))timevector=object@t_time
  if(is.null(genename)){
    if(stage=='transcription'){
      a=t(apply(pData(object@ratesPar.transcription), 1, function(x)pulseModel(x, timevector)))
      colnames(a) = timevector
      return(a)
    }else if(stage=='processing'){
      c = t(apply(pData(object@ratesPar.processing), 1, function(x)pulseModel(x, timevector)))
      colnames(c) = timevector
      return(c)
    }else if(stage=='degradation'){
      b = t(apply(pData(object@ratesPar.degradation), 1, function(x)pulseModel(x, timevector)))
      colnames(b) = timevector
      return(b)
    }else message('Please enter the correct stage parameters, you can try to use transcription, processing or degradation')
  }else{
    if(stage=='transcription'){
      a=t(apply(pData(object@ratesPar.transcription)[genename,], 1, function(x)pulseModel(x, timevector)))
      colnames(a) = timevector
      return(a)
    }else if(stage=='processing'){
      c = t(apply(pData(object@ratesPar.processing)[genename,], 1, function(x)pulseModel(x, timevector)))
      colnames(c) = timevector
      return(c)
    }else if(stage=='degradation'){
      b = t(apply(pData(object@ratesPar.degradation)[genename,], 1, function(x)pulseModel(x, timevector)))
      colnames(b) = timevector
      return(b)
    }else message('Please enter the correct stage parameters, you can try to use transcription, processing or degradation')
  }

})
