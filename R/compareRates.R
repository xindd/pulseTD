#' compareRates-Method
#' @description
#' compare result.
#' @param object1 An object of class 'pulseTDmodel'
#' @param object1 An object of class 'pulseTDmodel'
#' @return A list
#'
setMethod('compareRates', 'pulseTDmodel', function(object1, object2) {
  a1 = as.vector(pData(object1@ratesPar.transcription))
  c1 = as.vector(pData(object1@ratesPar.processing))
  b1 = as.vector(pData(object1@ratesPar.degradation))
  a2 = as.vector(pData(object2@ratesPar.transcription))
  c2 = as.vector(pData(object2@ratesPar.processing))
  b2 = as.vector(pData(object2@ratesPar.degradation))
  return(list(transcription = t.test(a1,a2), processing = t.test(c1, c2), degradation = t.test(b1, b2)))
})

