#' heatmapParams-Method
#' @description
#' Draw the image of the fitting result.
#' @param object An object of class 'pulseTDmodel'
#' @param genename The name of the gene you want to draw.
#' @return No return value
#'
setMethod('heatmapParams', 'pulseTDmodel', function(object, genename=NULL) {
  if(is.null(genename)){
    heatmapGene(a = as.matrix(pData(object@ratesPar.transcription)),
                c = as.matrix(pData(object@ratesPar.processing)),
                b = as.matrix(pData(object@ratesPar.degradation))
    )
  }else if(!genename %in% object@genenames && !genename  %in% 1:length(object@genenames)){
    stop('Gene name does not existe')
  }else{
    heatmapGene(a = as.matrix(pData(object@ratesPar.transcription)[genename,]),
                c = as.matrix(pData(object@ratesPar.processing)[genename,]),
                b = as.matrix(pData(object@ratesPar.degradation)[genename,])
    )
  }
})

