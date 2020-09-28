#' Get pulseModel params
#' @description
#' It is used to obtain the parameters of the pulse function J Comput Biol (2009 Feb) <doi:10.1089/cmb.2008.13TT>, including transcription parameters, degradation parameters and processing parameters.
#' You can use these parameters to analyze transcriptional characteristics, including steady-state analysis.
#' At the same time, you can use \code{\link{pulseModel}} to view the curve corresponding to the parameter.
#'
#' @param object a 'pulseTDmodel' that has been calculated with \code{\link{estimateParams}}
#' @param stage A character is one of three stages of transcriptional dynamics, transcription, processing and degradation.
#' @param genename a vector, default is NULL
#' \itemize{
#'  \item If it is NULL: Calculate the transcriptional dynamic rate of all genes
#'  \item If it is a gene vector: only calculate the transcriptional dynamic rate of a given gene vector
#'  \item If it is a numerical vector: only calculate the transcriptional dynamic rate of the gene corresponding to a given value.
#' }
#' @return A matrix or vector containing six parameters: h0,h1,h2,t1 t2,beta
#' @export
#' @examples
#' data('pulseRates', package='pulseTD')
#' TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
#' pulseRates_correct <- correctionParams(pulseRates)
#' transcription_params = getParams(pulseRates,'transcription')
#' degradation_params = getParams(pulseRates, 'degradation')
#' processing_params = getParams(pulseRates, 'processing')
#' head(transcription_params)
#' head(degradation_params)
#' head(processing_params)
#' transcription_params = getParams(pulseRates, 'transcription', genename=c(1,2,3))
#' head(transcription_params)
#' ###
#' transcription_pulse = pulseModel(as.matrix(transcription_params[1,]), TimeGrid)
#' degradation_pulse = pulseModel(as.matrix(degradation_params[1,]), TimeGrid)
#' processing_pulse = pulseModel(as.matrix(processing_params[1,]), TimeGrid)
#'
setGeneric('getParams',
           function(object, stage, genename=NULL) {
             standardGeneric('getParams')
           }
)


#' Get rate values
#' @description
#' It is used to calculate the transcriptional dynamic rate of a gene.
#' You can get the discrete or continuous rate values of the measurement time points.
#' At the same time, it has a predictive function that provides rate values for any future time node or any range of time.
#'
#' @param object a 'pulseTDmodel' that has been calculated with \code{\link{estimateParams}}
#' @param stage A character is one of three stages of transcriptional dynamics, transcription, processing and degradation.
#' @param timevector A vector of time, which can be any time, defaults to NULL
#' \itemize{
#'  \item If it is NULL: only calculate the transcriptional dynamics of the measurement time node, the time is obtained by t_time in pulseTDmodel
#'  \item If the time vector is positive: calculate the rate value for a given time, or the predicted rate if the time exceeds the measurement time node
#'  \item If the time vector is negative: Calculate the rate value at a given time, which represents the calculation of the transcriptional dynamic rate before the prediction of the experimental measurement.
#' }
#' @param genename a vector, default is NULL
#' \itemize{
#'  \item If it is NULL: Calculate the transcriptional dynamic rate of all genes
#'  \item If it is a gene vector: only calculate the transcriptional dynamic rate of a given gene vector
#'  \item If it is a numerical vector: only calculate the transcriptional dynamic rate of the gene corresponding to a given value.
#' }
#' @return A matrix, each row representing a gene, each column with a table to calculate the time node, and the value representing the rate of transcriptional dynamics.
#' @export
#' @examples
#' load(file.path(system.file(package="pulseTD"),"data","pulseRates.RData"))
#' pulseRates_correct <- correctionParams(pulseRates)
#' transcription = getRates(pulseRates_correct,'transcription')
#' degradation = getRates(pulseRates_correct, 'degradation')
#' processing = getRates(pulseRates_correct, 'processing')
#' head(transcription)
#' head(degradation)
#' head(processing)
#' trans = getRates(pulseRates_correct, 'transcription', timevector=c(0,1,2,3))
#' head(transcription)
#' trans=getRates(pulseRates_correct, 'transcription', timevector=c(2,3),genename=c(2,3))
#' head(transcription)
setGeneric('getRates',
           function(object, stage, timevector=NULL, genename=NULL) {
             standardGeneric('getRates')
           }
)

#' draw Method
#' @description
#' It is used to draw a rate image that contains six images.
#' The expression values of the real genes and the fitted expression values are plotted separately, and the rate of the change is varied with time.
#' Can also plot predicted expression values and rate values
#'
#' @param object a pulseTDmodel that has been calculated with \code{\link{estimateParams}}
#' @param genename The name of the gene to be drawn
#' @param predict The time interval you want to draw,default is FALSE, or a vector of three values,(start, end, step)
#' @return There is no return value, six pictures will be drawn separately.
#' @import ggplot2
#' @import grid
#' @export
#' @examples
#' data('pulseRates', package='pulseTD')
#' plotRates(pulseRates, 'NM_001001181')
#' plotRates(pulseRates, 'NM_001001181', predict=c(0,180,20))
setGeneric('plotRates',
           function(object, genename, predict=FALSE) {
             standardGeneric('plotRates')
           }
)


#' predict Expression
#' @description
#' It is used to predict the expression of all gene at a given time, including the expression of pre-mRNA (Precursor RNA) and the expression of total mRNA (Mature RNA).
#' End time and time interval can be arbitrarily defined
#'
#' @param object a 'pulseTDmodel' that has been calculated with \code{\link{estimateParams}}
#' @param tg A vector of points in time at which experimental data is collected, not allowed to be repeated
#' @return Returns a list containing predicted values for each gene and a 0.95 confidence interval
#' @export
#' @examples
#' data('pulseRates', package='pulseTD')
#' pulseRates_correct = correctionParams(pulseRates)
#' TimeGrid = seq(0,180,15)
#' \donttest{preExp = predictExpression(pulseRates_correct, tg=TimeGrid)}
#' data('preExp', packages='pulseTD')
#' df = data.frame(preExp[['NM_001002011']])
#' head(df)


setGeneric('predictExpression',
           function(object, tg) {
             standardGeneric('predictExpression')
           }
)

#' heatmapParams Method
#' @description
#' It is used to draw a paramaters image.
#' @param object a pulseTDmodel that has been calculated with \code{\link{estimateParams}}
#' @param genename The name of the gene to be drawn
#' @return There is no return value.
#' @import ComplexHeatmap
#' @export
#' @examples
#' data('pulseRates', package='pulseTD')
#' heatmapParams(pulseRates, 1:10)
setGeneric('heatmapParams',
           function(object, genename=NULL) {
             standardGeneric('heatmapParams')
           }
)

#' compareRates Method
#' @description
#' compare result.
#' @param object1 An object of class 'pulseTDmodel'
#' @param object1 An object of class 'pulseTDmodel'
#' @return A list
#' @export
#' @examples
#' data('pulseRates', package='pulseTD')
setGeneric('compareRates',
           function(object1, object2) {
             standardGeneric('compareRates')
           }
)

