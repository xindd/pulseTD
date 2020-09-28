#' pulse function
#'
#' @description Pulse function is a pulse function J Comput Biol (2009 Feb) <doi:10.1089/cmb.2008.13TT>.
#' This is a multiplication of two sigmoid functions, the parameter vector Theta = (h0, h1, h2, t1 t2, beta),
#' where h0, h1, h2, represent the initial state rate value, peak value and the steady state rate value is reached again,
#' t1 and t2 are the maximum times of the first and second rise or fall changes, respectively,
#' and beta is the slope of the two changes.
#'
#' @param xita Six parameters of the pulse function :h0, h1, h2, t1, t2, beta. Details \code{\link{getParams}}
#' @param x a vector or point of time
#'
#' @return a vector or point of pulse function value
#' @export
#' @examples
#' load(file.path(system.file(package="pulseTD"),"data","pulseRates.RData"))
#' TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
#' pulseRates_correct <- correctionParams(pulseRates)
#' transcription_params = getParams(pulseRates_correct,'transcription')
#' degradation_params = getParams(pulseRates_correct, 'degradation')
#' processing_params = getParams(pulseRates_correct, 'processing')
#' ###
#' transcription_pulse = pulseModel(as.matrix(transcription_params[1,]), TimeGrid)
#' degradation_pulse = pulseModel(as.matrix(degradation_params[1,]), TimeGrid)
#' processing_pulse = pulseModel(as.matrix(processing_params[1,]), TimeGrid)


pulseModel<-function(xita,x){
  h0 = xita[1]
  h1 = xita[2]
  h2 = xita[3]
  t1 = xita[4]
  t2 = xita[5] + t1
  beta=xita[6]
  ( h0 + (h1-h0)/(1+exp(-beta*(x-t1))) ) * ( h2 + (h1-h2)/(1+exp(beta*(x-t2))) ) / h1
}

#' integral function
#' @description Calculate the exon expression value of the label time.
#' @param xita Six parameters of the pulse function :h0, h1, h2, t1, t2, beta
#' @param timegrid a vector or point of time
#' @param tl label time
#' @return a vector or point of exon expression value
integrateA <- function(xita, timegrid, tl){
  res=sapply(timegrid,function(xt)(integrate(function(tt)pulseModel(xita,tt), lower=xt-tl, upper=xt,abs.tol=1e-100)$value))
  return(res)
}

#' integral function
#'@description Calculate the expression values at different time points.
#' @param xita Six parameters of the pulse function :h0, h1, h2, t1, t2, beta
#' @param timegrid a vector or point of time
#' @return a vector or point of exon expression value
integrateAC <- function(xita, timegrid){
  res=sapply(timegrid,function(xt)(integrate(function(tt)pulseModel(xita,tt), lower=0, upper=xt,abs.tol=1e-100)$value))
  return(res)
}
