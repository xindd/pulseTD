#' esitmate processing parameters
#'
#' @description It is used to solve the processing parameters, including the error function and the optimization function.
#' @param PT Intron expression value of total RNA
#' @param pre_a Parameter vector of transcription pulse function
#' @param TimeGrid a vector of times
#' @return a vector of processing parameters,contains pulse parameters and fitting error values

splitC5Params <- function(PT, pre_a, TimeGrid){

  ##################
  if(length(TimeGrid)<=5){
    upper_t1 = 2*TimeGrid[length(TimeGrid)]
    upper_t2 = 2*TimeGrid[length(TimeGrid)]
  }else{
    upper_t1 = Inf
    upper_t2 = Inf
  }
  res_C = list(par=rep(NA,5),objective=NA)
  tryCatch({
    res_C = nlminb(start=abs(rnorm(5,0,1)),
                   objective=function(params)errorBCFunc(params, list(TimeGrid=TimeGrid, pre_a=pre_a, y=PT)),
                   gradient = NULL,
                   #scale = 1, control = list(),
                   lower = c(0.1,0,0,0,0),
                   upper = c(Inf,pre_a[3],upper_t1, upper_t2,50))
  },
  error = function(e){}
  )
  return(c(res_C$par,res_C$objective))
}
