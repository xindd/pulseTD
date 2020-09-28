#' esitmate degradation parameters
#'
#' @description It is used to solve the degradation parameters, including the error function and the optimization function.
#' @param TT exon expression value of total RNA
#' @param pre_a Parameter vector of transcription pulse function
#' @param TimeGrid a vector  of times
#' @return a vector of degradation parameters,contains pulse parameters and fitting error values

splitB4Params <- function(TT, pre_a, TimeGrid){

  ##################
  if(length(TimeGrid)<=5){
    upper_t1 = 2*TimeGrid[length(TimeGrid)]
    upper_t2 = 2*TimeGrid[length(TimeGrid)]
  }else{
    upper_t1 = Inf
    upper_t2 = Inf
  }
  res_B = list(par=rep(NA,4),objective=NA)
  tryCatch({
    res_B = nlminb(start=abs(rnorm(4,0,1)),
                   objective=function(params)errorBFunc(params, list(TimeGrid=TimeGrid, pre_a=pre_a, y=TT)),
                   gradient = NULL,
                   #scale = 1, control = list(),
                   lower = c(0.1,0,0,0),
                   upper = c(Inf,upper_t1, upper_t2,50))
  },
  error = function(e){}
  )
  return(c(res_B$par,res_B$objective))
}
