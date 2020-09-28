#' esitmate transcription parameters
#'
#' @description It is used to solve the transcription parameters, including the error function and the optimization function.
#' @param TL Exon expression value of 4sU label RNA
#' @param tL label time
#' @param TimeGrid a vector  of times
#' @return a vector of transcription parameters,contains pulse parameters and fitting error values

splitA6Params <- function(TT, TL, tL, TimeGrid, w1, w2){
  ##################
  if(length(TimeGrid)<=5){
    upper_t1 = 2*TimeGrid[length(TimeGrid)]
    upper_t2 = 2*TimeGrid[length(TimeGrid)]
  }else{
    upper_t1 = Inf
    upper_t2 = Inf
  }
  #########################################
  res_A = list(par=rep(NA,6),objective=NA)
  tryCatch({
    res_A = nlminb(start=abs(rnorm(6,0,1)),
                    objective=function(params)errorAFunc(params,list(TimeGrid=TimeGrid, l=TL,tL=tL, TT=TT, w1=w1, w2=w2)),
                    gradient = NULL,
                    #scale = 1, control = list(),
                    lower = c(0, 0.1,0,0,0,0),
                   upper = c(2*TL[1]/tL, Inf,Inf,upper_t1, upper_t2,50))
    },
    error = function(e){}
  )
  return(c(res_A$par,res_A$objective))
}
