#' esitmate the global normalization factor
#'
#' @description It is used to the global normalization factor.
#' @param labexon Exon expression value of 4sU label RNA
#' @param totexon Exon expression value of total RNA
#' @param TimeGrid a vector  of times
#' @param tL label time
#' @param batch Number of genes evaluated simultaneously
#' @param samumber Number of samples randomly drawn per iteration
#' @param iternumber Number of iterations per sample
#' @return a vector of the global normalization factor

normFactors <- function(labexon, totexon, TimeGrid, tL, batch = 2, samumber=50, iternumber=50){
  geneNumbers = dim(labexon)[1]
  w_l = length(TimeGrid)
  lo = c(rep(0,w_l),
         rep(0.9, w_l),
         rep(c(0, 0.1,0,0,0,0), batch)
  )
  if(length(TimeGrid)<=5){
    upper_t1 = 2*TimeGrid[length(TimeGrid)]
    upper_t2 = 2*TimeGrid[length(TimeGrid)]
  }else{
    upper_t1 = Inf
    upper_t2 = Inf
  }

  final_w = c()

  for(s in 1:samumber){
    sam = sample(1:geneNumbers, batch, replace = TRUE)
    TL = as.matrix(labexon[sam,])
    TT = as.matrix(totexon[sam,])
    up_tl = c()
    for(i in 1:batch){
      up_tl = c(up_tl, c(2*TL[i,1]/tL, Inf,Inf,upper_t1, upper_t2,50))
    }
    up = c(rep(0.1, w_l),
           rep(1, w_l),
           up_tl
    )
    x=list(abc=list(tL=tL,TimeGrid=TimeGrid), TT=TT, TL=TL,BatchSize=batch)
    init_w = rnorm(2*w_l, 0, 1)
    for(i in 1:iternumber){
      W = c()
      st=c(init_w, rnorm(batch*6, 0,1))
      res = nlminb(start=st,
                   objective=function(params)factorError(x,params),
                   lower = lo,
                   upper = up)
      W = rbind(W, c(res$par[1:(2*w_l)], res$objective))
      init_w = W[which(W[,(2*w_l+1)]==min(W[,(2*w_l+1)]))[1], 1:(2*w_l)]
    }
    final_w = rbind(final_w, init_w)
    print(s)
  }
  return(final_w)
}
