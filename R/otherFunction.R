
filter.outers <- function(data, fold=1.5){
  filter_na = as.matrix(na.omit(data))

  filter_na[which(filter_na==Inf)] = 0
  eval_var=apply(filter_na, 1, var)
  fenwei = quantile(eval_var)
  filter_var <- filter_na[(eval_var < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) &
                             eval_var > (fenwei[2] - (fenwei[4]-fenwei[2])*fold)),]
  return(filter_var)
}

filterzeros <- function(data){
  data = na.omit(data)
  zeronumber = c()
  for(i in 1:dim(data)[1]){
    if(length(which(data[i,]==0)) >= 2/3*dim(data)[2]){
      zeronumber = c(zeronumber, i)
    }
  }
  if(length(zeronumber)!=0){
    data=data[-zeronumber,]
  }
  return(data)
}
