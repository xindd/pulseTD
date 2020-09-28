#' pulseTD-model
#'
#' @description 'pulseTD-model' is S4 class for storing calculation results.
#' It includes pulse model parameters corresponding to transcription rate, processing rate and degradation rate, and also stores some basic parameters of the calculation process, such as gene name, sequencing time, labeling time and so on.
#‘ It can be used as input to other functions, such as \code{\link{getRates}}.

#' @slot genenames All gene names entered by the user.
#' @slot t_time Storage time series.
#' @slot tL 4sU labeled duration
#' @slot fitfailure Gene vector for parameter optimization failure
#' @slot filterExpression.foursu_exons 'AnnotatedDataFrame' that stores 4sU-labeled exons expression values
#' @slot filterExpression.total_exons 'AnnotatedDataFrame' that stores total exons expression values
#' @slot filterExpression.total_introns 'AnnotatedDataFrame' that stores total introns expression values
#' @slot ratesPar.degradation 'AnnotatedDataFrame' that stores degradation pulse parameters
#' @slot ratesPar.processing 'AnnotatedDataFrame' that stores processing pulse parameters
#' @slot ratesPar.transcription 'AnnotatedDataFrame' that stores transcription pulse parameters
#' @slot ratesPar.score 'AnnotatedDataFrame' that stores the fit scores
#' @slot ratesPar.w 'AnnotatedDataFrame' that normalization factors
#' @exportClass pulseTDmodel
pulseTDmodel<-setClass('pulseTDmodel',
                       slots=c(
                         ratesPar=list(
                           transcription='AnnotatedDataFrame',
                           processing='AnnotatedDataFrame',
                           degradation='AnnotatedDataFrame',
                           score='AnnotatedDataFrame'

                         ),
                         filterExpression=list(
                           total_introns='AnnotatedDataFrame',
                           total_exons='AnnotatedDataFrame',
                           foursu_exons='AnnotatedDataFrame'
                         ),
                         genenames='vector',
                         t_time='vector',
                         tL='numeric',
                         fitfailure='vector',
                         w1='vector',
                         w2='vector',
                         chisq='vector'
                       ),
                       prototype=list(
                         fitfailure=NULL
                       )
)

