#' Estimated Expression
#'
#' This function is used to calculate the expression value of the gene.
#' It includs total exon expression values, total intron expression values, labeled exon expression values, and labeled intron expression values.
#'
#' @param txdb Genomic annotation file
#' @param filelist A vector of sam file paths containing 4sU labeled sam files and unlabeled files
#' @param by In terms of genes or transcripts, the default is the gene
#' @param format In terms of RPKM, FPKM or TPM, the default is the RPKM
#' @param mode mode can be one of the pre-defined count methods such as "Union", "IntersectionStrict", or "IntersectionNotEmpty" or it a user supplied count function.
#' @param singleEnd (Default TRUE) A logical indicating if reads are single or paired-end.
#' @param ignore.strand A logical indicating if strand should be considered when matching
#' @param fragments (Default FALSE) A logical; applied to paired-end data only. fragments control which function is used to read the data which subsequently affects which records are included in counting.
#' @param ... The parameters of the \code{\link{summarizeOverlaps}} method
#' @import GenomicAlignments
#' @import GenomicFeatures
#' @importFrom Rsamtools BamFileList
#' @import AnnotationDbi
#' @importFrom SummarizedExperiment assay
#' @export
#' @return A list of expression
#' @seealso Other parameters \code{\link{summarizeOverlaps}}
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' test_path <- file.path(system.file(package="pulseTD"),'extdata/test1.sorted.bam')
#' test_path2 <- file.path(system.file(package="pulseTD"),'extdata/test2.sorted.bam')
#' \donttest{rpkmres <- estimateExpression(txdb,c(test_path,test_path2), by='gene')}
#' data('rpkmres', package='pulseTD')
#' head(rpkmres$total_exp)
#' head(rpkmres$pre_exp)

estimateExpression<-function(txdb, filelist, by='gene', format='RPKM',
                             mode = "Union",singleEnd=TRUE,ignore.strand=TRUE,fragments=FALSE,...){
  message('Starting...')
  if(by=="gene"){
    ## Get the exons grouped by gene:
    exonsByGene = exonsBy(txdb, "gene")
    ## Get the introns grouped by gene:
    introns <- intronsByTranscript(txdb, use.names=TRUE)
    ulst <- unlist(introns)
    intronsNoDups <- ulst[!duplicated(ulst)]
    txnames <- names(intronsNoDups)
    maptg <- select(txdb, keys=txnames, keytype='TXNAME', columns='GENEID')
    idx <- maptg$GENEID[!is.na(maptg$GENEID)]
    intronsByGene <- S4Vectors::split(intronsNoDups[!is.na(maptg$GENEID)], idx)
    names(intronsByGene) <- unique(idx)
  }else if(by=="tx"){
    ## Get the exons grouped by gene:
    exonsByGene = exonsBy(txdb, "tx")
    ## Get the introns grouped by gene:
    intronsByGene <- intronsByTranscript(txdb, use.names=TRUE)
  }else{
    stop('chose \"by\" from  \"gene\" or \"tx\"')
  }

  ##
  message('import bam files')
  bamfiles <- BamFileList(filelist, yieldSize=200000)

  message("Calculating the number of statistical exons reads")
  exons_se <- summarizeOverlaps(features = exonsByGene, reads = bamfiles,
                                mode = "Union",
                                singleEnd=TRUE,
                                ignore.strand=TRUE,
                                fragments=FALSE,...)

  message("Calculating the number of statistical introns reads ")

  introns_se <- summarizeOverlaps(features = intronsByGene, reads = bamfiles,
                                  mode = "Union",
                                  singleEnd=TRUE,
                                  ignore.strand=TRUE,
                                  fragments=FALSE,...)
  intronsPerGene = as.matrix(assay(introns_se))
  introns_width=sapply(width(intronsByGene), sum)

  exonsPerGene = as.matrix(assay(exons_se))
  exons_width= sapply(width(exonsByGene), sum)

  total_reads = apply(intronsPerGene, 2, sum)+apply(exonsPerGene, 2, sum)

  if(format == 'RPKM'){
    #RPKM
    message("Calculating RPKM value")
    exp_introns = intronsPerGene*10^9/total_reads/introns_width
    exp_exons = exonsPerGene*10^9/total_reads/exons_width
  }else if(format == 'FPKM'){
    #FPKM
    message("Calculating FPKM value")
    exp_introns = intronsPerGene*10^9/total_reads/introns_width
    exp_exons = exonsPerGene*10^9/total_reads/exons_width

  }else if(format == 'TPM'){
    #TPM
    print('tpm')
    message("Calculating TPM value")
    exp_introns = intronsPerGene*10^6/introns_width / sum(intronsPerGene/introns_width)
    exp_exons = exonsPerGene*10^6/exons_width / sum(exonsPerGene/exons_width)
  }else{
    message("Please enter the correct format")
  }

  message("Running completed")
  return(list(total_exp = exp_exons, pre_exp=exp_exons))
}

