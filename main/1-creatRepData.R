####################
rpkms_rep1 = list(foursu_exons=c(),
                  total_exons=c(),
                  foursu_introns=c(),
                  total_introns=c()
)

rpkms_rep2 = list(foursu_exons=c(),
                  total_exons=c(),
                  foursu_introns=c(),
                  total_introns=c()
)

############读取表达数据####################
# 读取表达文件
dirpath = 'D:/output'
filename = list.files(dirpath)
rep1_geneid=c()
rep2_geneid=c()
for(name in filename){
  if(substr(name,12,12)=="1"){
    tmpfile = read.table(paste(dirpath,name,sep='/'), header = TRUE)
    tmpfile = tmpfile[order(tmpfile$id),]
    rpkms_rep1$foursu_exons=cbind(rpkms_rep1$foursu_exons,tmpfile$foursu_exons)
    rpkms_rep1$foursu_introns=cbind(rpkms_rep1$foursu_introns,tmpfile$foursu_introns)
    rpkms_rep1$total_exons=cbind(rpkms_rep1$total_exons,tmpfile$total_exons)
    rpkms_rep1$total_introns=cbind(rpkms_rep1$total_introns,tmpfile$total_introns)

  }
  if(substr(name,12,12)=="2"){
    tmpfile = read.table(paste(dirpath,name,sep='/'), header = TRUE)
    tmpfile = tmpfile[order(tmpfile$id),]
    rpkms_rep2$foursu_exons=cbind(rpkms_rep2$foursu_exons,tmpfile$foursu_exons)
    rpkms_rep2$foursu_introns=cbind(rpkms_rep2$foursu_introns,tmpfile$foursu_introns)
    rpkms_rep2$total_exons=cbind(rpkms_rep2$total_exons,tmpfile$total_exons)
    rpkms_rep2$total_introns=cbind(rpkms_rep2$total_introns,tmpfile$total_introns)
  }
}
rownames(rpkms_rep1$foursu_exons) = tmpfile$id
rownames(rpkms_rep1$foursu_introns) = tmpfile$id
rownames(rpkms_rep1$total_exons) = tmpfile$id
rownames(rpkms_rep1$total_introns) = tmpfile$id
rownames(rpkms_rep2$foursu_exons) = tmpfile$id
rownames(rpkms_rep2$foursu_introns) = tmpfile$id
rownames(rpkms_rep2$total_exons) = tmpfile$id
rownames(rpkms_rep2$total_introns) = tmpfile$id

colnames(rpkms_rep1$foursu_exons) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep1$foursu_introns) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep1$total_exons) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')
colnames(rpkms_rep1$total_introns) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')
colnames(rpkms_rep2$foursu_exons) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep2$foursu_introns) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep2$total_exons) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')
colnames(rpkms_rep2$total_introns) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')

############提取部分基因####################
## filter zeros expression; filter expression < 0.5
filterzeros <- function(data){
  rownum = c()
  for(i in 1:dim(data)[1]){
    if(length(which(data[i,]==0)) >= 5){
      rownum = c(rownum, i)
    }else if(length(which(data[i,]<=0.1)) >=5){
      rownum = c(rownum, i)
    }
  }
  data=data[-rownum,]
}

# fielter error data , must have TL<TT, PT<TT
filtererror <- function(TL, TT, PT){
  w = dim(TL)[2]
  h = dim(TL)[1]
  rownum = c()
  for(i in i:h){
    if(sum(TL[i,]<TT[i,]) < w){
      rownum = c(rownum, i)
    }else if(sum(PT[i,]<TT[i,]) < w){
      rownum = c(rownum, i)
    }
  }
  data = list()
  data$tl = TL[-rownum, ]
  data$tt = TT[-rownum, ]
  data$pt = PT[-rownum, ]
  return(data)
}

labexon = rpkms_rep1$foursu_exons
labintr = rpkms_rep1$foursu_introns
totexon = rpkms_rep1$total_exons
totintr = rpkms_rep1$total_introns

labexon = filterzeros(labexon)
labintr = filterzeros(labintr)
totexon = filterzeros(totexon)
totintr = filterzeros(totintr)

labexon2 = rpkms_rep2$foursu_exons
labintr2 = rpkms_rep2$foursu_introns
totexon2 = rpkms_rep2$total_exons
totintr2 = rpkms_rep2$total_introns

labexon2 = filterzeros(labexon2)
labintr2 = filterzeros(labintr2)
totexon2 = filterzeros(totexon2)
totintr2 = filterzeros(totintr2)


genelist2 = intersect(rownames(labexon2),rownames(labintr2))
genelist2 = intersect(genelist2,rownames(totexon2))
genelist2 = intersect(genelist2,rownames(totintr2))


genelist = intersect(rownames(labexon),rownames(labintr))
genelist = intersect(genelist,rownames(totexon))
genelist = intersect(genelist,rownames(totintr))

genelist = intersect(genelist, genelist2)

labexon = labexon[genelist,]
labintr = labintr[genelist,]
totexon = totexon[genelist,]
totintr = totintr[genelist,]
labexon = labexon[order(rownames(labexon),decreasing = TRUE),]
labintr = labintr[order(rownames(labintr),decreasing = TRUE),]
totexon = totexon[order(rownames(totexon),decreasing = TRUE),]
totintr = totintr[order(rownames(totintr),decreasing = TRUE),]

labexon2 = labexon2[genelist,]
labintr2 = labintr2[genelist,]
totexon2 = totexon2[genelist,]
totintr2 = totintr2[genelist,]
labexon2 = labexon2[order(rownames(labexon2),decreasing = TRUE),]
labintr2 = labintr2[order(rownames(labintr2),decreasing = TRUE),]
totexon2 = totexon2[order(rownames(totexon2),decreasing = TRUE),]
totintr2 = totintr2[order(rownames(totintr2),decreasing = TRUE),]

rpkmRep1=list()
rpkmRep2=list()
rpkmRep1$labexon=labexon[1:1000,]
rpkmRep1$labintr=labintr[1:1000,]
rpkmRep1$totexon=totexon[1:1000,]
rpkmRep1$totintr=totintr[1:1000,]

rpkmRep2$labexon2=labexon2[1:1000,]
rpkmRep2$labintr2=labintr2[1:1000,]
rpkmRep2$totexon2=totexon2[1:1000,]
rpkmRep2$totintr2=totintr2[1:1000,]
save(rpkmRep1, file = "rpkmRep1.RData")
save(rpkmRep2, file = "rpkmRep2.RData")

###############################################################################
t_time <- c(0, 2, 4, 6, 8)
tL <- 1
gnum = 10
rep1 = estimateParams(labexon[1:gnum,],totexon[1:gnum,],totintr[1:gnum,], t_time, tL, loopnumber=100)
rep2 = estimateParams(labexon2[1:gnum,],totexon2[1:gnum,],totintr2[1:gnum,], t_time, tL, loopnumber=100)
save(rep1, file = "ratesRep1.RData")
save(rep2, file = "ratesRep2.RData")

TT_PT = c()
for(i in 1:1000){
  if(sum(rpkmRep1$totintr[i,]<rpkmRep1$totexon[i,])==5){
    TT_PT = c(TT_PT,i)
  }
}
length(TT_PT)/1000

