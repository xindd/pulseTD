#######加载数据###############
dirpath = 'E:\\pluseTD\\RASG\\data'
rpkms_rep1=list()
rpkms_rep1$foursu_exons = read.table(paste(dirpath, '\\expression_gene.4sU.M.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$foursu_introns = read.table(paste(dirpath, '\\expression_gene.4sU.P.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$total_exons = read.table(paste(dirpath, '\\expression_gene.total.M.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$total_introns = read.table(paste(dirpath, '\\expression_gene.total.P.txt',sep=''),sep='\t', header = TRUE)

rownames(rpkms_rep1$foursu_exons) = rpkms_rep1$foursu_exons[,1]
rownames(rpkms_rep1$foursu_introns) = rpkms_rep1$foursu_introns[,1]
rownames(rpkms_rep1$total_exons) = rpkms_rep1$total_exons[,1]
rownames(rpkms_rep1$total_introns) = rpkms_rep1$total_introns[,1]

rpkms_rep1$foursu_exons = rpkms_rep1$foursu_exons[,-1]
rpkms_rep1$foursu_introns = rpkms_rep1$foursu_introns[,-1]
rpkms_rep1$total_exons = rpkms_rep1$total_exons[,-1]
rpkms_rep1$total_introns = rpkms_rep1$total_introns[,-1]

labexon = rpkms_rep1$foursu_exons
labintr = rpkms_rep1$foursu_introns
totexon = rpkms_rep1$total_exons
totintr = rpkms_rep1$total_introns

genelist = intersect(rownames(labexon),rownames(totexon))
genelist = intersect(genelist,rownames(totintr))
genelist = intersect(genelist,rownames(labintr))

labexon = labexon[genelist,]
labintr = labintr[genelist,]
totexon = totexon[genelist,]
totintr = totintr[genelist,]

rpkmRep=list()
rpkmRep$labexon=labexon[1:1000,]
rpkmRep$labintr=labintr[1:1000,]
rpkmRep$totexon=totexon[1:1000,]
rpkmRep$totintr=totintr[1:1000,]

save(rpkmRep, file = "rpkmTRUE.RData")
######生成仿真#######################################
TimeGrid <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
gnum = 1000
#计算速率
resmodel20 = estimateParams(labexon[1:gnum,],
                            totexon[1:gnum,],
                            totintr[1:gnum,], TimeGrid, tL, loopnumber=20)
plotRates(resmodel20,2)

################################
resmodel20 = correctionParams(resmodel20)

genename = resmodel20@genenames
a20 = getRates(resmodel20, 'transcription', timevector = TimeGrid)
b20 = getRates(resmodel20, 'degradation', timevector = TimeGrid)/as.matrix(totexon[genename,]-totintr[genename,])
c20 = getRates(resmodel20, 'processing', timevector = TimeGrid)/as.matrix(totintr[genename,])
ratesTrue = list(a=a20,b=b20, c=c20)
save(ratesTrue, file = "ratesTRUE.RData")
#龙哥库塔计算表达至
at = a20
bt = b20
ct = c20
# 计算T，P，TL初始值
sim_P0 = totintr[genename,1]
sim_T0 = totexon[genename,1]
sim_PT = c()
sim_TT = c()
sim_TL = at*tL
t = TimeGrid
for(i in 1:dim(at)[1]){
  parms  <- cbind(spline(TimeGrid,at[i,], length(t))$y,
                  spline(TimeGrid,bt[i,], length(t))$y,
                  spline(TimeGrid,ct[i,], length(t))$y,
                  t)
  P_T = deSolve::rk(c(P=sim_P0[i], C=sim_T0[i]), t, RungFunction2, parms)
  sim_PT = rbind(sim_PT, P_T[,2])
  sim_TT = rbind(sim_TT, P_T[,3])
}
rownames(sim_PT) = rownames(at)
rownames(sim_TT) = rownames(at)
rownames(sim_TL) = rownames(at)
colnames(sim_PT) = t
colnames(sim_TT) = t
colnames(sim_TL) = colnames(at)
sim_TL_sample = sim_TL
sim_PT_sample = sim_PT
sim_TT_sample = sim_TT

rpkmSim = list()
rpkmSim$labexon = sim_TL_sample
rpkmSim$totintr = sim_PT_sample
rpkmSim$totexon = sim_TT_sample
save(rpkmSim, file = "rpkmSim.RData")





