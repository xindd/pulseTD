plotSingleGene <- function(parlist, PT, TT, TL,w1, w2, tL, TimeGrid, predict=FALSE){
  if(length(predict)==3){
    pre_time = seq(predict[1],predict[2],predict[3])
  }else{
    pre_time = TimeGrid
  }
  #####################
  pre_a = parlist$pre_a
  pre_c = parlist$pre_c
  pre_b = parlist$pre_b

  par(mfrow = c(2, 3))
  pre_TL = as.numeric(w1 * TT + w2 * integrateA(pre_a,TimeGrid, tL))
  pre_PT = integrateAC(pre_a,pre_time) - integrateAC(pre_c,pre_time) + PT[1]
  pre_TT = integrateAC(pre_a,pre_time) - integrateAC(pre_b,pre_time) + TT[1]

  n=length(pre_time)
  up_PT = c()
  down_PT = c()
  up_TT = c()
  down_TT =c()
  for(i in 1:n){
    X = pre_time[i]
    Yj = pre_PT[i]
    tmp = Yj + qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_PT-Yj)^2) / (n-2)) * sqrt(1/n + (X -mean(pre_time))^2 / ((n-1)*var(pre_time)))
    tmp2 = Yj - qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_PT-Yj)^2) / (n-2)) * sqrt(1/n + (X -mean(pre_time))^2 / ((n-1)*var(pre_time)))
    tmp3 = pre_TT[i] + qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_TT-pre_TT[i])^2) / (n-2)) * sqrt(1/n + (X -mean(pre_time))^2 / ((n-1)*var(pre_time)))
    tmp4 = pre_TT[i] - qt(0.05, n-2, lower.tail = FALSE) * sqrt(sum((pre_TT-pre_TT[i])^2) / (n-2)) * sqrt(1/n + (X -mean(pre_time))^2 / ((n-1)*var(pre_time)))

    up_PT = c(up_PT, tmp)
    down_PT = c(down_PT, tmp2)
    up_TT = c(up_TT, tmp3)
    down_TT = c(down_TT, tmp4)
  }

  gTL = ggplot(show.legend = TRUE)+
    geom_point(aes(x=TimeGrid, y=as.vector(TL),color='true'), size=2)+
    geom_line(aes(x=TimeGrid, y=pre_TL, color='pre'), size=1)+
    scale_color_manual(values=c("#0000FF","black"))+
    labs(title="",y = "Labelexons", x='Time',color = "Type")+
    theme(legend.text = element_text(colour="black", size = 12, face = "bold"),
          legend.title =  element_text(colour="black", size = 12, face = "bold"),
          legend.position='top')+
    theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
          axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
    theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
          axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))

  gPT = ggplot()+
    geom_ribbon(aes(x=pre_time,ymin=down_PT,ymax=up_PT), fill="grey",alpha=0.5)+
    geom_line(aes(x=pre_time, y=pre_PT, color='pre'))+
    geom_line(aes(x=TimeGrid, y=as.vector(PT)))+
    geom_point(aes(x=TimeGrid, y=as.vector(PT), color='true'), size=2)+
    scale_color_manual(values=c("#0000FF","black"))+
    geom_vline(aes(xintercept=TimeGrid[length(TimeGrid)]),linetype=5,col="red")+
    labs(title="",y = "Totalintorns", x='Time',color = "Type")+
    theme(legend.text = element_text(colour="black", size = 12, face = "bold"),
          legend.title =  element_text(colour="black", size = 12, face = "bold"),
          legend.position='top',legend.key = element_blank(),legend.background = element_blank())+
    theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
          axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
    theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
          axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))

  gTT=ggplot()+
    geom_ribbon(aes(x=pre_time,ymin=down_TT,ymax=up_TT), fill="grey",alpha=0.5)+
    geom_line(aes(x=pre_time, y=pre_TT, color='pre'))+
    geom_line(aes(x=TimeGrid, y=as.vector(TT)))+
    geom_point(aes(x=TimeGrid, y=as.vector(TT), color='true'), size=2)+
    scale_color_manual(values=c("#0000FF","black"))+
    geom_vline(aes(xintercept=TimeGrid[length(TimeGrid)]),linetype=5,col="red")+
    labs(title="",y = "Totalexons", x='Time',color = "Type")+
    theme(legend.text = element_text(colour="black", size = 12, face = "bold"),
          legend.title =  element_text(colour="black", size = 12, face = "bold"),
          legend.position='top',legend.key = element_blank(),legend.background = element_blank())+
    theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
          axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
    theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
          axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))

  predict_a = pulseModel(pre_a, pre_time)
  ga = ggplot()+
    geom_line(aes(x=pre_time, y=predict_a, color='alpha'))+
    scale_color_manual(values=c("#0000FF"))+
    geom_vline(aes(xintercept=TimeGrid[length(TimeGrid)]),linetype=5,col="red")+
    labs(title="",y = "Transcription", x='Time',color = "Rates")+
    theme(legend.text = element_text(colour="black", size = 12, face = "bold"),
          legend.title =  element_text(colour="black", size = 12, face = "bold"),
          legend.position='top',legend.key = element_blank(),legend.background = element_blank())+
    theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
          axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
    theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
          axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))


  predict_c = pulseModel(pre_c, pre_time)
  gc = ggplot()+
    geom_line(aes(x=pre_time, y=predict_c, color='gamma'))+
    scale_color_manual(values=c("#0000FF"))+
    geom_vline(aes(xintercept=TimeGrid[length(TimeGrid)]),linetype=5,col="red")+
    labs(title="",y = "Processing", x='Time',color = "Rates")+
    theme(legend.text = element_text(colour="black", size = 12, face = "bold"),
          legend.title =  element_text(colour="black", size = 12, face = "bold"),
          legend.position='top',legend.key = element_blank(),legend.background = element_blank())+
    theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
          axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
    theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
          axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))

  predict_b = pulseModel(pre_b, pre_time)
  gb=ggplot()+
    geom_line(aes(x=pre_time, y=predict_b, color='beta'))+
    scale_color_manual(values=c("#0000FF"))+
    geom_vline(aes(xintercept=TimeGrid[length(TimeGrid)]),linetype=5,col="red")+
    labs(title="",y = "Degradation", x='Time',color = "Rates")+
    theme(legend.text = element_text(colour="black", size = 12, face = "bold"),
          legend.title =  element_text(colour="black", size = 12, face = "bold"),
          legend.position='top',legend.key = element_blank(),legend.background = element_blank())+
    theme(axis.text.x=element_text(face="bold",size=12,angle=0,color="Black"),
          axis.text.y=element_text(face="bold",size=12,angle=0,color="Black"))+
    theme(axis.title.y = element_text(size = 12, face="bold", color = "Black"),
          axis.title.x = element_text(size = 12, face="bold", color = "Black"))+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2,3)))
  print(gTL, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(gPT, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(gTT, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(ga, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(gc, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(gb, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
  popViewport()

}
