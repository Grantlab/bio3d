communities.comparison <- function(networks=networks,
                                   sse=NULL,
                                   dash.lines=FALSE){
  
  memb.matrix <- rbind(networks[[1]]$raw.communities$membership,
                       networks[[2]]$raw.communities$membership,
                       networks[[3]]$raw.communities$membership,
                       networks[[4]]$raw.communities$membership,
                       networks[[5]]$raw.communities$membership,
                       networks[[6]]$raw.communities$membership,
                       networks[[7]]$raw.communities$membership,
                       networks[[8]]$raw.communities$membership,
                       networks[[9]]$raw.communities$membership,
                       networks[[10]]$raw.communities$membership,
                       networks[[11]]$raw.communities$membership,
                       networks[[12]]$raw.communities$membership,
                       networks[[13]]$raw.communities$membership,
                       networks[[14]]$raw.communities$membership,
                       networks[[15]]$raw.communities$membership)

  rownames(memb.matrix) <- c("btwn.nocmap.cij0.4",
                           "btwn.nocmap.cij0.6",
                           "btwn.cmap.cij0.4",
                           "btwn.cmap.cij0.6",
                           "btwn.cmap.cij0.2",
                           "walk.nocmap.cij0.4",
                           "walk.nocmap.cij0.6",
                           "walk.cmap.cij0.4",
                           "walk.cmap.cij0.6",
                           "walk.cmap.cij0.2",
                           "greed.nocmap.cij0.4",
                           "greed.nocmap.cij0.6",
                           "greed.cmap.cij0.4",
                           "greed.cmap.cij0.6",
                           "vmd.cmap.cij0.4")
  
  alpha4.memb <- membership.comparison(memb.matrix, important.members=c(168:180,275:312))
 
  loop7.memb <- membership.comparison(memb.matrix,important.members=c(149:152))
  
  beta.strand.memb <- membership.comparison(memb.matrix,important.members=c(75:81,130:148,221:249))
   
  ploop.memb <- membership.comparison(memb.matrix,important.members=c(90:96))
  
  switch1.memb <- membership.comparison(memb.matrix,important.members=c(206:220))
    
  switch2.memb <- membership.comparison(memb.matrix,important.members=c(250:274))
  

  loop7.memb[loop7.memb >0] <- 2
  beta.strand.memb[beta.strand.memb >0] <- 3
  ploop.memb[ploop.memb >0] <- 4
  switch1.memb[switch1.memb >0] <- 5
  switch2.memb[switch2.memb >0] <- 6

  
  total.membership <- matrix(0, nrow=nrow(alpha4.memb), ncol=ncol(alpha4.memb))

  for(h in 1:dim(alpha4.memb)[1]){
    for(k in 1:dim(alpha4.memb)[2]){
      total.membership[h,k] <- max(alpha4.memb[h,k],
                                   loop7.memb[h,k],
                                   beta.strand.memb[h,k],
                                   ploop.memb[h,k],
                                   switch1.memb[h,k],
                                   switch2.memb[h,k])
    }
  }
  
  rownames(total.membership) <- rownames(alpha4.memb)
  col=c("white","green","yellow","cyan","red","blue","orange")

  par(mar=c(5,9,1,0.5))
  image(1:ncol(total.membership), 1:nrow(total.membership), t(total.membership), col=col, ylab="", axes=FALSE, xlab="Residue No", cex.lab=1.2)
  axis(BELOW<-1, at=seq(1,ncol(total.membership),25), labels=seq(1,ncol(total.membership),25), cex.axis=0.8)
  axis(LEFT <-2, at=seq(1,nrow(total.membership),1),labels=rownames(total.membership),
       las= HORIZONTAL<-1,
       cex.axis=0.8
       )

  if(!is.null(sse)){
    sse.draw(sse, ylim=c(0,dim(total.membership)[1]))
  }
  
  if(dash.lines){
    for(b in 1:(nrow(total.membership)-1)){
      abline(h=(b+0.5), lty=2)
    }
  }

  return(total.membership)
}
  
