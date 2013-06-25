#' Visualizes the effect of clustering on p-values
#'
#' See above
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
visualize <- function(rvobj,label.by,phen.label='',top=10) {

  # Prepare the cluster p-values
  temp.df <- rvobj$clusterinfo[,c('CLUSTERID','PVALUE')]
  temp.df <- temp.df[order(temp.df$PVALUE),]
  print('cluster p-values (ordered')
  print(temp.df$PVALUE)
  clus.sig.thresh <- 0.05 / nrow(temp.df)
  max <- min(top,nrow(temp.df))
  temp.df <- temp.df[1:max,]
  clus.labels <- sapply(temp.df$CLUSTERID,function(x){paste("Cluster #",x,sep='')})
  temp.df[,label.by] <- clus.labels
  temp.df$COLOR="orange"

  rank <- data.frame(CLUSTERID=temp.df$CLUSTERID,RANK=1:max)

  # Prepare the variant p-values
  df <- rvobj$variants[,c('CLUSTERID','PVALUE',label.by)]
  var.sig.thresh <- 0.05 / nrow(df)
  df <- df[df$CLUSTERID %in% temp.df$CLUSTERID,]
  df$COLOR <- "blue"
  df <- rbind(df,temp.df)

  print("Cluster significance threshold:")
  print(clus.sig.thresh)
  print("Variant significance threshold:")
  print(var.sig.thresh)

  # Order all clusters and variants by the wrt cluster ranking
  df$RANK <- sapply(df$CLUSTERID,function(x,rank){rank[rank$CLUSTERID==x,"RANK"]},rank=rank)
  
  # Order the p-values
  df <- df[order(df$RANK,df$PVALUE),]

  # Plot the p-values
  pdf(paste('barplot_',phen.label,'.pdf',sep=''),height=9,width=10)
  #bplot <- barplot(-log(df$PVALUE),horiz=TRUE,col=df$COLOR,names.arg=df[,label.by],main="The Effect of Clustering on Rare Variant Significance",xlab="Inv Log P-Value",ylab=label.by)
  bplot <- barplot(-log(df$PVALUE),col=df$COLOR,las=2,cex.names=0.5,names.arg=df[,label.by],main="The Effect of Clustering on Rare Variant Significance",xlab=label.by,ylab="Inv Log P-Value")
  lines(x=bplot,y=rep(-log(var.sig.thresh),nrow(df)),col="blue")
  lines(x=bplot,y=rep(-log(clus.sig.thresh),nrow(df)),col="orange")
  print(bplot)
  dev.off()

  return(rvobj)
}
