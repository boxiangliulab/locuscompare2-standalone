# input file first column should be "gene_id", rest columns should be p-value of different tools, no other columns are allowed
# results are ranking, the smaller the better
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments are rquired: input_file_path, output_file_path")
}
input_file_path = args[1]
output_file_path = args[2]
if (!file.exists(input_file_path) || !file.info(input_file_path)$size > 0) {
  print("input_file_path file does not exist or is empty!")
  q(save = "no")
}


BIRRA=function(data, prior=0.05, num.bin=50, num.iter=10, return.all=F, plot.legend=F, grp=NULL, cor.stop=1, ...){
  nr=nrow(data)
  nrp=floor(nrow(data)*prior)
  data=apply(data,2,rank)/nr
  
  nc=ncol(data)
  TPR=FPR=Bayes.factors=matrix(ncol=nc, nrow=num.bin)
  binned.data=ceiling(data*num.bin)
  
  
  bayes.data=matrix(nrow=nrow(data), ncol=ncol(data))
  
  guess=apply(data,1,mean)
  cprev=0
  #par(mfrow=c(floor(sqrt(num.iter)), ceiling(sqrt(num.iter))), mai=rep(0.7,4))
  for ( iter in 1:num.iter){
    if((cor.stop-cprev)>1e-15){
      guesslast=guess
      oo=order(guess)
      guess[oo[1:nrp]]=1
      guess[oo[(nrp+1):nr]]=0
      
      
      for (i in 1:nc){  
        for (bin in 1:num.bin){
          frac=bin/num.bin
          TPR=sum(guess[binned.data[,i]<=bin])
          FPR=sum((!guess)[binned.data[,i]<=bin])
          
          Bayes.factors[bin,i]=log((TPR+1)/(FPR+1)/(prior/(1-prior)))
          
        }
      }
      
      Bayes.factors=apply(Bayes.factors,2,smooth)
      Bayes.factors=apply(Bayes.factors,2,function(x){rev(cummax(rev(x)))})
      # Plot TPR vs bin for each data set
      # if(is.null(grp)){
      #   matplot(1:num.bin, Bayes.factors, type="l", lwd=2, ...)
      # }
      # else{
      #   matplot(1:num.bin, Bayes.factors, type="l", lwd=2, lty=grp, col=grp)
      # }
      
      # title(paste("Iteration", iter))
      # if (iter==1&plot.legend){
      #   legend("topright", col=1:5, lty=1:4, legend=colnames(data), lwd=2, ncol=2)
      # }
      for (bin in 1:num.bin){
        oo=order(Bayes.factors[bin,], decreasing=T)
        Bayes.factors[bin, oo[1]]=Bayes.factors[bin, oo[2]]
        
        
      }
      
      for (i in 1:nc){
        
        bayes.data[,i]=Bayes.factors[binned.data[,i],i]
        
      }
      
      
      bb=exp(apply(bayes.data,1, sum))
      f=prior/(1-prior)
      prob=bb*f/(1+bb*f)
      exp=sort(prob, decreasing=F)[nrp]
      
      guess=rank(-apply(bayes.data,1, sum))
      cprev=cor(guess, guesslast)
      message("correlation with pervious iteration=",cprev)
    }
    else{
      message("Converged");
      break
    }
  }
  if(return.all){
    return(list(result=guess, data=bayes.data, BF=Bayes.factors))
  }
  else{
    guess
  }
}

# row.names can not have duplicates
data = read.table(file = input_file_path, row.names="gene_id", sep = "\t", header = TRUE)
if(ncol(data) == 1) {
  # do not coerce dataframe dimension, see ?`[`
  data = data[order(data[1]),, drop=FALSE]
  data["gene_id"] = rownames(data)
  data["birra_ranking"] = 1:nrow(data)
} else {
  # run BIRRA
  birra_result = BIRRA(data, return.all=TRUE)
  data["birra_ranking"] = birra_result["result"]
  # add gene_id column and move it to the first column
  data["gene_id"] = rownames(data)
  cols = colnames(data)
  col_order = c("gene_id", cols[! cols %in% c("gene_id")])
  data = data[, col_order]
  # sort dataframe by birra result
  data = data[order(data["birra_ranking"]),]
}
# write result to output file
write.table(data, if (endsWith(output_file_path, ".gz")) gzfile(output_file_path) else output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
