### Plotting functions for FIV 16S study modified from code written by Dr. Zaid Abdo

### Species bar plot per sample using most abundant in sample split by treatment level
## This function requires an otu table (on any taxonomic level) in the form of a data frame
## It only plots the most abundant based on the cutoff given (the cutoff defaults to 0.1)
## title defaults to "Experiment", x label defaults to "Sample" but can be changed, and 
## y label defaults to "OTU" but can be changed. The function plots proportions and it is
## best to use the non-normalized data for this. 
## It splits the plot per treatment level as provided by Trt, which is required. Trt is a data frame 
## that includes 2 columns (the first is a treatment level and the second is the sample names or another factor

NC_plot_ftn = function(data,Trt, cutoff=0.01, llab="OTU", palette){
  lt = length(Trt[1,])
  nm = colnames(data)
  # reducing the plot df using the lowest treatment
  trt.ls = list()
  for(i in 1:lt) trt.ls[[i]] = levels(factor(Trt[,i]))
  {
    trt = p.df = c()
    for(i in 1:length(trt.ls[[1]])){
      for(j in 1:length(trt.ls[[2]])){
        trt1 = c(trt.ls[[1]][i],trt.ls[[2]][j])
        p.df = rbind(p.df, apply(data[Trt[,1]==trt.ls[[1]][i]&Trt[,2]==trt.ls[[2]][j],],2,sum)) 
        trt = rbind(trt,trt1)
      }
    }
    colnames(trt) = c("trt1","trt2")
    Trt = trt
    rs.vt = apply(p.df,1,sum)
    Trt = Trt[rs.vt > 0, ]
    p.df = p.df[rs.vt > 0, ]
    rs.vt = rs.vt[rs.vt>0]
    colnames(p.df) = nm
    p.df = p.df/rs.vt
    if(cutoff > 0){
      cut.b = apply(p.df,2,function(x)sum(x>cutoff))
      p.df = p.df[,cut.b > 0]
      p.rs = apply(p.df,1,sum)
      p.df = p.df/p.rs
    }
    nm.trt = colnames(Trt)
    hist.df = sp = w = c()
    nm = colnames(p.df)
    ln = length(p.df[,1])
    for(j in 1:length(p.df[1,])){
      sp = c(sp,rep(nm[j],ln))
      w = c(w,p.df[,j])
      hist.df = rbind(hist.df,Trt)
    }
    trt1 = factor(hist.df[,1])
    trt2 = factor(hist.df[,2])
    hist.df = data.frame(trt1,trt2,sp,w)
    hist.df = hist.df[order(hist.df$trt1,hist.df$trt2),]
    
    print(ggplot(data=hist.df,aes(x=trt2, fill = sp)) +
      geom_bar(aes(weight=w),color="white",linetype=1) +
      labs(color=llab, x = "", y="Proportion") +
      scale_fill_manual(name = llab, values = palette) + #, values = palette)         
      theme_minimal() +
      guides(fill = guide_legend(ncol=1)))
    
  }
}


cat_plot_ftn = function(data,Trt, cutoff=0.01, llab="OTU", palette){
  lt = length(Trt[1,])
  nm = colnames(data)
  # reducing the plot df using the lowest treatment
  trt.ls = list()
  for(i in 1:lt) trt.ls[[i]] = levels(factor(Trt[,i]))
  {
    trt = p.df = c()
    for(i in 1:length(trt.ls[[1]])){
      for(j in 1:length(trt.ls[[2]])){
        trt1 = c(trt.ls[[1]][i],trt.ls[[2]][j])
        p.df = rbind(p.df, apply(data[Trt[,1]==trt.ls[[1]][i]&Trt[,2]==trt.ls[[2]][j],],2,sum)) 
        trt = rbind(trt,trt1)
      }
    }
    colnames(trt) = c("trt1","trt2")
    Trt = trt
    rs.vt = apply(p.df,1,sum)
    Trt = Trt[rs.vt > 0, ]
    p.df = p.df[rs.vt > 0, ]
    rs.vt = rs.vt[rs.vt>0]
    colnames(p.df) = nm
    p.df = p.df/rs.vt
    if(cutoff > 0){
      cut.b = apply(p.df,2,function(x)sum(x>cutoff))
      p.df = p.df[,cut.b > 0]
      p.rs = apply(p.df,1,sum)
      p.df = p.df/p.rs
    }
    nm.trt = colnames(Trt)
    hist.df = sp = w = c()
    nm = colnames(p.df)
    ln = length(p.df[,1])
    for(j in 1:length(p.df[1,])){
      sp = c(sp,rep(nm[j],ln))
      w = c(w,p.df[,j])
      hist.df = rbind(hist.df,Trt)
    }
    trt1 = factor(hist.df[,1])
    trt2 = factor(hist.df[,2])
    hist.df = data.frame(trt1,trt2,sp,w)
    hist.df = hist.df[order(hist.df$trt1,hist.df$trt2),]
    
    print(ggplot(data=hist.df,aes(x=factor(trt2, levels = c("Cat1", "Cat2", "Cat3", "Cat4", "Cat5", "Cat6", #ctrl
                                                            "Cat7", "Cat8", "Cat9", "Cat10", "Cat11", "Cat12", #placebo
                                                            "Cat13", "Cat14", "Cat15", "Cat16", "Cat17", "Cat18") #cART
                                           ), fill = sp)) +
            geom_bar(aes(weight=w),color="white",linetype=1) +
            facet_grid(rows = vars(factor(trt1, 
                                          levels = c("Week -1", "Week 5", "Week 11", "Week 24")))) +
            labs(color = llab, x = "", y="Proportion") +
            scale_fill_manual(name = llab, values = palette) +
            theme_minimal() + 
            guides(fill = guide_legend(ncol=1)))
        
  }
}



