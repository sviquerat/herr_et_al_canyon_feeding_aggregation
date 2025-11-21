#### various ggplot functions

require(ggplot2)
require(gridExtra)
require(factoextra)
require(GGally)

gg_pairs<-function(num_dat,groups,columns=NULL,group_name='ID',data_scale=T,cor.method='pearson',force_single_plot=T,n_vars_per_plot=6,show_outliers=T,alpha=.3,colours=NULL,theme=NULL,point.size=10,diag.size=4,...){
  
  dat<-num_dat
  
  
  if (is.null(columns)){
    columns<-1:ncol(dat)
  }
  
  dat<-dat[,columns]
  
  if (data_scale){
    dat<-data.frame(scale(dat,center = T, scale=T))
  }
  
  dat<-cbind(groups,dat)
  names(dat)[1]<-group_name
  
  if (force_single_plot){
    p<-GGally::ggpairs(data=dat,
                       columns=1:ncol(dat),
                       mapping=ggplot2::aes(col=dat[[group_name]],alpha=alpha),
                       lower = list(continuous = wrap("points", size=point.size), combo = GGally::wrap("dot_no_facet",size=diag.size)),
                       upper=list(continuous = GGally::wrap("cor",method=cor.method,size=diag.size), combo = GGally::wrap("box_no_facet", outliers=show_outliers), discrete = "count", na = "na"))  
    if (!is.null(colours)){
      p<-p + scale_fill_manual(values=colours) + scale_colour_manual(values=colours)
    }
    if (!is.null(theme)){
      p<-p+theme
    }
    
    return(p)
  }
  
  out<-list()
  n_vars<-n_vars_per_plot
  i=1
  for (col in seq(1,ncol(dat),n_vars+1)){
    local_columns<-c(1,(col+1):(col+n_vars))
    local_columns<-local_columns[local_columns<=ncol(dat)]
    local_dat<-as.data.frame(dat[,c(local_columns)])
    p<-GGally::ggpairs(data=local_dat,
                       columns=1:ncol(local_dat),
                       mapping=ggplot2::aes(col=local_dat[[group_name]],alpha=alpha),
                       lower = list(continuous = "points", combo = "dot_no_facet"),
                       upper=list(continuous = GGally::wrap("cor",method=cor.method), combo = GGally::wrap("box_no_facet", outliers=show_outliers), discrete = "count", na = "na"))
    
    if (!is.null(colours)){
      p<-p + scale_fill_manual(values=colours) + scale_colour_manual(values=colours)
    }
    
    if (!is.null(theme)){
      p<-p+theme
    }
    
    out[[i]]<-p
    i<-i+1
  }
  return(out)
}

ggHeatmap<-function(cormat,legend_name='Heatmap',lower_limit=min(cormat,na.rm=T),upper_limit=max(cormat,na.rm=T),sig_digits=4,
                    mid_point=lower_limit+((upper_limit-lower_limit)/2),removeTri=T,force_symmetry=T,low_col='blue',mid_col='white',high_col='red'){
  require(reshape2)
  upper_limit=upper_limit
  lower_limit=lower_limit
  mid_point=mid_point
  if (force_symmetry){
    lower_limit<- -upper_limit
    mid_point<-0
  }
  
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
    return(cormat)
  }
  
  #cormat <- reorder_mat(cormat)
  if (removeTri){
    cormat <- get_upper_tri(cormat)
  }
  melted<-reshape2::melt(cormat,na.rm=T)
  melted$value<-round(melted$value,sig_digits)
  p<-ggplot(data = melted, aes(Var2, Var1, fill = value))
  p<-p+geom_tile(color = "white")
  p<-p+scale_fill_gradient2(low = low_col, high = high_col, mid = mid_col, midpoint = mid_point, 
                            limit = c(lower_limit,upper_limit), space = "Lab", name=legend_name)
  p<-p+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
  p<-p+coord_fixed()
  p<-p+geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
  p<-p+MAINTHEME + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  p<-p+xlab('')+ylab('')
  p<-p+guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, angle=45))
  return(p)
}