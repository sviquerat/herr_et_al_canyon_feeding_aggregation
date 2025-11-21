#### PCA related functions
require(ggplot2)
require(gridExtra)
# require(factoextra)
# require(GGally)

# add new line including all info neeed from a pca run
add_line<-function(pca_data,subset,step,text){
  pca_result<-prcomp(pca_data,scale=T,center=T)
  s<-data.frame(summary(pca_result)$importance[,1:3])
  s$var<-rownames(s)
  rownames(s)<-NULL
  print(paste0('current Cumulative Proportion of Variance on PC2:',100*s[3,3], ' %'))
  return(data.frame(step=step, subset=subset,text=text,s))
}

drop1_pca<-function(pca_data,axis=1:3,drop=T){
  columns<-1:ncol(pca_data)
  out<-data.frame()
  pca_result<-prcomp(pca_data,scale=T,center=T)
  s<-data.frame(summary(pca_result)$importance[,axis])
  s$var<-rownames(s)
  s$col_number<-0
  rownames(s)<-NULL
  out<-rbind(out,s)
  cp_o<-out[3,]
  if (drop){
    for (i in columns){
      if (ncol(pca_data)<=max(axis)){
        print('Not enough data to run the analysis')
        return (NULL)
      }
      pca_result<-prcomp(pca_data[,-i],scale=T,center=T)
      s<-data.frame(summary(pca_result)$importance[,axis])
      s$var<-rownames(s)
      rownames(s)<-NULL
      s$col_number<-i
      out<-rbind(out,s)
    }
  }
  out$exclusion<-'no exclusion'
  out$exclusion[out$col_number>0]<-names(pca_data)[out$col_number[out$col_number>0]]
  
  cp<-out[out$var=='Cumulative Proportion',]
  idx<-which(cp$PC2==max(cp$PC2))
  best<-cp[idx,]
  print(paste0('Cum Prop was: ',100*cp_o$PC2 ,'.'))
  print(paste0('Excluding column ',best$exclusion, ' increases Cum Prop on PC2 to:',100*best$PC2, ' %'))
  
  return(out[out$exclusion==best$exclusion,])
}

brute_force_dropper_pca<-function(pca_data){
  
  reduction<-drop1_pca(pca_data, drop=F)
  
  while(T){
    col<-drop1_pca(pca_data, drop=T)
    if (is.null(col)){break}
    reduction<-rbind(reduction,col)
    pca_data<-pca_data[,-col$col_number[]]
  }
  
  rownames(reduction)<-1:nrow(reduction)
  reduction$step<-ceiling(as.numeric(rownames(reduction)) / 3)
  
  reduction$pretty_excl<-paste0(formatC(reduction$step,flag='0',width=2), ' - ',reduction$exclusion)
  reduction$pretty_excl<-factor(reduction$pretty_excl,levels=sort(unique(reduction$pretty_exc)))
  
  env_col_names<-oce_col_names[env_col_names %in% names(pca_dat)]
  oce_col_names<-oce_col_names[oce_col_names %in% names(pca_dat)]
  
  df<-reshape2::melt(reduction, id.vars=c(4,7,8), measure.vars=c(1,2,3), variable.name = 'Axis')
  p<-ggplot(data=subset(df),aes(x=step,y=value,col=pretty_excl,group=Axis))+geom_line(aes(x=step,y=value,group=Axis))
  p<-p+geom_point(cex=3)+coord_flip() + scale_x_reverse() + scale_x_continuous(breaks=1:max(reduction$step), labels=1:max(reduction$step))
  p<-p+MAINTHEME+facet_grid(Axis~var,scales='free_x')+labs(col='order of exclusion')
  
  return(list(data=reduction,p=p))
}

# GGPLOTS

gg_pca_contrib<-function(pca_result,axes=c(1,2), gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),...){
  fv<-factoextra::fviz_pca_var(pca_result, axes=axes,col.var="contrib", gradient.cols = gradient.cols,repel = TRUE)+ggCircle(.5,lty=2)
  return(fv)
}

gg_pca_summary<-function(pca_dat, axes=c(1,2), main='PCA result contribution', add_cor=F, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),add_letters=T,...){
  pca_result<-prcomp(pca_dat,scale=T,center=T)
  p<-gg_pca_contrib(pca_result,axes=axes, gradient.cols = gradient.cols,...)+labs(title=NULL)
  l<-factoextra::fviz_screeplot(pca_result, addlabels = TRUE)+labs(title=NULL)
  g<-ggHeatmap(cor(pca_dat))+labs(title=NULL)
  if (add_letters){
    p<-p+labs(tag = "A")
    l<-l+labs(tag = "B")
    g<-g+labs(tag = "C")
  }
  if (add_cor){
    layout<-matrix(c(1,2,3,3,3,3),byrow=T,ncol=2)
    gridExtra::grid.arrange(p,l,g,layout_matrix=layout)
  }else{
    layout<-matrix(c(1,1,1,1,2,2),byrow=T,ncol=2)
    gridExtra::grid.arrange(p,l,layout_matrix=layout)
  }
  return(pca_result)
}

gg_pca_contrib_all<-function(pca_dat, axes=c(1,2), main='PCA result contribution', gradient.cols = c('#2166AC','#B2182B'),add_letters=T,grid_max_cols=2,...){
  pca_result<-prcomp(pca_dat,scale=T,center=T)
  comb<-combn(axes,2)
  P<-list()
  for (col in 1:ncol(comb)){
    ax<-comb[,col]
    p<-gg_pca_contrib(pca_result,axes=ax, gradient.cols = gradient.cols,...)+labs(title=NULL)
    if (add_letters){
      p<-p+labs(tag = LETTERS[col])
    }
    P<-c(P,list(p))
  }
  do.call("grid.arrange", c(P, ncol=min(c(grid_max_cols,length(axes)))))
  return(pca_result)
}

get_axis_arrows<-function(rotation,axis=c(1,2)){
  arrows<-data.frame()
  for (r in 1:nrow(rotation)){
    loadings<-x$rotation[r,]
    lbl<-row.names(x$rotation)[r]
    arrows<-rbind(arrows,data.frame(var=lbl,PC_x=axis[1], PC_y=axis[2],x=0,y=0,xend=loadings[axis[1]],yend=loadings[axis[2]]))
  }
  return(arrows)
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

ggCircle<-function(radius=1,...){
  p<-list()
  for (r in radius){
    dat<-circleFun(c(0,0),r*2,npoints=200)
    p<-append(p,geom_path(data=dat,aes(x,y),...))
  }
  
  return(p)
}

ggAxisLoads<-function(pca_dat, which='both',axis=NULL,colours=NULL,text_pos_contrib_col='black',text_neg_contrib_col='black'){
  if (is.null(axis)){
    axis<-1:ncol(pca_dat)
  }
  
  if (is.null(colours)){
    colours<-cm.colors(nrow(pca_dat))
  }
  
  
  loads<-data.frame()
  pca_res<-prcomp(pca_dat,scale=T,center=T)
  
  for (ax in 1:ncol(pca_res$rotation)){
    labels<-names(pca_res$rotation[,ax])
    values<-as.numeric(pca_res$rotation[,ax])
    values_round<-ceiling(round(values,3)*1000)/1000
    loads<-rbind(loads,data.frame(var=labels,values=values,values_round=values_round,axis=ax))
  }
  
  n_vars<-nrow(pca_res$rotation)
  var_colours<-colours
  dark_colours<-colorspace::darken(colours,.3)
  
  loadings<-loads
  loadings<-loadings[order(loadings$var,loadings$axis),]
  loadings$axis<-as.factor(loadings$axis)
  loadings$colour<-text_pos_contrib_col
  loadings$colour[loadings$values_round<0]<-text_neg_contrib_col
  loadings$values<-abs(loadings$values)
  
  max_axis<-max(axis)
  ax_loads<-loadings[as.numeric(loadings$axis)<=max_axis,]
  
  if (which %in%c('vars','both')){
    
    p<-ggplot(data=ax_loads,aes(x=var,y=values,fill=axis,label=values_round,col=var))
    p<-p+geom_bar(stat='identity',position='fill')
    p<-p+scale_fill_manual(values=var_colours)
    p<-p+scale_colour_manual(values=dark_colours)
    p<-p + geom_text(size = 3, col=ax_loads$colour,position = position_fill(vjust = 0.5))
    p<-p+guides(col='none')
    p<-p+MAINTHEME+labs(x='',y='') + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    
    if (which == 'vars'){
      return(p)
    }
    
  }
  
  if (which %in%c('axis','both')){
    g<-ggplot(data=ax_loads,aes(x=axis,y=values,fill=var,label=values_round,col=var))
    g<-g+geom_bar(stat='identity',position='fill')
    #g<-g+scale_fill_brewer(palette = "Set3")
    g<-g+scale_fill_manual(values=var_colours)
    g<-g+scale_colour_manual(values=dark_colours)
    g<-g+guides(col='none')
    
    g<-g + geom_text(size = 3,col=ax_loads$colour, position = position_fill(vjust = 0.5))
    g<-g+MAINTHEME+labs(x='',y='',fill='') + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    
    if (which == 'axis'){
      return(g)
    }
    
  }
  grobs<-gridExtra::arrangeGrob(grobs=list(p+labs(tag='A'),g+labs(tag='B')), ncol=2)
  return(grobs)
}

ggTotalLoading<-function(pca_dat,axis=NULL,sorted=F,reverse_sort=F,colours=NULL){
  if(is.null(axis)){
    axis<-1:ncol(pca_dat)
  }
  
  
  pca_res<-prcomp(pca_dat,scale=T,center=T)
  df<-as.data.frame(pca_res$rotation)
  df$var<-rownames(df)
  rownames(df)<-NULL
  df[,axis]<-abs(df[,axis])
  df$total_contribution<-apply(df[,axis],FUN = sum,MARGIN = 1)
  
  if (is.null(colours)){
    colours<-cm.colors(length(unique(df$var)))
  }
  var_colours<-colours
  
  if (sorted){
    idx<-order(df$total_contribution)
    if (reverse_sort){idx<-rev(idx)}
    df<-df[idx,]
    df$var<-factor(df$var,levels=df$var)
  }
  
  p<-ggplot(data=df,aes(x=var,y=total_contribution,fill=var))
  p<-p+geom_col() +  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  p<-p+labs(x='',y=paste0('Total contribution across axis ',min(axis),' - ', max(axis)))
  p<-p+scale_fill_manual(values=var_colours)+guides(fill='none')
  
  return(p)
  
}
