##### Distance sampling related functions

require(ggplot2)
require(RColorBrewer)
require(Distance)


bulk_ds<-function(ds_data, cov.mat, region_table=NA, obs_table=NA, sample_table=NA, convert_units=1,start=1,...){
  #### bulk distance sampling
  #### when adjustments are chosen, the model chosen within the adjustment series is based on AIC. 
  #### Some of these adjusted models may be not strictly monotonic
  out<-data.frame()
  perc_complete=0
  
  for (i in start:nrow(cov.mat)){
    step<-i %% floor(nrow(cov.mat)/10)
    if ( is.na(step) || step == 0){
      perc_complete=perc_complete+5
      print(paste0(perc_complete,' % complete'))
    }
    print(i)
    cov<-cov.mat[i,]
    if (is.na(cov$formula)){
      form<-'~1'
    }else{
      form<-paste0('~1+',cov$formula) 
    }
    
    modelname<-paste0(cov$key,'_',formatC(i,digits=2,width=2,flag='0'))
    adj<-cov$adjustment
    if (is.na(adj)){
      adj<-NULL
    }
    n_samples<-nrow(ds_data[ds_data$distance<=cov$truncation,])
    cov$formula<-form
    
    simple_run<-any(is.na(c(region_table,obs_table,sample_table)))
    m<-tryCatch(
      {
        if (simple_run){
          m<-Distance::ds(data=ds_data, key=cov$key, formula=as.formula(form), truncation=cov$truncation, 
                          adjustment=adj, convert_units=convert_units,...)
          
        }else{
          m<-Distance::ds(data=ds_data, region_table = region_table, sample_table=sample_table, obs_table=obs_table, 
                          key=cov$key, formula=as.formula(form), truncation=cov$truncation, adjustment=adj,
                          convert_units=convert_units,...)
        }
      },
      warning=function(w){
        print(w)
        m<-m
      },
      error=function(e){
        print(e)
        m<-NA
        
      },
      finally={
      }
    )
    
    
    if(is.na(m)[1]){
      cov$formula<-form
      n_samples<-nrow(m$ddf$data)
      
      new_data<-data.frame(i=i,modelname=modelname,cov,used_adj=NA, order=NA,
                           AIC=NA,CvM_p = NA, CvM_W = NA, pa_hat=NA, pa_hat_se=NA, 
                           N_sigs=n_samples,is_mono=NA,N_hat=NA)
    }else{
      p_hat<-m$ddf$fitted
      if (cov$key=='unif'){
        a<-data.frame(pa_hat=1,pa_hat_se=0)
      }else{
        a<-Distance::summarize_ds_models(m, output='plain', delta_only=T)[,5:6]
      }
      names(a)<-c('pa_hat','pa_hat_se')
      gof<-Distance::gof_ds(m,plot=F)$dsgof$CvM
      mono<-check.mono(m$ddf, plot=F)
      
      adjustment<-m$ddf$ds$aux$ddfobj$adjustment$series
      order_adjustment<-paste0('(',paste(m$ddf$ds$aux$ddfobj$adjustment$order,collapse=', '),')')
      if (is.null(adjustment)){
        adjustment<-order_adjustment<-NA
      }
      N_hat<-m$dht$clusters$N[4,2]
      if (is.null(N_hat)){
        N_hat<-NA
      }
      new_data<-data.frame(i=i,modelname=modelname,cov,used_adj=adjustment, order=order_adjustment,
                           AIC=AIC(m)$AIC,CvM_p = gof$p, CvM_W = gof$W, a, 
                           N_sigs=n_samples,is_mono=mono,
                           N_hat=N_hat)
    }
    out<-rbind(out,new_data)
  }
  row.names(out)<-1:nrow(out)
  return(out)
}

create_ds_settings<-function(truncations,adjustments=c(NA,'cos','herm','poly'),keys=c('hn','hr','unif'),covars=NA){
  
  cov.mat<-expand.grid(keys, adjustments, covars, truncations)
  names(cov.mat)<-c('key','adjustment','formula', 'truncation')
  cov.mat$formula<-as.character(cov.mat$formula)
  cov.mat$key<-as.character(cov.mat$key)
  cov.mat$adjustment<-as.character(cov.mat$adjustment)
  
  ### adjustments are not considered when using covars, so we can skip these
  idx<-which(!is.na(cov.mat$formula) & !is.na(cov.mat$adjustment))
  if (length(idx) > 0){
    cov.mat<-cov.mat[-idx,]
  }
  
  ### uniform cannot include covars, so we can skip these
  idx<-which(!is.na(cov.mat$formula) & cov.mat$key=='unif')
  if (length(idx) > 0){
    cov.mat<-cov.mat[-idx,]
  }
  
  cov.mat<-cov.mat[order(cov.mat$key,cov.mat$formula,cov.mat$truncation,cov.mat$adjustment,na.last=F),]
  rownames(cov.mat)<-1:nrow(cov.mat)
  return(cov.mat)
}

check_covar_levels<-function(ds_data, truncations,covars,distance_col='distance'){
  covar_df<-data.frame()
  for (trunc in truncations){
    df<-ds_data[ds_data$distance <= trunc,]
    for (cov in covars[!is.na(covars)]){
      x<-as.data.frame(table(df[,cov]))
      
      x<-cbind(rep(trunc,nrow(x)),rep(cov,nrow(x)),x)
      names(x)<-c('truncation','covar','covar_level','count')
      covar_df<-rbind(covar_df,x)
    }
  }
  return(covar_df)
}

##### plots

det.fct.plot<-function(ds.model,COLOURED=F,esw.col='black',lgd_position='topright',plot=T,ndist=250,...){
  options(warn=-1)
  o<-ds.model$ddf$data
  o$esw<-predict(ds.model$ddf,esw=T,newdata=o)$fitted
  esw_mean<-mean(unique(o$esw))
  
  
  left <- ds.model$ddf$meta.data$left
  width <- ds.model$ddf$meta.data$width
  
  covar<-ds.model$ddf$ds$aux$ddfobj$scale$formula
  if (covar=='~1'){
    covar<-NULL
    df<-sqldf::sqldf('select esw from o group by esw')
    
  }else{
    covar<-gsub('~1 \\+ ','',covar)
    df<-sqldf::sqldf(paste0('select esw, ', covar,' from o group by esw, ',covar))
  }
  
  plot.new()
  h_all<-h<-mrds::add_df_covar_line(ds.model,data=df,lty=0,ndist=ndist)
  graphics.off()
  
  #h_all<-h<-mrds::add_df_covar_line(ds.model,data=o,lty=0,ndist=ndist)
  h_mean<-apply(h,2,'mean')
  
  
  out<-list(o=o,esw_mean=esw_mean,h_all=h_all, h_mean = h_mean,left=left, right=width, covar=covar)
  
  
  if (plot){
    xx <- seq(left, width, length.out=250)
    x_steps<-pretty(seq(0,m$ddf$meta.data$width),10)
    plot(ds.model,showpoints=F,type='n',pl.col='#F1F1F1',border='darkgrey',ann=F,xaxt='n')
    title(ylab='Detection probability',xlab='Distance [m]')
    axis(1,at=x_steps,labels=format(x_steps,big.mark=','))
    lines(h_mean~xx,col=esw.col,lty=1,lwd=1.5)
    
    abline(v=esw_mean,col=esw.col,lwd=1.5)
    if(!is.null(covar)){
      txt<-NULL
      lvls<-sort(levels(as.factor(ds.model$ddf$data[[covar]])))
      cols<-hcl.colors(length(lvls),palette='Set 2')
      for (i in 1:length(lvls)){
        lvl<-lvls[i]
        col<-cols[i]
        df<-data.frame(covar=lvl)
        names(df)<-covar
        mrds::add_df_covar_line(ds.model,data=df,col=col,lwd=1)
        esw<-unique(o$esw[o[[covar]]==lvl])
        N<-sum(ds.model$ddf$data[[covar]]==lvl)
        txt<-c(txt,paste0(covar,': ',lvl,' (N: ',N,', esw: ',format(round(esw,0),big.mark=","),' m)'))
      }
      legend(lgd_position,legend=txt,lty=2,col=cols,lwd=2)
    }
    text(esw_mean,0,paste(format(round(esw_mean,0),big.mark=","),' m',sep=' '),srt=90,pos=4,col=esw.col)
    options(warn=0)
  }
  invisible(out)
}

get_ds_model_params<-function(ds.model){
  ddfobj<-ds.model$ddf$ds$aux$ddfobj
  type=ddfobj$type
  
  adjustment<-ddfobj$adjustment
  adjustment_series<-adjustment$series
  adjustment_order<-adjustment$order
  adjustment_estimate<-adjustment$parameters
  
  if (is.null(adjustment)){
    adjustment_series<-NA
    adjustment_order<-NA
    adjustment_estimate<-NA
  }
  
  left <- ds.model$ddf$meta.data$left
  right <- ds.model$ddf$meta.data$width
  AIC<-summary(ds.model$ddf)$aic
  summry<-Distance::summarize_ds_models(ds.model)
  key_function<-summry$`Key function`
  p_hat<-summry[,5]
  p_hat_se<-summry[,6]
  
  gof<-Distance::gof_ds(ds.model,plot=F)$dsgof
  
  form<-strsplit(as.character(ds.model$ddf$model),',')[[2]][2]
  form<-gsub(')','',strsplit(form,'~')[[1]][2])
  covars<-trimws(strsplit(form,'\\+')[[1]])
  
  if (length(covars)==1){
    covars<-NULL
  }else{
    covars<-covars[2:length(covars)]
  }
  
  results<-data.frame(left_truncation=left, right_truncation=right,
                      key_type=type, key_verbatim=key_function,series=adjustment_series,order=adjustment_order,
                      series_estimate=adjustment_estimate, p_hat = p_hat, p_hat_se = p_hat_se,
                      AIC=AIC)
  rownames(results)<-NULL
  out<-list(data=ds.model$ddf$data,covars=covars,model_result=results,gof=gof)
  
  return(out)
}

gg_DS_DetFct<-function(ds.model,
                       COLOURED=F,
                       ndist = 250,
                       limit_to_density=F,
                       add_average=F,
                       palette.name='Set1',
                       avg.lwd=1,
                       avg.lty=1,
                       cov.lwd=.5,
                       cov.lty=2,
                       col.det_fct='black',
                       col.esw='purple',
                       size.esw=2,
                       thousands_sep=NULL,
                       n_hist_bins=7,
                       auto_adjust_breaks=T,
                       n_digits_axis=0,
                       n_digits_esw=0,
                       ...){
  
  hist_scale<-function(x,left,right,n_breaks,max_x,p_hat){
    bp<-hist(x,plot=F,breaks=seq(left,right,length.out=n_breaks))
    binwidth=bp$breaks[2]
    area_hist<-sum(bp$counts*binwidth)
    area_total<-max_x
    n<-length(bp$counts)
    
    conv<-p_hat/area_hist*area_total
    
    hist_scaled<-data.frame(mids=bp$mids,x0 = bp$breaks[1:(length(bp$breaks)-1)],scaled_counts=bp$counts*conv, width=binwidth/2)
    return(hist_scaled)
  }
  
  n_hist_breaks=n_hist_bins+1
  out<-det.fct.plot(ds.model,plot=F,ndist=ndist)
  xx <- seq(out$left, out$right, length.out=length(out$h_mean))
  max_x<-max(xx)
  
  summry<-Distance::summarize_ds_models(ds.model)
  p_hat<-summry[,5]
  
  key<-ds.model$ddf$name.message
  det_fct_name<-key
  plot_data<-data.frame()
  if (!is.null(out$covar)){
    det_fct_name<-'Average'
    lvls<-sort(levels(as.factor(ds.model$ddf$data[[out$covar]])))
    for (r in 1:nrow(out$h_all)){
      lvl=lvls[r]
      plot_data<-rbind(plot_data,data.frame(y=out$h_all[r,],x=xx,lvl=lvl))
    }
  }
  plot_data<-rbind(plot_data,data.frame(x=xx,y=out$h_mean,lvl=det_fct_name))
  
  lvls<-unique(plot_data$lvl)
  n_vars<-length(lvls)
  
  palette=c(col.det_fct,RColorBrewer::brewer.pal(n_vars-1,palette.name))
  
  plot_data$lvl<-factor(plot_data$lvl, levels = c(det_fct_name,lvls[lvls!=det_fct_name]))
  
  # The scaling is such that the area under the detection curve is the same as the area of the histogram bars.
  if (auto_adjust_breaks){
    print('auto adjustment will be implemented later')
    optim<-data.frame()
    for (i in 2:20){
      X<-hist_scale(out$o$distance,left = out$left,right=out$right,n_breaks=i,max_x=max_x,p_hat=p_hat)
      areas<-distances<-c()
      for (m in 1:length(X$mids)){
        midpoint<-X$mids[m]
        x0<-midpoint-X$width[m]
        x1<-midpoint+X$width[m]
        l<-x0-xx
        r<-x1-xx
        idx1<-which(abs(l)==min(abs(l)))
        idx2<-which(abs(r)==min(abs(r)))
        y1<-out$h_mean[idx1][1]
        y2<-out$h_mean[idx2][1]
        hyp_area<-2*X$width[m]*max(y2,y1)
        bin_area<-2*X$width[m]*X$scaled_counts[m]
        areas<-c(areas,hyp_area-bin_area)
        distances<-c(distances,abs(plot_data$y-X$scaled_counts[m]))
      }
      optim<-rbind(optim,data.frame(breaks=i,total_dist=sum(distances),mean_dist=mean(distances),area=sum(areas))) 
    }
    n_optim<-optim$breaks[optim$W==min(optim$W)]
    print(n_optim)
  }
  
  bp<-hist(out$o$distance,plot=T,breaks=seq(out$left,out$right,length.out=n_hist_breaks))
  
  binwidth=bp$breaks[2]
  
  area_hist<-sum(bp$counts*binwidth)
  area_total<-max(xx)
  n<-length(bp$counts)
  
  conv<-p_hat/area_hist*area_total
  
  hist_scaled<-data.frame(mids=bp$mids,x0 = bp$breaks[1:(length(bp$breaks)-1)],scaled_counts=bp$counts*conv, width=binwidth/2)
  n_bins<-nrow(hist_scaled)
  
  max_y<-ceiling(max(hist_scaled$scaled_counts)*10)/10
  hist_palette<-colorRampPalette(c("darkgrey", "lightgrey"))
  
  labels<-seq(0,max_x,binwidth)
  labels=round(labels,n_digits_axis)
  esw_mean<-round(out$esw_mean,n_digits_esw)
  if (!is.null(thousands_sep)){
    labels=formatC(labels, format="f", big.mark=thousands_sep,drop0trailing=T)
    esw_mean=formatC(esw_mean, format="f", big.mark=thousands_sep,drop0trailing=T)
  }

    ### the actual plot ###
  p<-ggplot()
  p<-p+xlab('Distance [m]')+ylab('Detection probability')
  p<-p+scale_x_continuous(breaks=seq(0,max_x,binwidth), labels=labels)
  p<-p+ylim(0,max_y)#+xlim(0,max_x)
  
  p<-p+geom_bar(data=hist_scaled,stat='identity',inherit.aes=F,mapping=aes(x = mids, y = scaled_counts), width=binwidth, position=position_dodge(),
                fill=hist_palette(n_bins))
  p<-p+geom_vline(xintercept=out$esw_mean,col=col.esw,lty=1,alpha=.5,lwd=1)
  p<-p+annotate(geom='text', x=out$esw_mean, y=.1, label=paste0(esw_mean,' m'), angle=90, col=col.esw,vjust=1, size=size.esw)
  
  if (!is.null(out$covar)){
    p<-p+geom_line(data=subset(plot_data,lvl==det_fct_name),inherit.aes=F,mapping=aes(x=x,y=y,col=lvl),lty=avg.lty)
    p<-p+geom_line(data=subset(plot_data,lvl!=det_fct_name),inherit.aes=F,mapping=aes(x=x,y=y,col=lvl),lty=cov.lty,lwd=cov.lwd)
    p<-p+scale_colour_manual(name=out$covar,values=palette)
  }else{
    p<-p+geom_line(data=subset(plot_data,lvl==det_fct_name),inherit.aes=F,mapping=aes(x=x,y=y),col=col.det_fct,lty=avg.lty)
  }
  return(p)
}

gg_DS_GOF<-function(ds.model, add_points=T,add_poly=F,add_lines=F){
  params<-get_ds_model_params(ds.model)
  data=params$data
  
  df<-data.frame(object=data$object,distance=data$distance,size=data$size,edf1=params$gof$edf[,1],edf2=params$gof$edf[,2],cdf=params$gof$cdf)
  
  dist<-data.frame()
  for (r in 1:nrow(df)){
    dd<-data.frame(x0=df$edf2[r],y0=df$edf2[r],x1=df$edf2[r],y1=df$cdf[r])
    dd$dist<-df$cdf[r]-df$edf2[r]
    dist<-rbind(dist,dd)
  }
  
  covars<-params$covars
  if (!is.null(covars)){
    for (cov in covars){
      df[[cov]]  <- data[[cov]]
    }
  }
  
  p<-ggplot(data=df,aes(x=edf2,y=cdf))
  
  if (add_poly){
    polys<-data.frame(x=0,y=0)
    for (r in 1:nrow(dist)){
      x<-dist$x1[r]
      y<-dist$y1[r]
      polys<-rbind(polys,data.frame(x=x,y=y))
    }
    p<-p+ggnewscale::new_scale_colour() +geom_polygon(data=polys,inherit.aes=F,mapping=aes(x=x,y=y), fill='red',alpha=.3)
  }
  p<-p+geom_abline(intercept=0,slope=1,lwd=1,col='green')
  
  if (add_lines){
    p<-p+ggnewscale::new_scale_colour() +geom_segment(data=dist,inherit.aes=F,mapping=aes(x=x0,xend=x1,y=y0,yend=y1,col=dist),lwd=1)
    p<-p+scale_colour_continuous(low='red',high='purple',name = "Distance from center line", breaks=seq(min(df$size),max(df$size),1))
  }
  
  if (add_points){
    p<-p+ggnewscale::new_scale_colour()+geom_point()
  }  
  
  p<-p+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))
  p<-p+xlab('Empirical CDF')+ylab('Fitted CDF')
  return(p)
}
