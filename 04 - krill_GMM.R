rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Krill multimodality ---- #####
load(MSM115_KRILL_COMPLETE_DATA)
LF<-KRILL_FINAL$lf
SOI<-KRILL_FINAL$SOI
station_strata_table<-KRILL_FINAL$station_strata_table

df<-LF[LF$species %in% SOI,]
x<-range(df$length[!is.na(df$length)])
breaks<-seq(x[1],x[2],(x[2]-x[1])/60)

df$species<-factor(df$species,levels=c('Euphausia superba', 'Euphausia triacantha', 'Euphausia frigida', 'Thysanoessa macrura'), ordered=T)

p<-ggplot(df,aes(x=length,col=species,fill=species))
p<-p+scale_fill_manual(values=krill_colors)
p<-p+scale_colour_manual(values=krill_colors)
p<-p+geom_histogram(aes(y = after_stat(density)),breaks=breaks,col='grey',show.legend=F)
p<-p+geom_density(fill=NA,col='black')
p<-p+geom_rug(sides='t',na.rm=T,show.legend=F)
p<-p+labs(x='Body length [mm]',y='Density')
p<-p+facet_grid(species~.,scales='free_y')+MAINTHEME + theme_large_font()

df<-LF[LF$species =='Salpa thompsoni',]
df$state<-'solitary'
df$state[df$stage=='AGG']<-'colony'
df$state<-as.factor(df$state)
x<-range(df$length[!is.na(df$length)])
breaks<-seq(x[1],x[2],(x[2]-x[1])/100)

g<-ggplot(df,aes(x=length,fill=state),show.legend=F)
g<-g+scale_fill_manual(values=tunicate_stage_colors)
g<-g+scale_colour_manual(values=tunicate_stage_colors)
g<-g+geom_histogram(aes(y = after_stat(density)),breaks=breaks,show.legend=F,col='grey')
g<-g+geom_density(fill=NA,col='black')

g<-g+geom_rug(sides='t',na.rm=T,show.legend=F)
g<-g+labs(x='Body length [mm]',y='Density')
g<-g+ggh4x::facet_nested(state~species,scales='free_y')+MAINTHEME + theme_large_font()

m = matrix(c(1, 1, 2,2),ncol=2)

png(file.path(PUB_PACKAGE_DIR,'Figure 9.png'),6400,3200,res=300)
gridExtra::grid.arrange(grobs=list(p,g),layout_matrix=m,width=c(3,2))
graphics.off()

#### ----- Running Multimodality models ---- #####
out<-data.frame()
for (spp in SOI){
  length<-subset(LF,species==spp)
  
  for (st in unique(length$station_id)){
    
    l<-length$length[length$station_id==st]
    l<-l[!is.na(l)]
    if (length(l) < 5){
      next
    }
    
    mod<-mclust::densityMclust(log(l),modelNames=c('V'),plot=F)
    
    mod_table<-data.frame(model_type=mod$modelName,component=1:mod$G,mixing_prop=mod$parameters$pro,cumulative_prop=NA,mu=mod$parameters$mean, stdev=mod$parameters$variance$sigmasq)
    mod_table$cumulative_prop[1]<-mod_table$mixing_prop[1]
    if (nrow(mod_table)>1){
      for (r in 2:nrow(mod_table)){
        mod_table$cumulative_prop[r]<-sum(mod_table$mixing_prop[1:r])
      }
    }
    x<-diptest::dip.test(l)
    
    out<-rbind(out,data.frame(species=spp, station_id=st, n=mod$n, BIC=mod$bic, df=mod$df, loglik=mod$loglik,dip_test_p=x$p.value, n_components=mod$G,mod_table))
  }
}
normal_scale<-log_to_normal(out$mu, out$stdev)
out$mu<-normal_scale$normal_mean
out$stdev<-normal_scale$normal_sd

out$stratum<-sqldf::sqldf('select * from out as a left join station_strata_table as b on a.station_id=b.station_id')$stratum
out$x<-sqldf::sqldf('select * from out as a left join station_strata_table as b on a.station_id=b.station_id')$x
out$y<-sqldf::sqldf('select * from out as a left join station_strata_table as b on a.station_id=b.station_id')$y
out$dip_sig<-stars.pval(out$dip_test_p)
out$station<-paste0('Station ',formatC(out$station_id,width=2,flag='0'))
out$is_multimodal<-0
out$is_multimodal[out$n_components>1]<-1
multimodal<-out[out$n_components>1,]
table(out$model_type)

#### plotting ####
df<-multimodal
for (spp in unique(df$species)){
  idx<-df$species==spp
  
  df$ymin[idx]=min(df$mu[idx]-df$stdev[idx])
  df$ymax[idx]=max(df$mu[idx]+df$stdev[idx])
  
}

for (st in unique(df$station_id)){
  idx<-df$station_id==st
  for (cmp in unique(df$component[idx])){
    cmp_idx<-which(df$station_id==st & df$component==cmp)
    df$xmin[cmp_idx]=df$cumulative_prop[cmp_idx]-df$mixing_prop[cmp_idx]
    df$xmax[cmp_idx]=df$cumulative_prop[cmp_idx]
  }
}

df$component<-as.factor(df$component)

p<-ggplot(data=df,aes(x=cumulative_prop,y=mu, fill=component))
p<-p+geom_rect(data=df,mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)  )
p<-p+geom_segment(aes(x=xmin, y=mu, xend=xmax, yend=mu, group=component), colour="black")
p<-p+geom_segment(aes(x=xmin, y=mu+stdev, xend=xmax, yend=mu+stdev, group=component), colour="black",lty=2)
p<-p+geom_segment(aes(x=xmin, y=mu-stdev, xend=xmax, yend=mu-stdev, group=component), colour="black",lty=2)
p<-p+geom_text(data=df, x=0, y=Inf, hjust=0,vjust=1,mapping=aes(label=paste0('Dip test p: ', round(dip_test_p,2), dip_sig)))

p<-p+scale_fill_manual(values=component_colors[1:max(df$n_components)])
p<-p+labs(x='Cumulative proportion of component [%]',y='Body length [mm]')
p<-p+scale_x_continuous(breaks=c(0,.5,1),labels=c('0','50','100'))
p<-p+ggh4x::facet_nested(species~stratum+station, scale='free_y')

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_1.png'),5800,3500,res=300)
p+MAINTHEME+theme_large_font(14)
graphics.off()

spp<-'Euphausia superba'
p<-draw_components(spp, multimodal,LF,dip_test_size = 8) + MAINTHEME
p<-p+ggh4x::facet_nested(station+stratum~., scale='free_y')

g<-draw_components(spp, out,LF)
g<-g+ggh4x::facet_nested(station~stratum, scale='free_y')

n_x<-length(unique(p$data$stratum))
n_y<-length(unique(p$data$station_id))

if (all(c(n_x,n_y)>0)){
  png(file.path(PUB_PACKAGE_DIR,'Figure 10.png'),max(n_x,n_y)*800,max(n_x,n_y)*800,res=300)
  print(p+theme_large_font(14))
  graphics.off()
}

final_multimodal_table<-out[order(-out$is_multimodal,out$species,out$station),]
rownames(final_multimodal_table)<-1:nrow(final_multimodal_table)
openxlsx::write.xlsx(final_multimodal_table,MSM115_KRILL_MULTIMODALITY_FILE)

### spatial data
LL1<-sqldf::sqldf('select stratum, station, station_id,species,x,y,sum(n) as n, dip_test_p, dip_sig,n_components, is_multimodal from out group by stratum, station, station_id,species,x,y order by species, stratum, station_id')

### one column per species
transformed<-LL1[1,c(1:3,5,6)]
for (spp in unique(out$species)){
  transformed[[paste0('component_',spp)]]<-NA
  transformed[[paste0('is_multimodal_',spp)]]<-NA
  transformed[[paste0('dip_sig_',spp)]]<-NA
}
transformed<-transformed[-1,]
for (st in unique(LL1$station)){
  idx <- which( LL1$station==st )
  dat<-LL1[idx,]
  new_line<-transformed[1,]
  new_line$stratum<-unique(dat$stratum)
  new_line$station<-unique(dat$station)
  new_line$station_id<-unique(dat$station_id)
  new_line$x<-unique(dat$x)
  new_line$y<-unique(dat$y)
  
  for (r in 1:nrow(dat)){
    spp<-dat$species[r]
    new_line[[paste0('component_',spp)]]<-dat$n_components[r]
    new_line[[paste0('dip_sig_',spp)]]<-dat$dip_sig[r]
    new_line[[paste0('is_multimodal_',spp)]]<-as.numeric(dat$n_components[r]>1)
  }
  transformed<-rbind(transformed,new_line)
}
transformed$is_multimodal<-0
cols<-grep('is_multimodal.*',names(transformed))
transformed$is_multimodal[apply(transformed[,cols],1,sum)>0]<-1
transformed<-transformed[order(-transformed$is_multimodal,transformed$station),]

LL1<-sf::st_as_sf(LL1,coords=c('x','y'),crs=IBCSO)
LL2<-sf::st_as_sf(transformed,coords=c('x','y'),crs=IBCSO)

sf::st_write(LL1,MSM115_KRILL_MULTIMODALITY_GPKG, layer='summary',append=F, delete_dsn=T)
sf::st_write(LL1[LL1$dip_sig=='' & LL1$is_multimodal==1,],MSM115_KRILL_MULTIMODALITY_GPKG, layer='summary_dip_test',append=T)
sf::st_write(LL2,MSM115_KRILL_MULTIMODALITY_GPKG, layer='detailed',append=T)
sf::st_write(LL2[LL2$is_multimodal==1,],MSM115_KRILL_MULTIMODALITY_GPKG, layer='detailed_multimodality',append=T)
