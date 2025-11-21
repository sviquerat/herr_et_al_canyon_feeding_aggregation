require(ggplot2)
require(RColorBrewer)

##### -------- Plotting settings -------- #####
MAINTHEME<-ggplot2::theme_bw()


##### --------------- Colours ------------ ######

high_colour='#ce0000'
mid_colour='#ececec'
low_colour='#0092ce'

gg_pca_colours<-c('#4a4a4a','#2e3149','#ce0000')

canyon_colors=c("blind" = '', "on_shelf" = '', "incising"='')
cluster_colors = c("5*" = "#00AFBB", "1" = "#E7B800", "2" = "#911eb4","3" = "#de1fba","4" = "#98FB98","6"="#FC4E07")
cluster_name_colors = c("aggregations" = "#00AFBB", "nondescript" = "#E7B800")
cluster_source_colors = c("a priori classification" = "#00AFBB", "a priori classification and a posteriori classification" = "#E7B800",'a posteriori classification' = "#FC4E07")
component_colors = c("#ed6a5a", "#4DA8DA", "#FFD66B","#98FB98")
krill_colors = c('Euphausia superba' = "#f73532", 'Euphausia frigida' = "#b197fe", 'Euphausia triacantha' = "#add662", "Thysanoessa macrura" = '#F6BD60')
obs_presence_colors = c("#80c4fe", "#FFD66B", '#98FB98') # periods of observer presence; used in currents analysis
tunicate_colors = c("Salpa thompsoni"='#84A59D')
tunicate_stage_colors = c('#84A59D','#289d9e')
pal_depth<-colorRampPalette(c('darkblue','lightblue'))
pal_cluster<-colorRampPalette(cluster_colors)
pal_krill<-colorRampPalette(krill_colors)

CTD_COLORMAP=list(temp1_C_filter1='red',
              temp2_C_filter1='red',
              salinity_filter1='yellow',
              oxygen_ml_l_filter1='#80c4fe',
              oxygen_saturation_perc_filter1='#e449ff',
              fluorescence_mg_m3_filter1='green',
              turbidity_NTU_filter1='orange'
)

### ggplot themes ####
theme_large_font<-function(base_size=12){
  out<-theme(axis.text=element_text(size=base_size),
             axis.title=element_text(size=base_size+round(1/4*base_size,0),face="bold"),
             legend.title=element_text(size=base_size),
             legend.text=element_text(size=base_size)) + 
    theme(strip.text = element_text(size = base_size))
  return(out)
}

theme_ctd_print<-function(){
  out <- 
    theme_linedraw() + 
    theme(
      axis.title.x=element_blank(), 
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank()
    ) +
    theme(
      legend.position = c(0.5, 0),
      legend.justification = c("center", "bottom"),
      legend.margin = ggplot2::margin(2, 2, 2, 2)
    )
  return(out)
}

theme_ctd<-function(...){
  out<-
    theme_ctd_print() +
    theme(
      legend.key = element_rect(fill =  DARK_BG),
      legend.background = element_rect(fill = DARK_BG, color = "transparent"), 
      legend.text = element_text(color=MAJOR_GRID)
    ) +
    theme(
      panel.background = element_rect(fill = DARK_BG),
      panel.grid.major = element_line(color = adjustcolor(MAJOR_GRID,MAJOR_GRID_ALPHA)),
      panel.grid.minor = element_line(color = adjustcolor(MINOR_GRID,MINOR_GRID_ALPHA)),
      panel.border=element_rect(color = adjustcolor(MAJOR_GRID,MAJOR_GRID_ALPHA)),
      plot.background = element_rect(fill = DARK_BG),
      axis.title.y=element_text(color = adjustcolor(MAJOR_GRID,MAJOR_GRID_ALPHA)), 
      axis.text.y=element_text(color = adjustcolor(MAJOR_GRID,MAJOR_GRID_ALPHA)), 
      axis.ticks.y=element_line(color = adjustcolor(MAJOR_GRID,MAJOR_GRID_ALPHA))
    )
  return (out)
}

theme_axis_print<-function(current_color){
  out<-
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          plot.background=element_blank()
    ) +
    theme(
      axis.text.x = element_text(color=current_color, face=3),
      axis.title.x = element_text(color=current_color, face=3,hjust = 0),
      axis.ticks.x = element_line(size = 1.5, color=current_color),
      axis.line = element_line(size = 1.5, color = current_color, linetype=1)
    )
  return(out)
}

theme_axis<-function(current_color){
  out<-
    theme_axis_print(current_color) +
    theme(
      plot.background = element_rect(fill = DARK_BG)
      )
  return (out)
}

