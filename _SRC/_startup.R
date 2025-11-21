#### Startup
SRC_DIR<-file.path(getwd(),'_SRC')
print(SRC_DIR)

set.seed(666)

.required_packages<-c('sqldf','terra','sf','ggplot2','dplyr', 'oce', 'signal','vegan','tidyterra','factoextra','e1071','mgcv','NbClust','ggh4x','mclust','RColorBrewer')

package_check<-function(){
  print('Checking package requirements...')
  
  for (pkg in .required_packages){
    if (! (pkg %in% installed.packages()) ) {
      print(paste0('Need to install package: ',pkg))
      install.packages(pkg,dep=T)
    }
  }
  rm(pkg)
  
  print('Done!')
}

srcs<-list.files(SRC_DIR,'*.R', full.names = T)
for (src in srcs){
  if (basename(src) == '_startup.R'){next}
  source(src)
}
