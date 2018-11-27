#This is the analysis file. The functions used in this file are cointained in synthetic_control_functions.R
#There are two model variants: 
# *_full - Full synthetic control model with all covariates (excluding user-specified covariates).
# *_time - Trend adjustment using the specified variable (e.g., non-respiratory hospitalization or population size) as the denominator.

#############################
#                           #
#    System Preparations    #
#                           #
#############################

source('./main analysis components/synthetic_control_functions.R', local = FALSE)
#source('synthetic_control_functions.R', local = FALSE)

#############################
#Automatically set working directory to desktop
#setwd('~/synthetic-control-master/main analysis components')  #directory where .Rmd file is saved
#Set working directory: default to desktop--different path for windows vs Mac
if(.Platform$OS.type == "windows") {
  desktop<-file.path(Sys.getenv("USERPROFILE"),"Desktop")
  desktop<-gsub(pattern='\\',replacement='/', desktop, fixed=TRUE)
} else {
  desktop<- "~/Desktop"
}
auto.wd<-file.path(paste0(desktop,'/synthetic-control-poisson-master/main analysis components/'))
#

packages <- c('parallel', 'splines', 'lubridate','loo','MASS', 'RcppRoll','pomp','lme4', 'ggplot2', 'reshape','dummies')
packageHandler(packages, update_packages, install_packages)
sapply(packages, library, quietly = TRUE, character.only = TRUE)

#Detect if pogit package installed; if not download archive (no longer on cran)
if("BayesLogit" %in% rownames(installed.packages())==FALSE){
  if(.Platform$OS.type == "windows") {
  #url_BayesLogit<- "https://mran.microsoft.com/snapshot/2017-02-04/src/contrib/BayesLogit_0.6.tar.gz"
  install_github("jwindle/BayesLogit")
  }else{
    url_BayesLogit<- "https://github.com/weinbergerlab/synthetic-control-poisson/blob/master/packages/BayesLogit_0.6_mac.tgz?raw=true"
  }
  pkgFile_BayesLogit <- "BayesLogit.tar.gz"
  download.file(url = url_BayesLogit, destfile = pkgFile_BayesLogit)
  install.packages(url_BayesLogit, type="source", repos=NULL)
}
if("pogit" %in% rownames(installed.packages())==FALSE){
  url_pogit <- "https://cran.r-project.org/src/contrib/Archive/pogit/pogit_1.1.0.tar.gz"
  pkgFile_pogit <- "pogit_1.1.0.tar.gz"
  download.file(url = url_pogit, destfile = pkgFile_pogit)
  install.packages(pkgs=pkgFile_pogit, type="source", repos=NULL)
  install.packages('logistf')
}
library(pogit)

#Detects number of available cores on computers. Used for parallel processing to speed up analysis.
n_cores <- detectCores()
set.seed(1)

###################################################
#                                                 #
# Directory setup and initialization of constants #
#                                                 #
###################################################

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
groups <- as.character(unique(unlist(prelog_data[, group_name], use.names = FALSE)))
if (exists('exclude_group')) {groups <- groups[!(groups %in% exclude_group)]}

###############################################
#                                             #
# Data and covariate preparation for analysis #
#                                             #
###############################################

#Make sure we are in right format
prelog_data[,date_name]<-as.Date(as.character(prelog_data[,date_name]), tryFormats=c("%m/%d/%Y",'%Y-%m-%d' ))

#test<-split(prelog_data, factor(prelog_data[,group_name]))
#outcome.na<-sapply(test, function(x) sum(is.na(x[,outcome_name])))
prelog_data[, date_name] <- formatDate(prelog_data[, date_name])
prelog_data <- setNames(lapply(groups, FUN = splitGroup, ungrouped_data = prelog_data, group_name = group_name, date_name = date_name, start_date = start_date, end_date = end_date, no_filter = c(group_name, date_name, outcome_name, denom_name)), groups)
#if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

#Log-transform all variables, adding 0.5 to counts of 0.
ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name,outcome_name)), groups)
time_points <- unique(ds[[1]][, date_name])

#Monthly dummies
if(n_seasons==4){dt<-quarter(as.Date(time_points))}
if(n_seasons==12){dt<-month(as.Date(time_points))}
if(n_seasons==3){
  dt.m<-month(as.Date(time_points))
  dt<-dt.m
  dt[dt.m %in% c(1,2,3,4)]<-1
  dt[dt.m %in% c(5,6,7,8)]<-2
  dt[dt.m %in% c(9,10,11,12)]<-3
    }
season.dummies<-dummy(dt)
season.dummies<-as.data.frame(season.dummies)
names(season.dummies)<-paste0('s', 1:n_seasons)
season.dummies<-season.dummies[,-n_seasons]

ds <- lapply(ds, function(ds) {
	if (!(denom_name %in% colnames(ds))) {
		ds[denom_name] <- 0
	}
	return(ds)
})

sparse_groups <- sapply(ds, function(ds) {
	return(ncol(ds[!(colnames(ds) %in% c(date_name, group_name, denom_name, outcome_name, exclude_covar))]) == 0)
})
ds <- ds[!sparse_groups]
groups <- groups[!sparse_groups]

#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars_full <- setNames(lapply(ds, makeCovars), groups)
covars_full <- lapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), groups)
covars_null <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies)))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {data[, outcome_name]})
outcome_plot=outcome
offset<- sapply(ds, FUN=function(data) exp(data[, denom_name]) )  #offset term on original scale; 1 column per age group

post.start.index<-which(time_points==post_period[1])

#Combine the outcome, covariates, and time point information.
#data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,       covars = covars_full), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=TRUE), groups)

rr.its1<-lapply(data_time,its_func)
rr.t<-sapply(rr.its1, `[[`, "rr.q.t", simplify='array')
rr.end<-t(sapply(rr.its1, `[[`, "rr.q.post", simplify='array')) 

matplot(rr.t[,,10], bty='l', type='l', lty=c(2,1,2), col='gray')
abline(h=1)

write.csv(rr.end, paste(output_directory, country, 'rr_classic_its.csv', sep = ''))
