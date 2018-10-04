#This is the analysis file. The functions used in this file are cointained in synthetic_control_functions.R
#There are two model variants: 
# *_full - Full synthetic control model with all covariates (excluding user-specified covariates).
# *_time - Trend adjustment using the specified variable (e.g., non-respiratory hospitalization or population size) as the denominator.

#############################
#                           #
#    System Preparations    #
#                           #
#############################

source('synthetic_control_functions.R', local = TRUE)

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

packages <- c('parallel', 'splines', 'lubridate','loo', 'RcppRoll','pomp','lme4', 'rstanarm', 'ggplot2', 'reshape','dummies')
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
################################
#set up for STL+PCA
################################
##SECTION 1: CREATING SMOOTHED VERSIONS OF CONTROL TIME SERIES AND APPENDING THEM ONTO ORIGINAL DATAFRAME OF CONTROLS
#EXTRACT LONG TERM TREND WITH DIFFERENT LEVELS OF SMOOTHNESS USING STL
# Set a list of parameters for STL
stl.covars<-mapply(smooth_func,ds.list=ds,covar.list=covars_full, SIMPLIFY=FALSE) 
post.start.index<-which(time_points==post_period[1])

if (length(groups)>1){ 
  stl.data.setup<-mapply(stl_data_fun,covars=stl.covars, ds.sub=ds ,SIMPLIFY=FALSE)  #list of lists that has covariates for each regression for each strata
}else{
  stl.data.setup <- list(mapply(stl_data_fun,covars=stl.covars, ds.sub=ds ))
}

##SECTION 2: run first stage models
n_cores <- detectCores()-1
glm.results<- vector("list",  length=length(stl.data.setup)) #combine models into a list
cl1 <- makeCluster(n_cores)
clusterEvalQ(cl1, {library(lme4, quietly = TRUE)})
clusterExport(cl1, c('stl.data.setup',  'glm.fun', 'time_points', 'n_seasons','post.start.index'), environment())
for(i in 1:length(stl.data.setup)){
  print(i)
  glm.results[[i]]<-parLapply(cl=cl1 ,     stl.data.setup[[i]], fun=glm.fun )
}
stopCluster(cl1)
######################

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,       covars = covars_full), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=TRUE), groups)
data_pca<-mapply(FUN=pca_top_var,glm.results.in=glm.results, covars=stl.covars,ds.in=ds, SIMPLIFY=FALSE)
names(data_pca)<-groups
#Null model where we only include seasonal terms but no covariates
data_null <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_null, trend=FALSE), groups)
#Time trend model but without a denominator
data_time_no_offset <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=FALSE), groups)
data_allvars <- vector("list", length(data_full)) 
for(i in 1: length(data_full)){
  data_allvars[[i]]<-   cbind.data.frame(data_full[[i]], pc1=scale(data_pca[[i]][,ncol(data_pca[[i]])]), time.index=scale(1:nrow(data_pca[[i]])) )
}
  
  
###############################
#                             #
#        Main analysis        #
#                             #
###############################

impact_full <- setNames(lapply( data_full, doCausalImpact, intervention_date = intervention_date, adapt_delta=0.9999,n_iter=3000, var.select.on=TRUE, time_points = time_points,n_cores=n_cores), groups)
impact_time <- setNames(lapply( data_time, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points,n_cores=n_cores, trend = TRUE), groups)
impact_time_no_offset <- setNames(lapply( data_time_no_offset, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points,n_cores=n_cores,  trend = FALSE), groups)
impact_pca <- setNames(lapply( data_pca, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points,n_cores=n_cores), groups)
impact_allvars <- setNames(lapply( data_allvars, doCausalImpact, intervention_date = intervention_date,adapt_delta=0.9999,n_iter=3000, var.select.on=TRUE, time_points = time_points,n_cores=n_cores), groups)


#Test1
#groups2.test<-groups[1:2]
#data_full.test<-data_full[1:2]
#impact_full_test1 <- setNames(lapply( data_full.test, doCausalImpact, intervention_date = intervention_date, adapt_delta=0.9999,n_iter=3000, var.select.on=TRUE, time_points = time_points,n_cores=n_cores), groups2.test)
# #impact_full_test1 <- setNames(lapply( data_full.test, doCausalImpact, intervention_date = intervention_date, adapt_delta=0.9999,n_iter=3000, var.select.on=TRUE, time_points = time_points,n_cores=n_cores), groups2.test)
# n_iter=2000
# test1<-stan_glmer(form1, data = data.fit, family = poisson()  ,prior=hs ( df =1 , global_df =1 , global_scale = tau0 ) , QR=TRUE,   chains = 4, cores = n_cores, seed = 123 ,iter=n_iter)
# test2<-stan_glmer(form1, data = data.fit, family = poisson()  ,prior=hs() , QR=TRUE,   chains = 4, cores = n_cores, seed = 123 ,iter=n_iter)
# test3<-stan_glmer(form1, data = data.fit, family = poisson()  ,prior=hs() ,prior_intercept =normal(autoscale=FALSE),   QR=TRUE,   chains = 4, cores = n_cores, seed = 123 ,iter=n_iter)
# test4<-stan_glmer(form1, data = data.fit, family = poisson()  ,prior=hs() ,prior_intercept =normal(autoscale=FALSE),     chains = 4, cores = n_cores, seed = 123 ,iter=n_iter)

#test4<-stan_glmer(form1, data = data.fit, family = poisson()  ,prior=hs_plus() ,prior_intercept =normal(autoscale=FALSE),   QR=TRUE,   chains = 4, cores = n_cores, seed = 123 ,iter=n_iter)
#test5<-stan_glmer(form1, data = data.fit, family = poisson()  ,prior=hs_plus() ,prior_intercept =normal(autoscale=TRUE),   QR=TRUE,   chains = 4, cores = n_cores, seed = 123 ,iter=n_iter)

##########################################################################
##########################################################################
#Save the inclusion probabilities from each of the models
inclusion_prob_full <- setNames(lapply(impact_full, inclusionProb), groups)
inclusion_prob_allvars <- setNames(lapply(impact_allvars, inclusionProb), groups)
inclusion_prob_time <- setNames(lapply(impact_time, inclusionProb), groups)

#All model results combined
quantiles_full <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_full[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time[[group]], denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time_no_offset <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time_no_offset[[group]], denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)}), groups)
quantiles_pca <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_pca[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
quantiles_allvars <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_allvars[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)


#Model predicitons
pred_quantiles_full <- sapply(quantiles_full, getPred, simplify = 'array')
pred_quantiles_time <- sapply(quantiles_time, getPred, simplify = 'array')
pred_quantiles_time_no_offset <- sapply(quantiles_time_no_offset, getPred, simplify = 'array')
pred_quantiles_pca <- sapply(quantiles_pca, getPred, simplify = 'array')
pred_quantiles_allvars <- sapply(quantiles_allvars, getPred, simplify = 'array')

#Predictions, aggregated by year
ann_pred_quantiles_full <- sapply(quantiles_full, getAnnPred, simplify = FALSE)
ann_pred_quantiles_time <- sapply(quantiles_time, getAnnPred, simplify = FALSE)
ann_pred_quantiles_time_no_offset <- sapply(quantiles_time_no_offset, getAnnPred, simplify = FALSE)
ann_pred_quantiles_pca <- sapply(quantiles_pca, getAnnPred, simplify = FALSE)
ann_pred_quantiles_allvars <- sapply(quantiles_allvars, getAnnPred, simplify = FALSE)

#Pointwise RR and uncertainty for second stage meta analysis
log_rr_quantiles   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
dimnames(log_rr_quantiles)[[1]] <- time_points
log_rr_sd   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
log_rr_full_t_samples.prec<-sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
saveRDS(log_rr_quantiles, file=paste0(output_directory, country, "_log_rr_quantiles.rds"))
saveRDS(log_rr_sd, file=paste0(output_directory, country, "_log_rr_sd.rds"))
saveRDS(log_rr_full_t_samples.prec, file=paste0(output_directory, country, "_log_rr_full_t_samples.prec.rds"))

#Rolling rate ratios
rr_roll_full <- sapply(quantiles_full, FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array')
rr_roll_time <- sapply(quantiles_time, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')
rr_roll_time_no_offset <- sapply(quantiles_time_no_offset, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')
rr_roll_pca <- sapply(quantiles_pca, FUN = function(quantiles_pca) {quantiles_pca$roll_rr}, simplify = 'array')
rr_roll_allvars <- sapply(quantiles_allvars, FUN = function(quantiles_allvars) {quantiles_allvars$roll_rr}, simplify = 'array')


#Rate ratios for evaluation period.
rr_mean_full <- t(sapply(quantiles_full, getRR))
rr_mean_time <- t(sapply(quantiles_time, getRR))
rr_mean_time_no_offset <- t(sapply(quantiles_time_no_offset, getRR))
rr_mean_pca <- t(sapply(quantiles_pca, getRR))
rr_mean_allvars <- t(sapply(quantiles_allvars, getRR))


rr_mean_full_intervals <- data.frame('SC Estimate (95% CI)'     = makeInterval(rr_mean_full[, 2], rr_mean_full[, 3], rr_mean_full[, 1]), check.names = FALSE, row.names = groups)
rr_mean_time_intervals <- data.frame('Time trend Estimate (95% CI)' = makeInterval(rr_mean_time[, 2], rr_mean_time[, 3], rr_mean_time[, 1]), check.names = FALSE, row.names = groups)
rr_mean_time_no_offset_intervals <- data.frame('Time trend (no offset) Estimate (95% CI)' = makeInterval(rr_mean_time_no_offset[, 2], rr_mean_time_no_offset[, 3], rr_mean_time_no_offset[, 1]), check.names = FALSE, row.names = groups)
rr_mean_pca_intervals <- data.frame('STL+PCA Estimate (95% CI)'     = makeInterval(rr_mean_pca[, 2], rr_mean_pca[, 3], rr_mean_pca[, 1]), check.names = FALSE, row.names = groups)
rr_mean_allvars_intervals <- data.frame('All vars Estimate (95% CI)'     = makeInterval(rr_mean_allvars[, 2], rr_mean_allvars[, 3], rr_mean_allvars[, 1]), check.names = FALSE, row.names = groups)

colnames(rr_mean_time) <- paste('Time_trend', colnames(rr_mean_time))

#Combine RRs into 1 file for plotting
rr_mean_combo<- as.data.frame(rbind( cbind(rep(1, nrow(rr_mean_full)),groups,  seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_full),
                       cbind(rep(2, nrow(rr_mean_time)),groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)), rr_mean_time),
                       cbind(rep(3, nrow(rr_mean_time_no_offset)),groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)), rr_mean_time_no_offset),
                    cbind(rep(4, nrow(rr_mean_pca)), groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_pca),
                    cbind(rep(5, nrow(rr_mean_allvars)), groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_allvars)) )
            
        names(rr_mean_combo)<-c('Model', 'groups', 'group.index','lcl','mean.rr','ucl')
       # if(crossval){
          #point.weights2<-stacking_weights.all.m
        # }else{
          point.weights2<-as.data.frame(matrix(rep(1,nrow(rr_mean_combo)), ncol=1))
           names(point.weights2)<-'value'
        # }
        rr_mean_combo$point.weights<-point.weights2$value
        rr_mean_combo$group.index<-as.numeric(as.character(rr_mean_combo$group.index))
        rr_mean_combo$mean.rr<-as.numeric(as.character(rr_mean_combo$mean.rr))
        rr_mean_combo$lcl<-as.numeric(as.character(rr_mean_combo$lcl))
        rr_mean_combo$ucl<-as.numeric(as.character(rr_mean_combo$ucl))
        rr_mean_combo$group.index[rr_mean_combo$Model==2]<-rr_mean_combo$group.index[rr_mean_combo$Model==2]+0.15
        rr_mean_combo$group.index[rr_mean_combo$Model==3]<-rr_mean_combo$group.index[rr_mean_combo$Model==3]+0.3
        rr_mean_combo$group.index[rr_mean_combo$Model==4]<-rr_mean_combo$group.index[rr_mean_combo$Model==4]+0.45
        rr_mean_combo$group.index[rr_mean_combo$Model==5]<-rr_mean_combo$group.index[rr_mean_combo$Model==5]+0.55
        
        rr_mean_combo$Model<-as.character(rr_mean_combo$Model)
        rr_mean_combo$Model[rr_mean_combo$Model=='1']<-"Synthetic Controls"
        rr_mean_combo$Model[rr_mean_combo$Model=='2']<-"Time trend"
        rr_mean_combo$Model[rr_mean_combo$Model=='3']<-"Time trend (No offset)"
        rr_mean_combo$Model[rr_mean_combo$Model=='4']<-"STL+PCA"
        rr_mean_combo$Model[rr_mean_combo$Model=='5']<-"Combined"
        cbPalette <- c("#1b9e77", "#d95f02", "#7570b3",'#e7298a','#66a61e')
        rr_mean_combo$est.index<-as.factor(1:nrow(rr_mean_combo))
        #Fix order for axis
        rr_mean_combo$Model<-as.factor(rr_mean_combo$Model)
        rr_mean_combo$Model = factor(rr_mean_combo$Model,levels(rr_mean_combo$Model)[c(3,4,5,2,1)])
        #print(levels(rr_mean_combo$Model))

cumsum_prevented <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_full, simplify = 'array')
cumsum_prevented_pca <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_pca, simplify = 'array')
cumsum_prevented_time <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_time, simplify = 'array')
cumsum_prevented_allvars <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_allvars, simplify = 'array')


##Save output from allvars analysis
save.allvars.est<-list(post_period,outcome_plot, time_points,ann_pred_quantiles_allvars, pred_quantiles_allvars,rr_roll_allvars,rr_mean_allvars,rr_mean_allvars_intervals,cumsum_prevented_allvars)
names(save.allvars.est)<-c('post_period','outcome_plot','time_points', 'ann_pred_quantiles_allvars', 'pred_quantiles_allvars','rr_roll_allvars','rr_mean_allvars','rr_mean_allvars_intervals','cumsum_prevented_allvars')
saveRDS(save.allvars.est, file=paste0(output_directory, country, "allvars estimates.rds"))
log_rr_quantiles_allvars   <- sapply(quantiles_allvars,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
dimnames(log_rr_quantiles_allvars)[[1]] <- time_points
log_rr_full_t_samples.allvars.prec<-sapply(quantiles_allvars,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec.post}, simplify = 'array')
saveRDS(log_rr_quantiles_allvars, file=paste0(output_directory, country, "_log_rr_quantiles_allvars.rds"))
saveRDS(log_rr_full_t_samples.allvars.prec, file=paste0(output_directory, country, "_log_rr_full_t_samples.allvars.prec.rds"))
#

################################
#                              #
#     Sensitivity Analyses     #
#                              #
################################

  sensitivity_table_intervals <- NA

