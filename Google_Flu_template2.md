# Google Flu: CDC: ILI.log regression:: template2
bdanalytics  

**  **    
**Date: (Mon) Jun 22, 2015**    

# Introduction:  

Data: 
Source: 
    Training:   https://courses.edx.org/asset-v1:MITx+15.071x_2a+2T2015+type@asset+block/FluTrain.csv  
    New:        https://courses.edx.org/asset-v1:MITx+15.071x_2a+2T2015+type@asset+block/FluTest.csv  
Time period: 



# Synopsis:

Based on analysis utilizing <> techniques, <conclusion heading>:  

Regression results:
First run:
    <glb_sel_mdl_id>: 
        OOB_RMSE=<0.4f>; new_RMSE=<0.4f>; <feat1>=<imp>; <feat2>=<imp>

Classification results:
First run:
    <glb_sel_mdl_id>: Leaderboard: <accuracy>
        newobs_tbl=[0=, 1=]; submit_filename=
        OOB_conf_mtrx=[YN=, NY=]=; max.Accuracy.OOB=; opt.prob.threshold.OOB=
            <feat1>=<imp>; <feat1>=<imp>; <feat1>=<imp>; 
            <txt.feat1>=<imp>; <txt.feat1>=<imp>; <txt.feat1>=<imp>; 

### Prediction Accuracy Enhancement Options:
- import.data chunk:
    - which obs should be in fit vs. OOB (currently dirty.0 vs .1 is split 50%)
    
- inspect.data chunk:
    - For date variables
        - Appropriate factors ?
        - Different / More last* features ?
        
- scrub.data chunk:        
- transform.data chunk:
    - derive features from multiple features
    
- manage.missing.data chunk:
    - Not fill missing vars
    - Fill missing numerics with a different algorithm
    - Fill missing chars with data based on clusters 
    
- extract.features chunk:
    - Text variables: move to date extraction chunk ???
        - Mine acronyms
        - Mine places

- Review set_global_options chunk after features are finalized

### ![](<filename>.png)

## Potential next steps include:
- Organization:
    - Categorize by chunk
    - Priority criteria:
        0. Ease of change
        1. Impacts report
        2. Cleans innards
        3. Bug report
        
- all chunks:
    - at chunk-end rm(!glb_<var>)
    
- manage.missing.data chunk:
    - cleaner way to manage re-splitting of training vs. new entity

- extract.features chunk:
    - Add n-grams for glb_txt_vars
        - "RTextTools", "tau", "RWeka", and "textcat" packages
    - Convert user-specified mutate code to config specs
    
- fit.models chunk:
    - Prediction accuracy scatter graph:
    -   Add tiles (raw vs. PCA)
    -   Use shiny for drop-down of "important" features
    -   Use plot.ly for interactive plots ?
    
    - Change .fit suffix of model metrics to .mdl if it's data independent (e.g. AIC, Adj.R.Squared - is it truly data independent ?, etc.)
    - move model_type parameter to myfit_mdl before indep_vars_vctr (keep all model_* together)
    - create a custom model for rpart that has minbucket as a tuning parameter
    - varImp for randomForest crashes in caret version:6.0.41 -> submit bug report

- Probability handling for multinomials vs. desired binomial outcome
-   ROCR currently supports only evaluation of binary classification tasks (version 1.0.7)
-   extensions toward multiclass classification are scheduled for the next release

- Skip trControl.method="cv" for dummy classifier ?
- Add custom model to caret for a dummy (baseline) classifier (binomial & multinomial) that generates proba/outcomes which mimics the freq distribution of glb_rsp_var values; Right now glb_dmy_glm_mdl always generates most frequent outcome in training data
- glm_dmy_mdl should use the same method as glm_sel_mdl until custom dummy classifer is implemented

- fit.all.training chunk:
    - myplot_prediction_classification: displays 'x' instead of '+' when there are no prediction errors 
- Compare glb_sel_mdl vs. glb_fin_mdl:
    - varImp
    - Prediction differences (shd be minimal ?)

- Move glb_analytics_diag_plots to mydsutils.R: (+) Easier to debug (-) Too many glb vars used
- Add print(ggplot.petrinet(glb_analytics_pn) + coord_flip()) at the end of every major chunk
- Parameterize glb_analytics_pn
- Move glb_impute_missing_data to mydsutils.R: (-) Too many glb vars used; glb_<>_df reassigned
- Replicate myfit_mdl_classification features in myfit_mdl_regression
- Do non-glm methods handle interaction terms ?
- f-score computation for classifiers should be summation across outcomes (not just the desired one ?)
- Add accuracy computation to glb_dmy_mdl in predict.data.new chunk
- Why does splitting fit.data.training.all chunk into separate chunks add an overhead of ~30 secs ? It's not rbind b/c other chunks have lower elapsed time. Is it the number of plots ?
- Incorporate code chunks in print_sessionInfo
- Test against 
    - projects in github.com/bdanalytics
    - lectures in jhu-datascience track

# Analysis: 

```r
rm(list=ls())
set.seed(12345)
options(stringsAsFactors=FALSE)
source("~/Dropbox/datascience/R/myscript.R")
source("~/Dropbox/datascience/R/mydsutils.R")
```

```
## Loading required package: caret
## Loading required package: lattice
## Loading required package: ggplot2
```

```r
source("~/Dropbox/datascience/R/myplot.R")
source("~/Dropbox/datascience/R/mypetrinet.R")
source("~/Dropbox/datascience/R/myplclust.R")
# Gather all package requirements here
suppressPackageStartupMessages(require(doMC))
registerDoMC(4) # max(length(glb_txt_vars), glb_n_cv_folds) + 1
#packageVersion("snow")
#require(sos); findFn("cosine", maxPages=2, sortby="MaxScore")

# Analysis control global variables
glb_trnng_url <- "https://courses.edx.org/asset-v1:MITx+15.071x_2a+2T2015+type@asset+block/FluTrain.csv"
glb_newdt_url <- "https://courses.edx.org/asset-v1:MITx+15.071x_2a+2T2015+type@asset+block/FluTest.csv"
glb_out_pfx <- "template2_"
glb_save_envir <- FALSE # or TRUE

glb_is_separate_newobs_dataset <- TRUE    # or TRUE
    glb_split_entity_newobs_datasets <- TRUE   # or FALSE
    glb_split_newdata_method <- "sample"          # "condition" or "sample" or "copy"
    glb_split_newdata_condition <- NULL # or "is.na(<var>)"; "<var> <condition_operator> <value>"
    glb_split_newdata_size_ratio <- 0.3               # > 0 & < 1
    glb_split_sample.seed <- 123               # or any integer

glb_max_fitobs <- NULL # or any integer                         
glb_is_regression <- TRUE; glb_is_classification <- !glb_is_regression; 
    glb_is_binomial <- NULL # or TRUE or FALSE

glb_rsp_var_raw <- "ILI"

# for classification, the response variable has to be a factor
glb_rsp_var <- "ILI.log"

# if the response factor is based on numbers/logicals e.g (0/1 OR TRUE/FALSE vs. "A"/"B"), 
#   or contains spaces (e.g. "Not in Labor Force")
#   caret predict(..., type="prob") crashes
glb_map_rsp_raw_to_var <- function(raw) {
    return(log(raw))
#     ret_vals <- rep_len(NA, length(raw)); ret_vals[!is.na(raw)] <- ifelse(raw[!is.na(raw)] == 1, "Y", "N"); return(relevel(as.factor(ret_vals), ref="N"))
#     #as.factor(paste0("B", raw))
#     #as.factor(gsub(" ", "\\.", raw))    
}
glb_map_rsp_raw_to_var(c(5.660867, 6.339272, 6.815222, 7.388359, 7.618892))
```

```
## [1] 1.733577 1.846764 1.919159 1.999906 2.030631
```

```r
glb_map_rsp_var_to_raw <- function(var) {
    return(exp(var))
#     as.numeric(var) - 1
#     #as.numeric(var)
#     #gsub("\\.", " ", levels(var)[as.numeric(var)])
#     c("<=50K", " >50K")[as.numeric(var)]
#     #c(FALSE, TRUE)[as.numeric(var)]
}
glb_map_rsp_var_to_raw(glb_map_rsp_raw_to_var(c(5.660867, 6.339272, 6.815222, 7.388359, 7.618892)))
```

```
## [1] 5.660867 6.339272 6.815222 7.388359 7.618892
```

```r
if ((glb_rsp_var != glb_rsp_var_raw) & is.null(glb_map_rsp_raw_to_var))
    stop("glb_map_rsp_raw_to_var function expected")
glb_rsp_var_out <- paste0(glb_rsp_var, ".predict.") # model_id is appended later

# List info gathered for various columns
# <col_name>:   <description>; <notes>
# "Week" - The range of dates represented by this observation, in year/month/day format.
# 
# "ILI" - This column lists the percentage of ILI-related physician visits for the corresponding week.
# 
# "Queries" - This column lists the fraction of queries that are ILI-related for the corresponding week, adjusted to be between 0 and 1 (higher values correspond to more ILI-related search queries).

# If multiple vars are parts of id, consider concatenating them to create one id var
# If glb_id_var == NULL, ".rownames <- row.names()" is the default
glb_id_var <- c("Week")
glb_category_vars <- NULL # or c("<var1>", "<var2>")
glb_drop_vars <- c(NULL) # or c("<col_name>")

glb_map_vars <- NULL # or c("<var1>", "<var2>")
glb_map_urls <- list();
# glb_map_urls[["<var1>"]] <- "<var1.url>"

glb_assign_pairs_lst <- NULL; 
# glb_assign_pairs_lst[["<var1>"]] <- list(from=c(NA),
#                                            to=c("NA.my"))
glb_assign_vars <- names(glb_assign_pairs_lst)

# Derived features
glb_derive_lst <- NULL;
glb_derive_lst[["Week.bgn"]] <- list(
    mapfn=function(Week) { return(substr(Week, 1, 10)) }
    , args=c("Week"))
glb_derive_lst[["Week.end"]] <- list(
    mapfn=function(Week) { return(substr(Week, 14, 23)) }
    , args=c("Week"))

require(zoo)
```

```
## Loading required package: zoo
## 
## Attaching package: 'zoo'
## 
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```r
# # If glb_allobs_df is not sorted in the desired manner
# glb_derive_lst[["ILI.2.lag"]] <- list(
#     mapfn=function(Week) { return(coredata(lag(zoo(orderBy(~Week, glb_allobs_df)$ILI), -2, na.pad=TRUE))) }
#     , args=c("Week"))
glb_derive_lst[["ILI.2.lag"]] <- list(
    mapfn=function(ILI) { return(coredata(lag(zoo(ILI), -2, na.pad=TRUE))) }
    , args=c("ILI"))
glb_derive_lst[["ILI.2.lag.log"]] <- list(
    mapfn=function(ILI.2.lag) { return(log(ILI.2.lag)) }
    , args=c("ILI.2.lag"))

#     mapfn=function(PTS, oppPTS) { return(PTS - oppPTS) }
#     , args=c("PTS", "oppPTS"))

# Add logs of numerics that are not distributed normally ->  do automatically ???

#     mapfn=function(raw) { tfr_raw <- as.character(cut(raw, 5)); 
#                           tfr_raw[is.na(tfr_raw)] <- "NA.my";
#                           return(as.factor(tfr_raw)) }

# glb_derive_lst[["<txt_var>.niso8859.log"]] <- list(
#     mapfn=function(<txt_var>) { match_lst <- gregexpr("&#[[:digit:]]{3};", <txt_var>)
#                         match_num_vctr <- unlist(lapply(match_lst, 
#                                                         function(elem) length(elem)))
#                         return(log(1 + match_num_vctr)) }
#     , args=c("<txt_var>"))

#     mapfn=function(raw) { mod_raw <- raw;
#         mod_raw <- gsub("&#[[:digit:]]{3};", " ", mod_raw);
#         # Modifications for this exercise only
#         mod_raw <- gsub("\\bgoodIn ", "good In", mod_raw);
#                           return(mod_raw)

#         # Create user-specified pattern vectors 
# #sum(mycount_pattern_occ("Metropolitan Diary:", glb_allobs_df$Abstract) > 0)
#         if (txt_var %in% c("Snippet", "Abstract")) {
#             txt_X_df[, paste0(txt_var_pfx, ".P.metropolitan.diary.colon")] <-
#                 as.integer(0 + mycount_pattern_occ("Metropolitan Diary:", 
#                                                    glb_allobs_df[, txt_var]))
#summary(glb_allobs_df[ ,grep("P.on.this.day", names(glb_allobs_df), value=TRUE)])

# args_lst <- NULL; for (arg in glb_derive_lst[["Week.bgn"]]$args) args_lst[[arg]] <- glb_allobs_df[, arg]; do.call(mapfn, args_lst)

# glb_derive_lst[["<var1>"]] <- glb_derive_lst[["<var2>"]]
glb_derive_vars <- names(glb_derive_lst)

glb_date_vars <- c("Week.bgn", "Week.end")
glb_date_fmts <- list(); glb_date_fmts[["Week.bgn"]] <- glb_date_fmts[["Week.end"]] <- "%Y-%m-%d"; 
glb_date_tzs <- list();  glb_date_tzs[["Week.bgn"]] <- glb_date_tzs[["Week.end"]] <- "America/New_York"
#grep("America/New", OlsonNames(), value=TRUE)

glb_txt_vars <- NULL # or c("<txt_var1>", "<txt_var2>")   
#Sys.setlocale("LC_ALL", "C") # For english

glb_append_stop_words <- list()
# Remember to use unstemmed words
#orderBy(~ -cor.y.abs, subset(glb_feats_df, grepl("[HSA]\\.T\\.", id) & !is.na(cor.high.X)))
#dsp_obs(Headline.contains="polit")
#subset(glb_allobs_df, H.T.compani > 0)[, c("UniqueID", "Headline", "H.T.compani")]
# glb_append_stop_words[["<txt_var1>"]] <- c(NULL
# #                             ,"<word1>" # <reason1>
#                             )
#subset(glb_allobs_df, S.T.newyorktim > 0)[, c("UniqueID", "Snippet", "S.T.newyorktim")]
#glb_txt_lst[["Snippet"]][which(glb_allobs_df$UniqueID %in% c(8394, 8317, 8339, 8350, 8307))]

glb_important_terms <- list()
# Remember to use stemmed terms 

glb_sprs_thresholds <- NULL # or c(0.988, 0.970, 0.970) # Generates 29, 22, 22 terms
# Properties:
#   numrows(glb_feats_df) << numrows(glb_fitobs_df)
#   Select terms that appear in at least 0.2 * O(FP/FN(glb_OOBobs_df))
#       numrows(glb_OOBobs_df) = 1.1 * numrows(glb_newobs_df)
names(glb_sprs_thresholds) <- glb_txt_vars

# User-specified exclusions  
glb_exclude_vars_as_features <- c("Week", "ILI.2.lag") 
if (glb_rsp_var_raw != glb_rsp_var)
    glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
                                            glb_rsp_var_raw)

# List feats that shd be excluded due to known causation by prediction variable
glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
                                      c(NULL)) # or c("<col_name>")

glb_impute_na_data <- TRUE # or TRUE
glb_mice_complete.seed <- 144 # or any integer

glb_cluster <- FALSE # or TRUE

glb_interaction_only_features <- NULL # or ???

glb_models_lst <- list(); glb_models_df <- data.frame()
# Regression
if (glb_is_regression)
    glb_models_method_vctr <- c("lm", "glm", "bayesglm", "rpart", "rf") else
# Classification
    if (glb_is_binomial)
        glb_models_method_vctr <- c("glm", "bayesglm", "rpart", "rf") else  
        glb_models_method_vctr <- c("rpart", "rf")

# Baseline prediction model feature(s)
glb_Baseline_mdl_var <- NULL # or c("<col_name>")

glb_model_metric_terms <- NULL # or matrix(c(
#                               0,1,2,3,4,
#                               2,0,1,2,3,
#                               4,2,0,1,2,
#                               6,4,2,0,1,
#                               8,6,4,2,0
#                           ), byrow=TRUE, nrow=5)
glb_model_metric <- NULL # or "<metric_name>"
glb_model_metric_maximize <- NULL # or FALSE (TRUE is not the default for both classification & regression) 
glb_model_metric_smmry <- NULL # or function(data, lev=NULL, model=NULL) {
#     confusion_mtrx <- t(as.matrix(confusionMatrix(data$pred, data$obs)))
#     #print(confusion_mtrx)
#     #print(confusion_mtrx * glb_model_metric_terms)
#     metric <- sum(confusion_mtrx * glb_model_metric_terms) / nrow(data)
#     names(metric) <- glb_model_metric
#     return(metric)
# }

glb_tune_models_df <- 
   rbind(
    #data.frame(parameter="cp", min=0.00005, max=0.00005, by=0.000005),
                            #seq(from=0.01,  to=0.01, by=0.01)
    #data.frame(parameter="mtry",  min=080, max=100, by=10),
    #data.frame(parameter="mtry",  min=08, max=10, by=1),    
    data.frame(parameter="dummy", min=2, max=4, by=1)
        ) 
# or NULL
glb_n_cv_folds <- 3 # or NULL

glb_clf_proba_threshold <- NULL # 0.5

# Model selection criteria
if (glb_is_regression)
    glb_model_evl_criteria <- c("min.RMSE.OOB", "max.R.sq.OOB", "max.Adj.R.sq.fit")
if (glb_is_classification) {
    if (glb_is_binomial)
        glb_model_evl_criteria <- 
            c("max.Accuracy.OOB", "max.auc.OOB", "max.Kappa.OOB", "min.aic.fit") else
        glb_model_evl_criteria <- c("max.Accuracy.OOB", "max.Kappa.OOB")
}

glb_sel_mdl_id <- NULL # or "<model_id_prefix>.<model_method>"
glb_fin_mdl_id <- glb_sel_mdl_id # or "Final"

# Depict process
glb_analytics_pn <- petrinet(name="glb_analytics_pn",
                        trans_df=data.frame(id=1:6,
    name=c("data.training.all","data.new",
           "model.selected","model.final",
           "data.training.all.prediction","data.new.prediction"),
    x=c(   -5,-5,-15,-25,-25,-35),
    y=c(   -5, 5,  0,  0, -5,  5)
                        ),
                        places_df=data.frame(id=1:4,
    name=c("bgn","fit.data.training.all","predict.data.new","end"),
    x=c(   -0,   -20,                    -30,               -40),
    y=c(    0,     0,                      0,                 0),
    M0=c(   3,     0,                      0,                 0)
                        ),
                        arcs_df=data.frame(
    begin=c("bgn","bgn","bgn",        
            "data.training.all","model.selected","fit.data.training.all",
            "fit.data.training.all","model.final",    
            "data.new","predict.data.new",
            "data.training.all.prediction","data.new.prediction"),
    end  =c("data.training.all","data.new","model.selected",
            "fit.data.training.all","fit.data.training.all","model.final",
            "data.training.all.prediction","predict.data.new",
            "predict.data.new","data.new.prediction",
            "end","end")
                        ))
#print(ggplot.petrinet(glb_analytics_pn))
print(ggplot.petrinet(glb_analytics_pn) + coord_flip())
```

```
## Loading required package: grid
```

![](Google_Flu_template2_files/figure-html/set_global_options-1.png) 

```r
glb_analytics_avl_objs <- NULL

glb_chunks_df <- myadd_chunk(NULL, "import.data")
```

```
##         label step_major step_minor   bgn end elapsed
## 1 import.data          1          0 6.863  NA      NA
```

## Step `1.0: import data`
#### chunk option: eval=<r condition>

```r
#glb_chunks_df <- myadd_chunk(NULL, "import.data")

glb_trnobs_df <- myimport_data(url=glb_trnng_url, comment="glb_trnobs_df", 
                                force_header=TRUE)
```

```
## [1] "Reading file ./data/FluTrain.csv..."
## [1] "dimensions of data in ./data/FluTrain.csv: 417 rows x 3 cols"
##                      Week      ILI   Queries
## 1 2004-01-04 - 2004-01-10 2.418331 0.2377158
## 2 2004-01-11 - 2004-01-17 1.809056 0.2204515
## 3 2004-01-18 - 2004-01-24 1.712024 0.2257636
## 4 2004-01-25 - 2004-01-31 1.542495 0.2377158
## 5 2004-02-01 - 2004-02-07 1.437868 0.2244356
## 6 2004-02-08 - 2004-02-14 1.324274 0.2071713
##                        Week      ILI    Queries
## 15  2004-04-11 - 2004-04-17 0.836130 0.07569721
## 63  2005-03-13 - 2005-03-19 2.673201 0.26560425
## 213 2008-01-27 - 2008-02-02 4.433810 0.41434263
## 303 2009-10-18 - 2009-10-24 7.618892 1.00000000
## 304 2009-10-25 - 2009-10-31 7.388359 0.92695883
## 411 2011-11-13 - 2011-11-19 1.462212 0.45551129
##                        Week      ILI   Queries
## 412 2011-11-20 - 2011-11-26 1.655415 0.4130146
## 413 2011-11-27 - 2011-12-03 1.465723 0.4780876
## 414 2011-12-04 - 2011-12-10 1.518106 0.4648074
## 415 2011-12-11 - 2011-12-17 1.663954 0.4794157
## 416 2011-12-18 - 2011-12-24 1.852736 0.5378486
## 417 2011-12-25 - 2011-12-31 2.124130 0.6188579
## 'data.frame':	417 obs. of  3 variables:
##  $ Week   : chr  "2004-01-04 - 2004-01-10" "2004-01-11 - 2004-01-17" "2004-01-18 - 2004-01-24" "2004-01-25 - 2004-01-31" ...
##  $ ILI    : num  2.42 1.81 1.71 1.54 1.44 ...
##  $ Queries: num  0.238 0.22 0.226 0.238 0.224 ...
##  - attr(*, "comment")= chr "glb_trnobs_df"
## NULL
```

```r
# glb_trnobs_df <- read.delim("data/hygiene.txt", header=TRUE, fill=TRUE, sep="\t",
#                             fileEncoding='iso-8859-1')
# glb_trnobs_df <- read.table("data/hygiene.dat.labels", col.names=c("dirty"),
#                             na.strings="[none]")
# glb_trnobs_df$review <- readLines("data/hygiene.dat", n =-1)
# comment(glb_trnobs_df) <- "glb_trnobs_df"                                

# glb_trnobs_df <- data.frame()
# for (symbol in c("Boeing", "CocaCola", "GE", "IBM", "ProcterGamble")) {
#     sym_trnobs_df <- 
#         myimport_data(url=gsub("IBM", symbol, glb_trnng_url), comment="glb_trnobs_df", 
#                                     force_header=TRUE)
#     sym_trnobs_df$Symbol <- symbol
#     glb_trnobs_df <- myrbind_df(glb_trnobs_df, sym_trnobs_df)
# }
                                
# glb_trnobs_df <- 
#     glb_trnobs_df %>% dplyr::filter(Year >= 1999)
                                
if (glb_is_separate_newobs_dataset) {
    glb_newobs_df <- myimport_data(url=glb_newdt_url, comment="glb_newobs_df", 
                                   force_header=TRUE)
    
    # To make plots / stats / checks easier in chunk:inspectORexplore.data
    glb_allobs_df <- myrbind_df(glb_trnobs_df, glb_newobs_df); 
    comment(glb_allobs_df) <- "glb_allobs_df"
} else {
    glb_allobs_df <- glb_trnobs_df; comment(glb_allobs_df) <- "glb_allobs_df"
    if (!glb_split_entity_newobs_datasets) {
        stop("Not implemented yet") 
        glb_newobs_df <- glb_trnobs_df[sample(1:nrow(glb_trnobs_df),
                                          max(2, nrow(glb_trnobs_df) / 1000)),]                    
    } else      if (glb_split_newdata_method == "condition") {
            glb_newobs_df <- do.call("subset", 
                list(glb_trnobs_df, parse(text=glb_split_newdata_condition)))
            glb_trnobs_df <- do.call("subset", 
                list(glb_trnobs_df, parse(text=paste0("!(", 
                                                      glb_split_newdata_condition,
                                                      ")"))))
        } else if (glb_split_newdata_method == "sample") {
                require(caTools)
                
                set.seed(glb_split_sample.seed)
                split <- sample.split(glb_trnobs_df[, glb_rsp_var_raw], 
                                      SplitRatio=(1-glb_split_newdata_size_ratio))
                glb_newobs_df <- glb_trnobs_df[!split, ] 
                glb_trnobs_df <- glb_trnobs_df[split ,]
        } else if (glb_split_newdata_method == "copy") {  
            glb_trnobs_df <- glb_allobs_df
            comment(glb_trnobs_df) <- "glb_trnobs_df"
            glb_newobs_df <- glb_allobs_df
            comment(glb_newobs_df) <- "glb_newobs_df"
        } else stop("glb_split_newdata_method should be %in% c('condition', 'sample', 'copy')")   

    comment(glb_newobs_df) <- "glb_newobs_df"
    myprint_df(glb_newobs_df)
    str(glb_newobs_df)

    if (glb_split_entity_newobs_datasets) {
        myprint_df(glb_trnobs_df)
        str(glb_trnobs_df)        
    }
}         
```

```
## [1] "Reading file ./data/FluTest.csv..."
## [1] "dimensions of data in ./data/FluTest.csv: 52 rows x 3 cols"
##                      Week      ILI   Queries
## 1 2012-01-01 - 2012-01-07 1.766707 0.5936255
## 2 2012-01-08 - 2012-01-14 1.543401 0.4993360
## 3 2012-01-15 - 2012-01-21 1.647615 0.5006640
## 4 2012-01-22 - 2012-01-28 1.684297 0.4794157
## 5 2012-01-29 - 2012-02-04 1.863542 0.4714475
## 6 2012-02-05 - 2012-02-11 1.864079 0.5033201
##                       Week      ILI   Queries
## 1  2012-01-01 - 2012-01-07 1.766707 0.5936255
## 9  2012-02-26 - 2012-03-03 2.095549 0.4608234
## 20 2012-05-13 - 2012-05-19 1.266919 0.3027888
## 24 2012-06-10 - 2012-06-16 1.086121 0.2509960
## 49 2012-12-02 - 2012-12-08 2.978047 0.6719788
## 51 2012-12-16 - 2012-12-22 4.547268 0.7875166
##                       Week      ILI   Queries
## 47 2012-11-18 - 2012-11-24 2.304625 0.5112882
## 48 2012-11-25 - 2012-12-01 2.225997 0.6095618
## 49 2012-12-02 - 2012-12-08 2.978047 0.6719788
## 50 2012-12-09 - 2012-12-15 3.600230 0.7051793
## 51 2012-12-16 - 2012-12-22 4.547268 0.7875166
## 52 2012-12-23 - 2012-12-29 6.033614 0.8054209
## 'data.frame':	52 obs. of  3 variables:
##  $ Week   : chr  "2012-01-01 - 2012-01-07" "2012-01-08 - 2012-01-14" "2012-01-15 - 2012-01-21" "2012-01-22 - 2012-01-28" ...
##  $ ILI    : num  1.77 1.54 1.65 1.68 1.86 ...
##  $ Queries: num  0.594 0.499 0.501 0.479 0.471 ...
##  - attr(*, "comment")= chr "glb_newobs_df"
## NULL
```

```r
if ((num_nas <- sum(is.na(glb_trnobs_df[, glb_rsp_var_raw]))) > 0)
    stop("glb_trnobs_df$", glb_rsp_var_raw, " contains NAs for ", num_nas, " obs")

if (nrow(glb_trnobs_df) == nrow(glb_allobs_df))
    warning("glb_trnobs_df same as glb_allobs_df")
if (nrow(glb_newobs_df) == nrow(glb_allobs_df))
    warning("glb_newobs_df same as glb_allobs_df")

if (length(glb_drop_vars) > 0) {
    warning("dropping vars: ", paste0(glb_drop_vars, collapse=", "))
    glb_allobs_df <- glb_allobs_df[, setdiff(names(glb_allobs_df), glb_drop_vars)]
    glb_trnobs_df <- glb_trnobs_df[, setdiff(names(glb_trnobs_df), glb_drop_vars)]    
    glb_newobs_df <- glb_newobs_df[, setdiff(names(glb_newobs_df), glb_drop_vars)]    
}

#stop(here"); sav_allobs_df <- glb_allobs_df # glb_allobs_df <- sav_allobs_df
# Combine trnent & newobs into glb_allobs_df for easier manipulation
glb_trnobs_df$.src <- "Train"; glb_newobs_df$.src <- "Test"; 
glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, ".src")
glb_allobs_df <- myrbind_df(glb_trnobs_df, glb_newobs_df)
comment(glb_allobs_df) <- "glb_allobs_df"

# Check for duplicates in glb_id_var
if (length(glb_id_var) == 0) {
    warning("using .rownames as identifiers for observations")
    glb_allobs_df$.rownames <- rownames(glb_allobs_df)
    glb_trnobs_df$.rownames <- rownames(subset(glb_allobs_df, .src == "Train"))
    glb_newobs_df$.rownames <- rownames(subset(glb_allobs_df, .src == "Test"))    
    glb_id_var <- ".rownames"
}
if (sum(duplicated(glb_allobs_df[, glb_id_var, FALSE])) > 0)
    stop(glb_id_var, " duplicated in glb_allobs_df")
glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, glb_id_var)

glb_allobs_df <- orderBy(reformulate(glb_id_var), glb_allobs_df)
glb_trnobs_df <- glb_newobs_df <- NULL

glb_chunks_df <- myadd_chunk(glb_chunks_df, "inspect.data", major.inc=TRUE)
```

```
##          label step_major step_minor   bgn   end elapsed
## 1  import.data          1          0 6.863 7.173    0.31
## 2 inspect.data          2          0 7.173    NA      NA
```

## Step `2.0: inspect data`

```r
#print(str(glb_allobs_df))
#View(glb_allobs_df)

dsp_class_dstrb <- function(var) {
    xtab_df <- mycreate_xtab_df(glb_allobs_df, c(".src", var))
    rownames(xtab_df) <- xtab_df$.src
    xtab_df <- subset(xtab_df, select=-.src)
    print(xtab_df)
    print(xtab_df / rowSums(xtab_df, na.rm=TRUE))    
}    

# Performed repeatedly in other chunks
glb_chk_data <- function() {
    # Histogram of predictor in glb_trnobs_df & glb_newobs_df
    print(myplot_histogram(glb_allobs_df, glb_rsp_var_raw) + facet_wrap(~ .src))
    
    if (glb_is_classification) 
        dsp_class_dstrb(var=ifelse(glb_rsp_var %in% names(glb_allobs_df), 
                                   glb_rsp_var, glb_rsp_var_raw))
    mycheck_problem_data(glb_allobs_df)
}
glb_chk_data()
```

![](Google_Flu_template2_files/figure-html/inspect.data-1.png) 

```
## [1] "numeric data missing in glb_allobs_df: "
## named integer(0)
## [1] "numeric data w/ 0s in glb_allobs_df: "
## named integer(0)
## [1] "numeric data w/ Infs in glb_allobs_df: "
## named integer(0)
## [1] "numeric data w/ NaNs in glb_allobs_df: "
## named integer(0)
## [1] "string data missing in glb_allobs_df: "
## Week 
##    0
```

```r
# Create new features that help diagnostics
if (!is.null(glb_map_rsp_raw_to_var)) {
    glb_allobs_df[, glb_rsp_var] <- 
        glb_map_rsp_raw_to_var(glb_allobs_df[, glb_rsp_var_raw])
    mycheck_map_results(mapd_df=glb_allobs_df, 
                        from_col_name=glb_rsp_var_raw, to_col_name=glb_rsp_var)
        
    if (glb_is_classification) dsp_class_dstrb(glb_rsp_var)
}
```

```
##                      Week      ILI   Queries  .src   ILI.log
## 1 2004-01-04 - 2004-01-10 2.418331 0.2377158 Train 0.8830777
## 2 2004-01-11 - 2004-01-17 1.809056 0.2204515 Train 0.5928051
## 3 2004-01-18 - 2004-01-24 1.712024 0.2257636 Train 0.5376762
## 4 2004-01-25 - 2004-01-31 1.542495 0.2377158 Train 0.4334013
## 5 2004-02-01 - 2004-02-07 1.437868 0.2244356 Train 0.3631616
## 6 2004-02-08 - 2004-02-14 1.324274 0.2071713 Train 0.2808644
##                        Week       ILI   Queries  .src     ILI.log
## 153 2006-12-03 - 2006-12-09 1.8596834 0.3559097 Train  0.62040627
## 213 2008-01-27 - 2008-02-02 4.4338100 0.4143426 Train  1.48925927
## 300 2009-09-27 - 2009-10-03 4.6036164 0.6786189 Train  1.52684217
## 329 2010-04-18 - 2010-04-24 1.1620668 0.2602922 Train  0.15020011
## 447 2012-07-22 - 2012-07-28 0.9160412 0.2509960  Test -0.08769394
## 450 2012-08-12 - 2012-08-18 0.9017871 0.2695883  Test -0.10337676
##                        Week      ILI   Queries .src   ILI.log
## 464 2012-11-18 - 2012-11-24 2.304625 0.5112882 Test 0.8349182
## 465 2012-11-25 - 2012-12-01 2.225997 0.6095618 Test 0.8002047
## 466 2012-12-02 - 2012-12-08 2.978047 0.6719788 Test 1.0912677
## 467 2012-12-09 - 2012-12-15 3.600230 0.7051793 Test 1.2809977
## 468 2012-12-16 - 2012-12-22 4.547268 0.7875166 Test 1.5145266
## 469 2012-12-23 - 2012-12-29 6.033614 0.8054209 Test 1.7973462
```

![](Google_Flu_template2_files/figure-html/inspect.data-2.png) 

```r
# check distribution of all numeric data
dsp_numeric_feats_dstrb <- function(feats_vctr) {
    for (feat in feats_vctr) {
        print(sprintf("feat: %s", feat))
        if (glb_is_regression)
            gp <- myplot_scatter(df=glb_allobs_df, ycol_name=glb_rsp_var, xcol_name=feat,
                                 smooth=TRUE)
        if (glb_is_classification)
            gp <- myplot_box(df=glb_allobs_df, ycol_names=feat, xcol_name=glb_rsp_var)
        if (inherits(glb_allobs_df[, feat], "factor"))
            gp <- gp + facet_wrap(reformulate(feat))
        print(gp)
    }
}
# dsp_numeric_vars_dstrb(setdiff(names(glb_allobs_df), 
#                                 union(myfind_chr_cols_df(glb_allobs_df), 
#                                       c(glb_rsp_var_raw, glb_rsp_var))))                                      

add_new_diag_feats <- function(obs_df, ref_df=glb_allobs_df) {
    require(plyr)
    
    obs_df <- mutate(obs_df,
#         <col_name>.NA=is.na(<col_name>),

#         <col_name>.fctr=factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))), 
#         <col_name>.fctr=relevel(factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))),
#                                   "<ref_val>"), 
#         <col2_name>.fctr=relevel(factor(ifelse(<col1_name> == <val>, "<oth_val>", "<ref_val>")), 
#                               as.factor(c("R", "<ref_val>")),
#                               ref="<ref_val>"),

          # This doesn't work - use sapply instead
#         <col_name>.fctr_num=grep(<col_name>, levels(<col_name>.fctr)), 
#         
#         Date.my=as.Date(strptime(Date, "%m/%d/%y %H:%M")),
#         Year=year(Date.my),
#         Month=months(Date.my),
#         Weekday=weekdays(Date.my)

#         <col_name>=<table>[as.character(<col2_name>)],
#         <col_name>=as.numeric(<col2_name>),

#         <col_name> = trunc(<col2_name> / 100),

        .rnorm = rnorm(n=nrow(obs_df))
                        )

    # If levels of a factor are different across obs_df & glb_newobs_df; predict.glm fails  
    # Transformations not handled by mutate
#     obs_df$<col_name>.fctr.num <- sapply(1:nrow(obs_df), 
#         function(row_ix) grep(obs_df[row_ix, "<col_name>"],
#                               levels(obs_df[row_ix, "<col_name>.fctr"])))
    
    #print(summary(obs_df))
    #print(sapply(names(obs_df), function(col) sum(is.na(obs_df[, col]))))
    return(obs_df)
}
glb_allobs_df <- add_new_diag_feats(glb_allobs_df)
```

```
## Loading required package: plyr
```

```r
require(dplyr)
```

```
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:plyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
#stop(here"); sav_allobs_df <- glb_allobs_df # glb_allobs_df <- sav_allobs_df
# Merge some <descriptor>
# glb_allobs_df$<descriptor>.my <- glb_allobs_df$<descriptor>
# glb_allobs_df[grepl("\\bAIRPORT\\b", glb_allobs_df$<descriptor>.my),
#               "<descriptor>.my"] <- "AIRPORT"
# glb_allobs_df$<descriptor>.my <-
#     plyr::revalue(glb_allobs_df$<descriptor>.my, c(
#         "ABANDONED BUILDING" = "OTHER",
#         "##"                      = "##"
#     ))
# print(<descriptor>_freq_df <- mycreate_sqlxtab_df(glb_allobs_df, c("<descriptor>.my")))
# # print(dplyr::filter(<descriptor>_freq_df, grepl("(MEDICAL|DENTAL|OFFICE)", <descriptor>.my)))
# # print(dplyr::filter(dplyr::select(glb_allobs_df, -<var.zoo>), 
# #                     grepl("STORE", <descriptor>.my)))
# glb_exclude_vars_as_features <- c(glb_exclude_vars_as_features, "<descriptor>")

# Check distributions of newly transformed / extracted vars
#   Enhancement: remove vars that were displayed ealier
dsp_numeric_feats_dstrb(feats_vctr=setdiff(names(glb_allobs_df), 
        c(myfind_chr_cols_df(glb_allobs_df), glb_rsp_var_raw, glb_rsp_var, 
          glb_exclude_vars_as_features)))
```

```
## [1] "feat: Queries"
```

![](Google_Flu_template2_files/figure-html/inspect.data-3.png) 

```
## [1] "feat: .rnorm"
```

![](Google_Flu_template2_files/figure-html/inspect.data-4.png) 

```r
#   Convert factors to dummy variables
#   Build splines   require(splines); bsBasis <- bs(training$age, df=3)

#pairs(subset(glb_trnobs_df, select=-c(col_symbol)))
# Check for glb_newobs_df & glb_trnobs_df features range mismatches

# Other diagnostics:
# print(subset(glb_trnobs_df, <col1_name> == max(glb_trnobs_df$<col1_name>, na.rm=TRUE) & 
#                         <col2_name> <= mean(glb_trnobs_df$<col1_name>, na.rm=TRUE)))

# print(glb_trnobs_df[which.max(glb_trnobs_df$<col_name>),])

# print(<col_name>_freq_glb_trnobs_df <- mycreate_tbl_df(glb_trnobs_df, "<col_name>"))
# print(which.min(table(glb_trnobs_df$<col_name>)))
# print(which.max(table(glb_trnobs_df$<col_name>)))
# print(which.max(table(glb_trnobs_df$<col1_name>, glb_trnobs_df$<col2_name>)[, 2]))
# print(table(glb_trnobs_df$<col1_name>, glb_trnobs_df$<col2_name>))
# print(table(is.na(glb_trnobs_df$<col1_name>), glb_trnobs_df$<col2_name>))
# print(table(sign(glb_trnobs_df$<col1_name>), glb_trnobs_df$<col2_name>))
# print(mycreate_xtab_df(glb_trnobs_df, <col1_name>))
# print(mycreate_xtab_df(glb_trnobs_df, c(<col1_name>, <col2_name>)))
# print(<col1_name>_<col2_name>_xtab_glb_trnobs_df <- 
#   mycreate_xtab_df(glb_trnobs_df, c("<col1_name>", "<col2_name>")))
# <col1_name>_<col2_name>_xtab_glb_trnobs_df[is.na(<col1_name>_<col2_name>_xtab_glb_trnobs_df)] <- 0
# print(<col1_name>_<col2_name>_xtab_glb_trnobs_df <- 
#   mutate(<col1_name>_<col2_name>_xtab_glb_trnobs_df, 
#             <col3_name>=(<col1_name> * 1.0) / (<col1_name> + <col2_name>))) 
# print(mycreate_sqlxtab_df(glb_allobs_df, c("<col1_name>", "<col2_name>")))

# print(<col2_name>_min_entity_arr <- 
#    sort(tapply(glb_trnobs_df$<col1_name>, glb_trnobs_df$<col2_name>, min, na.rm=TRUE)))
# print(<col1_name>_na_by_<col2_name>_arr <- 
#    sort(tapply(glb_trnobs_df$<col1_name>.NA, glb_trnobs_df$<col2_name>, mean, na.rm=TRUE)))

# Other plots:
# print(myplot_box(df=glb_trnobs_df, ycol_names="<col1_name>"))
# print(myplot_box(df=glb_trnobs_df, ycol_names="<col1_name>", xcol_name="<col2_name>"))
# print(myplot_line(subset(glb_trnobs_df, Symbol %in% c("CocaCola", "ProcterGamble")), 
#                   "Date.POSIX", "StockPrice", facet_row_colnames="Symbol") + 
#     geom_vline(xintercept=as.numeric(as.POSIXlt("2003-03-01"))) +
#     geom_vline(xintercept=as.numeric(as.POSIXlt("1983-01-01")))        
#         )
# print(myplot_line(subset(glb_trnobs_df, Date.POSIX > as.POSIXct("2004-01-01")), 
#                   "Date.POSIX", "StockPrice") +
#     geom_line(aes(color=Symbol)) + 
#     coord_cartesian(xlim=c(as.POSIXct("1990-01-01"),
#                            as.POSIXct("2000-01-01"))) +     
#     coord_cartesian(ylim=c(0, 250)) +     
#     geom_vline(xintercept=as.numeric(as.POSIXlt("1997-09-01"))) +
#     geom_vline(xintercept=as.numeric(as.POSIXlt("1997-11-01")))        
#         )
# print(myplot_scatter(glb_allobs_df, "<col1_name>", "<col2_name>", smooth=TRUE))
# print(myplot_scatter(glb_allobs_df, "<col1_name>", "<col2_name>", colorcol_name="<Pred.fctr>") + 
#         geom_point(data=subset(glb_allobs_df, <condition>), 
#                     mapping=aes(x=<x_var>, y=<y_var>), color="red", shape=4, size=5) +
#         geom_vline(xintercept=84))

glb_chunks_df <- myadd_chunk(glb_chunks_df, "scrub.data", major.inc=FALSE)
```

```
##          label step_major step_minor   bgn   end elapsed
## 2 inspect.data          2          0 7.173 9.516   2.343
## 3   scrub.data          2          1 9.516    NA      NA
```

### Step `2.1: scrub data`

```r
mycheck_problem_data(glb_allobs_df)
```

```
## [1] "numeric data missing in glb_allobs_df: "
## named integer(0)
## [1] "numeric data w/ 0s in glb_allobs_df: "
## named integer(0)
## [1] "numeric data w/ Infs in glb_allobs_df: "
## named integer(0)
## [1] "numeric data w/ NaNs in glb_allobs_df: "
## named integer(0)
## [1] "string data missing in glb_allobs_df: "
## Week 
##    0
```

```r
dsp_catgs <- function() {
    print("NewsDesk:")
    print(table(glb_allobs_df$NewsDesk))
    print("SectionName:")    
    print(table(glb_allobs_df$SectionName))
    print("SubsectionName:")        
    print(table(glb_allobs_df$SubsectionName))
}

# sel_obs <- function(Popular=NULL, 
#                     NewsDesk=NULL, SectionName=NULL, SubsectionName=NULL,
#         Headline.contains=NULL, Snippet.contains=NULL, Abstract.contains=NULL,
#         Headline.pfx=NULL, NewsDesk.nb=NULL, .clusterid=NULL, myCategory=NULL,
#         perl=FALSE) {
sel_obs <- function(vars_lst) {
    tmp_df <- glb_allobs_df
    # Does not work for Popular == NAs ???
    if (!is.null(Popular)) {
        if (is.na(Popular))
            tmp_df <- tmp_df[is.na(tmp_df$Popular), ] else   
            tmp_df <- tmp_df[tmp_df$Popular == Popular, ]    
    }    
    if (!is.null(NewsDesk)) 
        tmp_df <- tmp_df[tmp_df$NewsDesk == NewsDesk, ]
    if (!is.null(SectionName)) 
        tmp_df <- tmp_df[tmp_df$SectionName == SectionName, ]
    if (!is.null(SubsectionName)) 
        tmp_df <- tmp_df[tmp_df$SubsectionName == SubsectionName, ]
    if (!is.null(Headline.contains))
        tmp_df <- 
            tmp_df[grep(Headline.contains, tmp_df$Headline, perl=perl), ]
    if (!is.null(Snippet.contains))
        tmp_df <- 
            tmp_df[grep(Snippet.contains, tmp_df$Snippet, perl=perl), ]
    if (!is.null(Abstract.contains))
        tmp_df <- 
            tmp_df[grep(Abstract.contains, tmp_df$Abstract, perl=perl), ]
    if (!is.null(Headline.pfx)) {
        if (length(grep("Headline.pfx", names(tmp_df), fixed=TRUE, value=TRUE))
            > 0) tmp_df <- 
                tmp_df[tmp_df$Headline.pfx == Headline.pfx, ] else
        warning("glb_allobs_df does not contain Headline.pfx; ignoring that filter")                    
    }    
    if (!is.null(NewsDesk.nb)) {
        if (any(grepl("NewsDesk.nb", names(tmp_df), fixed=TRUE)) > 0) 
            tmp_df <- 
                tmp_df[tmp_df$NewsDesk.nb == NewsDesk.nb, ] else
        warning("glb_allobs_df does not contain NewsDesk.nb; ignoring that filter")                    
    }    
    if (!is.null(.clusterid)) {
        if (any(grepl(".clusterid", names(tmp_df), fixed=TRUE)) > 0) 
            tmp_df <- 
                tmp_df[tmp_df$clusterid == clusterid, ] else
        warning("glb_allobs_df does not contain clusterid; ignoring that filter")                       }
    if (!is.null(myCategory)) {    
        if (!(myCategory %in% names(glb_allobs_df)))
            tmp_df <-
                tmp_df[tmp_df$myCategory == myCategory, ] else
        warning("glb_allobs_df does not contain myCategory; ignoring that filter")                    
    }    
    
    return(glb_allobs_df$UniqueID %in% tmp_df$UniqueID)
}

dsp_obs <- function(..., cols=c(NULL), all=FALSE) {
    tmp_df <- glb_allobs_df[sel_obs(...), 
                            union(c("UniqueID", "Popular", "myCategory", "Headline"), cols), FALSE]
    if(all) { print(tmp_df) } else { myprint_df(tmp_df) }
}
#dsp_obs(Popular=1, NewsDesk="", SectionName="", Headline.contains="Boehner")
# dsp_obs(Popular=1, NewsDesk="", SectionName="")
# dsp_obs(Popular=NA, NewsDesk="", SectionName="")

dsp_tbl <- function(...) {
    tmp_entity_df <- glb_allobs_df[sel_obs(...), ]
    tmp_tbl <- table(tmp_entity_df$NewsDesk, 
                     tmp_entity_df$SectionName,
                     tmp_entity_df$SubsectionName, 
                     tmp_entity_df$Popular, useNA="ifany")
    #print(names(tmp_tbl))
    #print(dimnames(tmp_tbl))
    print(tmp_tbl)
}

dsp_hdlxtab <- function(str) 
    print(mycreate_sqlxtab_df(glb_allobs_df[sel_obs(Headline.contains=str), ],
                           c("Headline.pfx", "Headline", glb_rsp_var)))
#dsp_hdlxtab("(1914)|(1939)")

dsp_catxtab <- function(str) 
    print(mycreate_sqlxtab_df(glb_allobs_df[sel_obs(Headline.contains=str), ],
        c("Headline.pfx", "NewsDesk", "SectionName", "SubsectionName", glb_rsp_var)))
# dsp_catxtab("1914)|(1939)")
# dsp_catxtab("19(14|39|64):")
# dsp_catxtab("19..:")

# Create myCategory <- NewsDesk#SectionName#SubsectionName
#   Fix some data before merging categories
# glb_allobs_df[sel_obs(Headline.contains="Your Turn:", NewsDesk=""),
#               "NewsDesk"] <- "Styles"
# glb_allobs_df[sel_obs(Headline.contains="School", NewsDesk="", SectionName="U.S.",
#                       SubsectionName=""),
#               "SubsectionName"] <- "Education"
# glb_allobs_df[sel_obs(Headline.contains="Today in Small Business:", NewsDesk="Business"),
#               "SectionName"] <- "Business Day"
# glb_allobs_df[sel_obs(Headline.contains="Today in Small Business:", NewsDesk="Business"),
#               "SubsectionName"] <- "Small Business"
# glb_allobs_df[sel_obs(Headline.contains="Readers Respond:"),
#               "SectionName"] <- "Opinion"
# glb_allobs_df[sel_obs(Headline.contains="Readers Respond:"),
#               "SubsectionName"] <- "Room For Debate"

# glb_allobs_df[sel_obs(NewsDesk="Business", SectionName="", SubsectionName="", Popular=NA),
#               "SubsectionName"] <- "Small Business"
# print(glb_allobs_df[glb_allobs_df$UniqueID %in% c(7973), 
#     c("UniqueID", "Headline", "myCategory", "NewsDesk", "SectionName", "SubsectionName")])
# 
# glb_allobs_df[sel_obs(NewsDesk="Business", SectionName="", SubsectionName=""),
#               "SectionName"] <- "Technology"
# print(glb_allobs_df[glb_allobs_df$UniqueID %in% c(5076, 5736, 5924, 5911, 6532), 
#     c("UniqueID", "Headline", "myCategory", "NewsDesk", "SectionName", "SubsectionName")])
# 
# glb_allobs_df[sel_obs(SectionName="Health"),
#               "NewsDesk"] <- "Science"
# glb_allobs_df[sel_obs(SectionName="Travel"),
#               "NewsDesk"] <- "Travel"
# 
# glb_allobs_df[sel_obs(SubsectionName="Fashion & Style"),
#               "SectionName"] <- ""
# glb_allobs_df[sel_obs(SubsectionName="Fashion & Style"),
#               "SubsectionName"] <- ""
# glb_allobs_df[sel_obs(NewsDesk="Styles", SectionName="", SubsectionName="", Popular=1),
#               "SectionName"] <- "U.S."
# print(glb_allobs_df[glb_allobs_df$UniqueID %in% c(5486), 
#     c("UniqueID", "Headline", "myCategory", "NewsDesk", "SectionName", "SubsectionName")])
# 
# glb_allobs_df$myCategory <- paste(glb_allobs_df$NewsDesk, 
#                                   glb_allobs_df$SectionName,
#                                   glb_allobs_df$SubsectionName,
#                                   sep="#")

# dsp_obs( Headline.contains="Music:"
#         #,NewsDesk=""
#         #,SectionName=""  
#         #,SubsectionName="Fashion & Style"
#         #,Popular=1 #NA
#         ,cols= c("UniqueID", "Headline", "Popular", "myCategory", 
#                 "NewsDesk", "SectionName", "SubsectionName"),
#         all=TRUE)
# dsp_obs( Headline.contains="."
#         ,NewsDesk=""
#         ,SectionName="Opinion"  
#         ,SubsectionName=""
#         #,Popular=1 #NA
#         ,cols= c("UniqueID", "Headline", "Popular", "myCategory", 
#                 "NewsDesk", "SectionName", "SubsectionName"),
#         all=TRUE)
                                        
# Merge some categories
# glb_allobs_df$myCategory <-
#     plyr::revalue(glb_allobs_df$myCategory, c(      
#         "#Business Day#Dealbook"            = "Business#Business Day#Dealbook",
#         "#Business Day#Small Business"      = "Business#Business Day#Small Business",
#         "#Crosswords/Games#"                = "Business#Crosswords/Games#",
#         "Business##"                        = "Business#Technology#",
#         "#Open#"                            = "Business#Technology#",
#         "#Technology#"                      = "Business#Technology#",
#         
#         "#Arts#"                            = "Culture#Arts#",        
#         "Culture##"                         = "Culture#Arts#",        
#         
#         "#World#Asia Pacific"               = "Foreign#World#Asia Pacific",        
#         "Foreign##"                         = "Foreign#World#",    
#         
#         "#N.Y. / Region#"                   = "Metro#N.Y. / Region#",  
#         
#         "#Opinion#"                         = "OpEd#Opinion#",                
#         "OpEd##"                            = "OpEd#Opinion#",        
# 
#         "#Health#"                          = "Science#Health#",
#         "Science##"                         = "Science#Health#",        
#         
#         "Styles##"                          = "Styles##Fashion",                        
#         "Styles#Health#"                    = "Science#Health#",                
#         "Styles#Style#Fashion & Style"      = "Styles##Fashion",        
# 
#         "#Travel#"                          = "Travel#Travel#",                
#         
#         "Magazine#Magazine#"                = "myOther",
#         "National##"                        = "myOther",
#         "National#U.S.#Politics"            = "myOther",        
#         "Sports##"                          = "myOther",
#         "Sports#Sports#"                    = "myOther",
#         "#U.S.#"                            = "myOther",        
#         
# 
# #         "Business##Small Business"        = "Business#Business Day#Small Business",        
# #         
# #         "#Opinion#"                       = "#Opinion#Room For Debate",        
#         "##"                                = "##"
# #         "Business##" = "Business#Business Day#Dealbook",
# #         "Foreign#World#" = "Foreign##",
# #         "#Open#" = "Other",
# #         "#Opinion#The Public Editor" = "OpEd#Opinion#",
# #         "Styles#Health#" = "Styles##",
# #         "Styles#Style#Fashion & Style" = "Styles##",
# #         "#U.S.#" = "#U.S.#Education",
#     ))

# ctgry_xtab_df <- orderBy(reformulate(c("-", ".n")),
#                           mycreate_sqlxtab_df(glb_allobs_df,
#     c("myCategory", "NewsDesk", "SectionName", "SubsectionName", glb_rsp_var)))
# myprint_df(ctgry_xtab_df)
# write.table(ctgry_xtab_df, paste0(glb_out_pfx, "ctgry_xtab.csv"), 
#             row.names=FALSE)

# ctgry_cast_df <- orderBy(~ -Y -NA, dcast(ctgry_xtab_df, 
#                        myCategory + NewsDesk + SectionName + SubsectionName ~ 
#                            Popular.fctr, sum, value.var=".n"))
# myprint_df(ctgry_cast_df)
# write.table(ctgry_cast_df, paste0(glb_out_pfx, "ctgry_cast.csv"), 
#             row.names=FALSE)

# print(ctgry_sum_tbl <- table(glb_allobs_df$myCategory, glb_allobs_df[, glb_rsp_var], 
#                              useNA="ifany"))

dsp_chisq.test <- function(...) {
    sel_df <- glb_allobs_df[sel_obs(...) & 
                            !is.na(glb_allobs_df$Popular), ]
    sel_df$.marker <- 1
    ref_df <- glb_allobs_df[!is.na(glb_allobs_df$Popular), ]
    mrg_df <- merge(ref_df[, c(glb_id_var, "Popular")],
                    sel_df[, c(glb_id_var, ".marker")], all.x=TRUE)
    mrg_df[is.na(mrg_df)] <- 0
    print(mrg_tbl <- table(mrg_df$.marker, mrg_df$Popular))
    print("Rows:Selected; Cols:Popular")
    #print(mrg_tbl)
    print(chisq.test(mrg_tbl))
}
# dsp_chisq.test(Headline.contains="[Ee]bola")
# dsp_chisq.test(Snippet.contains="[Ee]bola")
# dsp_chisq.test(Abstract.contains="[Ee]bola")

# print(mycreate_sqlxtab_df(glb_allobs_df[sel_obs(Headline.contains="[Ee]bola"), ], 
#                           c(glb_rsp_var, "NewsDesk", "SectionName", "SubsectionName")))

# print(table(glb_allobs_df$NewsDesk, glb_allobs_df$SectionName))
# print(table(glb_allobs_df$SectionName, glb_allobs_df$SubsectionName))
# print(table(glb_allobs_df$NewsDesk, glb_allobs_df$SectionName, glb_allobs_df$SubsectionName))

# glb_allobs_df$myCategory.fctr <- as.factor(glb_allobs_df$myCategory)
# glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
#                                       c("myCategory", "NewsDesk", "SectionName", "SubsectionName"))

# Copy Headline into Snipper & Abstract if they are empty
# print(glb_allobs_df[nchar(glb_allobs_df[, "Snippet"]) == 0, c("Headline", "Snippet")])
# print(glb_allobs_df[glb_allobs_df$Headline == glb_allobs_df$Snippet, 
#                     c("UniqueID", "Headline", "Snippet")])
# glb_allobs_df[nchar(glb_allobs_df[, "Snippet"]) == 0, "Snippet"] <- 
#     glb_allobs_df[nchar(glb_allobs_df[, "Snippet"]) == 0, "Headline"]
# 
# print(glb_allobs_df[nchar(glb_allobs_df[, "Abstract"]) == 0, c("Headline", "Abstract")])
# print(glb_allobs_df[glb_allobs_df$Headline == glb_allobs_df$Abstract, 
#                     c("UniqueID", "Headline", "Abstract")])
# glb_allobs_df[nchar(glb_allobs_df[, "Abstract"]) == 0, "Abstract"] <- 
#     glb_allobs_df[nchar(glb_allobs_df[, "Abstract"]) == 0, "Headline"]

# WordCount_0_df <- subset(glb_allobs_df, WordCount == 0)
# table(WordCount_0_df$Popular, WordCount_0_df$WordCount, useNA="ifany")
# myprint_df(WordCount_0_df[, 
#                 c("UniqueID", "Popular", "WordCount", "Headline")])
```

### Step `2.1: scrub data`

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "transform.data", major.inc=FALSE)
```

```
##            label step_major step_minor    bgn    end elapsed
## 3     scrub.data          2          1  9.516 10.161   0.645
## 4 transform.data          2          2 10.161     NA      NA
```

```r
### Mapping dictionary
#sav_allobs_df <- glb_allobs_df; glb_allobs_df <- sav_allobs_df
if (!is.null(glb_map_vars)) {
    for (feat in glb_map_vars) {
        map_df <- myimport_data(url=glb_map_urls[[feat]], 
                                            comment="map_df", 
                                           print_diagn=TRUE)
        glb_allobs_df <- mymap_codes(glb_allobs_df, feat, names(map_df)[2], 
                                     map_df, map_join_col_name=names(map_df)[1], 
                                     map_tgt_col_name=names(map_df)[2])
    }
    glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, glb_map_vars)
}

### Forced Assignments
#stop(here"); sav_allobs_df <- glb_allobs_df; glb_allobs_df <- sav_allobs_df
for (feat in glb_assign_vars) {
    new_feat <- paste0(feat, ".my")
    print(sprintf("Forced Assignments for: %s -> %s...", feat, new_feat))
    glb_allobs_df[, new_feat] <- glb_allobs_df[, feat]
    
    pairs <- glb_assign_pairs_lst[[feat]]
    for (pair_ix in 1:length(pairs$from)) {
        if (is.na(pairs$from[pair_ix]))
            nobs <- nrow(filter(glb_allobs_df, 
                                is.na(eval(parse(text=feat),
                                            envir=glb_allobs_df)))) else
            nobs <- sum(glb_allobs_df[, feat] == pairs$from[pair_ix])
        #nobs <- nrow(filter(glb_allobs_df, is.na(Married.fctr)))    ; print(nobs)
        
        if ((is.na(pairs$from[pair_ix])) && (is.na(pairs$to[pair_ix])))
            stop("what are you trying to do ???")
        if (is.na(pairs$from[pair_ix]))
            glb_allobs_df[is.na(glb_allobs_df[, feat]), new_feat] <- 
                pairs$to[pair_ix] else
            glb_allobs_df[glb_allobs_df[, feat] == pairs$from[pair_ix], new_feat] <- 
                pairs$to[pair_ix]
                    
        print(sprintf("    %s -> %s for %s obs", 
                      pairs$from[pair_ix], pairs$to[pair_ix], format(nobs, big.mark=",")))
    }

    glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, glb_assign_vars)
}

### Derivations using mapping functions
#stop(here"); sav_allobs_df <- glb_allobs_df; glb_allobs_df <- sav_allobs_df
for (new_feat in glb_derive_vars) {
    print(sprintf("Creating new feature: %s...", new_feat))
    args_lst <- NULL 
    for (arg in glb_derive_lst[[new_feat]]$args) 
        args_lst[[arg]] <- glb_allobs_df[, arg]
    glb_allobs_df[, new_feat] <- do.call(glb_derive_lst[[new_feat]]$mapfn, args_lst)
}
```

```
## [1] "Creating new feature: Week.bgn..."
## [1] "Creating new feature: Week.end..."
## [1] "Creating new feature: ILI.2.lag..."
## [1] "Creating new feature: ILI.2.lag.log..."
```

## Step `2.2: transform data`

```r
#```{r extract_features, cache=FALSE, eval=!is.null(glb_txt_vars)}
glb_chunks_df <- myadd_chunk(glb_chunks_df, "extract.features", major.inc=TRUE)
```

```
##              label step_major step_minor    bgn    end elapsed
## 4   transform.data          2          2 10.161 10.224   0.063
## 5 extract.features          3          0 10.224     NA      NA
```

```r
extract.features_chunk_df <- myadd_chunk(NULL, "extract.features_bgn")
```

```
##                  label step_major step_minor   bgn end elapsed
## 1 extract.features_bgn          1          0 10.23  NA      NA
```

```r
# Options:
#   Select Tf, log(1 + Tf), Tf-IDF or BM25Tf-IDf

# Create new features that help prediction
# <col_name>.lag.2 <- lag(zoo(glb_trnobs_df$<col_name>), -2, na.pad=TRUE)
# glb_trnobs_df[, "<col_name>.lag.2"] <- coredata(<col_name>.lag.2)
# <col_name>.lag.2 <- lag(zoo(glb_newobs_df$<col_name>), -2, na.pad=TRUE)
# glb_newobs_df[, "<col_name>.lag.2"] <- coredata(<col_name>.lag.2)
# 
# glb_newobs_df[1, "<col_name>.lag.2"] <- glb_trnobs_df[nrow(glb_trnobs_df) - 1, 
#                                                    "<col_name>"]
# glb_newobs_df[2, "<col_name>.lag.2"] <- glb_trnobs_df[nrow(glb_trnobs_df), 
#                                                    "<col_name>"]
                                                   
# glb_allobs_df <- mutate(glb_allobs_df,
#     A.P.http=ifelse(grepl("http",Added,fixed=TRUE), 1, 0)
#                     )
# 
# glb_trnobs_df <- mutate(glb_trnobs_df,
#                     )
# 
# glb_newobs_df <- mutate(glb_newobs_df,
#                     )

#   Convert dates to numbers 
#       typically, dates come in as chars; 
#           so this must be done before converting chars to factors

#stop(here"); sav_allobs_df <- glb_allobs_df #; glb_allobs_df <- sav_allobs_df
if (!is.null(glb_date_vars)) {
    glb_allobs_df <- cbind(glb_allobs_df, 
        myextract_dates_df(df=glb_allobs_df, vars=glb_date_vars, 
                           id_vars=glb_id_var, rsp_var=glb_rsp_var))
    for (sfx in c("", ".POSIX"))
        glb_exclude_vars_as_features <- 
            union(glb_exclude_vars_as_features, 
                    paste(glb_date_vars, sfx, sep=""))

    for (feat in glb_date_vars) {
        glb_allobs_df <- orderBy(reformulate(paste0(feat, ".POSIX")), glb_allobs_df)
#         print(myplot_scatter(glb_allobs_df, xcol_name=paste0(feat, ".POSIX"),
#                              ycol_name=glb_rsp_var, colorcol_name=glb_rsp_var))
        print(myplot_scatter(glb_allobs_df[glb_allobs_df[, paste0(feat, ".POSIX")] >=
                                               strptime("2012-12-01", "%Y-%m-%d"), ], 
                             xcol_name=paste0(feat, ".POSIX"),
                             ycol_name=glb_rsp_var, colorcol_name=paste0(feat, ".wkend")))

        # Create features that measure the gap between previous timestamp in the data
        require(zoo)
        z <- zoo(as.numeric(as.POSIXlt(glb_allobs_df[, paste0(feat, ".POSIX")])))
        glb_allobs_df[, paste0(feat, ".zoo")] <- z
        print(head(glb_allobs_df[, c(glb_id_var, feat, paste0(feat, ".zoo"))]))
        print(myplot_scatter(glb_allobs_df[glb_allobs_df[,  paste0(feat, ".POSIX")] >
                                            strptime("2012-10-01", "%Y-%m-%d"), ], 
                            xcol_name=paste0(feat, ".zoo"), ycol_name=glb_rsp_var,
                            colorcol_name=glb_rsp_var))
        b <- zoo(, seq(nrow(glb_allobs_df)))
        
        last1 <- as.numeric(merge(z-lag(z, -1), b, all=TRUE)); last1[is.na(last1)] <- 0
        glb_allobs_df[, paste0(feat, ".last1.log")] <- log(1 + last1)
        print(gp <- myplot_box(df=glb_allobs_df[glb_allobs_df[, 
                                                    paste0(feat, ".last1.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last1.log"), 
                               xcol_name=glb_rsp_var))
        
        last2 <- as.numeric(merge(z-lag(z, -2), b, all=TRUE)); last2[is.na(last2)] <- 0
        glb_allobs_df[, paste0(feat, ".last2.log")] <- log(1 + last2)
        print(gp <- myplot_box(df=glb_allobs_df[glb_allobs_df[, 
                                                    paste0(feat, ".last2.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last2.log"), 
                               xcol_name=glb_rsp_var))
        
        last10 <- as.numeric(merge(z-lag(z, -10), b, all=TRUE)); last10[is.na(last10)] <- 0
        glb_allobs_df[, paste0(feat, ".last10.log")] <- log(1 + last10)
        print(gp <- myplot_box(df=glb_allobs_df[glb_allobs_df[, 
                                                    paste0(feat, ".last10.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last10.log"), 
                               xcol_name=glb_rsp_var))
        
        last100 <- as.numeric(merge(z-lag(z, -100), b, all=TRUE)); last100[is.na(last100)] <- 0
        glb_allobs_df[, paste0(feat, ".last100.log")] <- log(1 + last100)
        print(gp <- myplot_box(df=glb_allobs_df[glb_allobs_df[, 
                                                    paste0(feat, ".last100.log")] > 0, ], 
                               ycol_names=paste0(feat, ".last100.log"), 
                               xcol_name=glb_rsp_var))
        
        glb_allobs_df <- orderBy(reformulate(glb_id_var), glb_allobs_df)
        glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
                                                c(paste0(feat, ".zoo")))
        # all2$last3 = as.numeric(merge(z-lag(z, -3), b, all = TRUE))
        # all2$last5 = as.numeric(merge(z-lag(z, -5), b, all = TRUE))
        # all2$last10 = as.numeric(merge(z-lag(z, -10), b, all = TRUE))
        # all2$last20 = as.numeric(merge(z-lag(z, -20), b, all = TRUE))
        # all2$last50 = as.numeric(merge(z-lag(z, -50), b, all = TRUE))
        # 
        # 
        # # order table
        # all2 = all2[order(all2$id),]
        # 
        # ## fill in NAs
        # # count averages
        # na.avg = all2 %>% group_by(weekend, hour) %>% dplyr::summarise(
        #     last1=mean(last1, na.rm=TRUE),
        #     last3=mean(last3, na.rm=TRUE),
        #     last5=mean(last5, na.rm=TRUE),
        #     last10=mean(last10, na.rm=TRUE),
        #     last20=mean(last20, na.rm=TRUE),
        #     last50=mean(last50, na.rm=TRUE)
        # )
        # 
        # # fill in averages
        # na.merge = merge(all2, na.avg, by=c("weekend","hour"))
        # na.merge = na.merge[order(na.merge$id),]
        # for(i in c("last1", "last3", "last5", "last10", "last20", "last50")) {
        #     y = paste0(i, ".y")
        #     idx = is.na(all2[[i]])
        #     all2[idx,][[i]] <- na.merge[idx,][[y]]
        # }
        # rm(na.avg, na.merge, b, i, idx, n, pd, sec, sh, y, z)
    }
}
```

```
## Loading required package: XML
```

```
## Warning in myplot_scatter(glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".POSIX")] >= : converting Week.bgn.wkend to class:factor
```

```
##                      Week   Week.bgn Week.bgn.zoo
## 1 2004-01-04 - 2004-01-10 2004-01-04   1073192400
## 2 2004-01-11 - 2004-01-17 2004-01-11   1073797200
## 3 2004-01-18 - 2004-01-24 2004-01-18   1074402000
## 4 2004-01-25 - 2004-01-31 2004-01-25   1075006800
## 5 2004-02-01 - 2004-02-07 2004-02-01   1075611600
## 6 2004-02-08 - 2004-02-14 2004-02-08   1076216400
```

```
## Warning in myplot_scatter(glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".POSIX")] > : converting ILI.log to class:factor
```

![](Google_Flu_template2_files/figure-html/extract.features-1.png) 

```
## Don't know how to automatically pick scale for object of type zoo. Defaulting to continuous
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last1.log")] > : xcol_name:ILI.log is not a factor; creating ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-2.png) 

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last2.log")] > : xcol_name:ILI.log is not a factor; creating ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-3.png) 

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last10.log")] > : xcol_name:ILI.log is not a factor; creating
## ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-4.png) 

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last100.log")] > : xcol_name:ILI.log is not a factor; creating
## ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-5.png) 

```
## Warning in myplot_scatter(glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".POSIX")] >= : converting Week.end.wkend to class:factor
```

![](Google_Flu_template2_files/figure-html/extract.features-6.png) 

```
##                      Week   Week.end Week.end.zoo
## 1 2004-01-04 - 2004-01-10 2004-01-10   1073192400
## 2 2004-01-11 - 2004-01-17 2004-01-17   1073797200
## 3 2004-01-18 - 2004-01-24 2004-01-24   1074402000
## 4 2004-01-25 - 2004-01-31 2004-01-31   1075006800
## 5 2004-02-01 - 2004-02-07 2004-02-07   1075611600
## 6 2004-02-08 - 2004-02-14 2004-02-14   1076216400
```

```
## Warning in myplot_scatter(glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".POSIX")] > : converting ILI.log to class:factor
```

![](Google_Flu_template2_files/figure-html/extract.features-7.png) 

```
## Don't know how to automatically pick scale for object of type zoo. Defaulting to continuous
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last1.log")] > : xcol_name:ILI.log is not a factor; creating ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-8.png) 

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last2.log")] > : xcol_name:ILI.log is not a factor; creating ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-9.png) 

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last10.log")] > : xcol_name:ILI.log is not a factor; creating
## ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-10.png) 

```
## Warning in myplot_box(df = glb_allobs_df[glb_allobs_df[, paste0(feat,
## ".last100.log")] > : xcol_name:ILI.log is not a factor; creating
## ILI.log_fctr
```

![](Google_Flu_template2_files/figure-html/extract.features-11.png) ![](Google_Flu_template2_files/figure-html/extract.features-12.png) 

```r
rm(last1, last10, last100)

#   Create factors of string variables
extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "factorize.str.vars"), major.inc=TRUE)
```

```
##                                 label step_major step_minor     bgn
## 1                extract.features_bgn          1          0  10.230
## 2 extract.features_factorize.str.vars          2          0 108.047
##       end elapsed
## 1 108.046  97.816
## 2      NA      NA
```

```r
#stop(here"); sav_allobs_df <- glb_allobs_df; #glb_allobs_df <- sav_allobs_df
print(str_vars <- myfind_chr_cols_df(glb_allobs_df))
```

```
##       Week       .src   Week.bgn   Week.end 
##     "Week"     ".src" "Week.bgn" "Week.end"
```

```r
if (length(str_vars <- setdiff(str_vars, 
                               c(glb_exclude_vars_as_features, glb_txt_vars))) > 0) {
    for (var in str_vars) {
        warning("Creating factors of string variable: ", var, 
                ": # of unique values: ", length(unique(glb_allobs_df[, var])))
        glb_allobs_df[, paste0(var, ".fctr")] <- 
            relevel(factor(glb_allobs_df[, var]),
                    names(which.max(table(glb_allobs_df[, var], useNA = "ifany"))))
    }
    glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, str_vars)
}

if (!is.null(glb_txt_vars)) {
    require(foreach)
    require(gsubfn)
    require(stringr)
    require(tm)
    
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "process.text"), major.inc=TRUE)
    
    chk_pattern_freq <- function(rex_str, ignore.case=TRUE) {
        match_mtrx <- str_extract_all(txt_vctr, regex(rex_str, ignore_case=ignore.case), 
                                      simplify=TRUE)
        match_df <- as.data.frame(match_mtrx[match_mtrx != ""])
        names(match_df) <- "pattern"
        return(mycreate_sqlxtab_df(match_df, "pattern"))        
    }

#     match_lst <- gregexpr("\\bok(?!ay)", txt_vctr[746], ignore.case = FALSE, perl=TRUE); print(match_lst)
    dsp_pattern <- function(rex_str, ignore.case=TRUE, print.all=TRUE) {
        match_lst <- gregexpr(rex_str, txt_vctr, ignore.case = ignore.case, perl=TRUE)
        match_lst <- regmatches(txt_vctr, match_lst)
        match_df <- data.frame(matches=sapply(match_lst, 
                                              function (elems) paste(elems, collapse="#")))
        match_df <- subset(match_df, matches != "")
        if (print.all)
            print(match_df)
        return(match_df)
    }
    
    dsp_matches <- function(rex_str, ix) {
        print(match_pos <- gregexpr(rex_str, txt_vctr[ix], perl=TRUE))
        print(str_sub(txt_vctr[ix], (match_pos[[1]] / 100) *  99 +   0, 
                                    (match_pos[[1]] / 100) * 100 + 100))        
    }

    myapply_gsub <- function(...) {
        if ((length_lst <- length(names(gsub_map_lst))) == 0)
            return(txt_vctr)
        for (ptn_ix in 1:length_lst) {
            if ((ptn_ix %% 10) == 0)
                print(sprintf("running gsub for %02d (of %02d): #%s#...", ptn_ix, 
                                length(names(gsub_map_lst)), names(gsub_map_lst)[ptn_ix]))
            txt_vctr <- gsub(names(gsub_map_lst)[ptn_ix], gsub_map_lst[[ptn_ix]], 
                               txt_vctr, ...)
        }
        return(txt_vctr)
    }    

    myapply_txtmap <- function(txt_vctr, ...) {
        nrows <- nrow(glb_txt_map_df)
        for (ptn_ix in 1:nrows) {
            if ((ptn_ix %% 10) == 0)
                print(sprintf("running gsub for %02d (of %02d): #%s#...", ptn_ix, 
                                nrows, glb_txt_map_df[ptn_ix, "rex_str"]))
            txt_vctr <- gsub(glb_txt_map_df[ptn_ix, "rex_str"], 
                             glb_txt_map_df[ptn_ix, "rpl_str"], 
                               txt_vctr, ...)
        }
        return(txt_vctr)
    }    

    chk.equal <- function(bgn, end) {
        print(all.equal(sav_txt_lst[["Headline"]][bgn:end], 
                        glb_txt_lst[["Headline"]][bgn:end]))
    }    
    dsp.equal <- function(bgn, end) {
        print(sav_txt_lst[["Headline"]][bgn:end])
        print(glb_txt_lst[["Headline"]][bgn:end])
    }    
#sav_txt_lst <- glb_txt_lst; all.equal(sav_txt_lst, glb_txt_lst)
#all.equal(sav_txt_lst[["Headline"]][1:4200], glb_txt_lst[["Headline"]][1:4200])
#chk.equal( 1, 100)
#dsp.equal(86, 90)
    
    glb_txt_map_df <- read.csv("mytxt_map.csv", comment.char="#", strip.white=TRUE)
    glb_txt_lst <- list(); 
    print(sprintf("Building glb_txt_lst..."))
    glb_txt_lst <- foreach(txt_var=glb_txt_vars) %dopar% {   
#     for (txt_var in glb_txt_vars) {
        txt_vctr <- glb_allobs_df[, txt_var]
        
        # myapply_txtmap shd be created as a tm_map::content_transformer ?
        #print(glb_txt_map_df)
        #txt_var=glb_txt_vars[3]; txt_vctr <- glb_txt_lst[[txt_var]]
        #print(rex_str <- glb_txt_map_df[163, "rex_str"])
        #print(rex_str <- glb_txt_map_df[glb_txt_map_df$rex_str == "\\bWall St\\.", "rex_str"])
        #print(rex_str <- glb_txt_map_df[grepl("du Pont", glb_txt_map_df$rex_str), "rex_str"])        
        #print(rex_str <- glb_txt_map_df[glb_txt_map_df$rpl_str == "versus", "rex_str"])             
        #print(tmp_vctr <- grep(rex_str, txt_vctr, value=TRUE, ignore.case=FALSE))
        #ret_lst <- regexec(rex_str, txt_vctr, ignore.case=FALSE); ret_lst <- regmatches(txt_vctr, ret_lst); ret_vctr <- sapply(1:length(ret_lst), function(pos_ix) ifelse(length(ret_lst[[pos_ix]]) > 0, ret_lst[[pos_ix]], "")); print(ret_vctr <- ret_vctr[ret_vctr != ""])
        #gsub(rex_str, glb_txt_map_df[glb_txt_map_df$rex_str == rex_str, "rpl_str"], tmp_vctr, ignore.case=FALSE)
        #grep("Hong Hong", txt_vctr, value=TRUE)
    
        txt_vctr <- myapply_txtmap(txt_vctr, ignore.case=FALSE)    
    }
    names(glb_txt_lst) <- glb_txt_vars

    for (txt_var in glb_txt_vars) {
        print(sprintf("Remaining OK in %s:", txt_var))
        txt_vctr <- glb_txt_lst[[txt_var]]
        
        print(chk_pattern_freq(rex_str <- "(?<!(BO|HO|LO))OK(?!(E\\!|ED|IE|IN|S ))",
                               ignore.case=FALSE))
        match_df <- dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
        for (row in row.names(match_df))
            dsp_matches(rex_str, ix=as.numeric(row))

        print(chk_pattern_freq(rex_str <- "Ok(?!(a\\.|ay|in|ra|um))", ignore.case=FALSE))
        match_df <- dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
        for (row in row.names(match_df))
            dsp_matches(rex_str, ix=as.numeric(row))

        print(chk_pattern_freq(rex_str <- "(?<!( b| B| c| C| g| G| j| M| p| P| w| W| r| Z|\\(b|ar|bo|Bo|co|Co|Ew|gk|go|ho|ig|jo|kb|ke|Ke|ki|lo|Lo|mo|mt|no|No|po|ra|ro|sm|Sm|Sp|to|To))ok(?!(ay|bo|e |e\\)|e,|e\\.|eb|ed|el|en|er|es|ey|i |ie|in|it|ka|ke|ki|ly|on|oy|ra|st|u |uc|uy|yl|yo))",
                               ignore.case=FALSE))
        match_df <- dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
        for (row in row.names(match_df))
            dsp_matches(rex_str, ix=as.numeric(row))
    }    
    # txt_vctr <- glb_txt_lst[[glb_txt_vars[1]]]
    # print(chk_pattern_freq(rex_str <- "(?<!( b| c| C| p|\\(b|bo|co|lo|Lo|Sp|to|To))ok(?!(ay|e |e\\)|e,|e\\.|ed|el|en|es|ey|ie|in|on|ra))", ignore.case=FALSE))
    # print(chk_pattern_freq(rex_str <- "ok(?!(ay|el|on|ra))", ignore.case=FALSE))
    # dsp_pattern(rex_str, ignore.case=FALSE, print.all=FALSE)
    # dsp_matches(rex_str, ix=8)
    # substr(txt_vctr[86], 5613, 5620)
    # substr(glb_allobs_df[301, "review"], 550, 650)

#stop(here"); sav_txt_lst <- glb_txt_lst    
    for (txt_var in glb_txt_vars) {
        print(sprintf("Remaining Acronyms in %s:", txt_var))
        txt_vctr <- glb_txt_lst[[txt_var]]
        
        print(chk_pattern_freq(rex_str <- "([[:upper:]]\\.( *)){2,}", ignore.case=FALSE))
        
        # Check for names
        print(subset(chk_pattern_freq(rex_str <- "(([[:upper:]]+)\\.( *)){1}",
                                      ignore.case=FALSE),
                     .n > 1))
        # dsp_pattern(rex_str="(OK\\.( *)){1}", ignore.case=FALSE)
        # dsp_matches(rex_str="(OK\\.( *)){1}", ix=557)
        #dsp_matches(rex_str="\\bR\\.I\\.P(\\.*)(\\B)", ix=461)
        #dsp_matches(rex_str="\\bR\\.I\\.P(\\.*)", ix=461)        
        #print(str_sub(txt_vctr[676], 10100, 10200))
        #print(str_sub(txt_vctr[74], 1, -1))        
    }

    for (txt_var in glb_txt_vars) {
        re_str <- "\\b(Fort|Ft\\.|Hong|Las|Los|New|Puerto|Saint|San|St\\.)( |-)(\\w)+"
        print(sprintf("Remaining #%s# terms in %s: ", re_str, txt_var))
        txt_vctr <- glb_txt_lst[[txt_var]]        
        print(orderBy(~ -.n +pattern, subset(chk_pattern_freq(re_str, ignore.case=FALSE), 
                                             grepl("( |-)[[:upper:]]", pattern))))
        print("    consider cleaning if relevant to problem domain; geography name; .n > 1")
        #grep("New G", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("St\\. Wins", txt_vctr, value=TRUE, ignore.case=FALSE)
    }        
        
#stop(here"); sav_txt_lst <- glb_txt_lst    
    for (txt_var in glb_txt_vars) {
        re_str <- "\\b(N|S|E|W|C)( |\\.)(\\w)+"
        print(sprintf("Remaining #%s# terms in %s: ", re_str, txt_var))        
        txt_vctr <- glb_txt_lst[[txt_var]]                
        print(orderBy(~ -.n +pattern, subset(chk_pattern_freq(re_str, ignore.case=FALSE), 
                                             grepl(".", pattern))))
        #grep("N Weaver", txt_vctr, value=TRUE, ignore.case=FALSE)        
    }    

    for (txt_var in glb_txt_vars) {
        re_str <- "\\b(North|South|East|West|Central)( |\\.)(\\w)+"
        print(sprintf("Remaining #%s# terms in %s: ", re_str, txt_var))        
        txt_vctr <- glb_txt_lst[[txt_var]]                        
        print(orderBy(~ -.n +pattern, subset(chk_pattern_freq(re_str, ignore.case=FALSE), 
                                             grepl(".", pattern))))
        #grep("Central (African|Bankers|Cast|Italy|Role|Spring)", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("East (Africa|Berlin|London|Poland|Rivals|Spring)", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("North (American|Korean|West)", txt_vctr, value=TRUE, ignore.case=FALSE)        
        #grep("South (Pacific|Street)", txt_vctr, value=TRUE, ignore.case=FALSE)
        #grep("St\\. Martins", txt_vctr, value=TRUE, ignore.case=FALSE)
    }    

    find_cmpnd_wrds <- function(txt_vctr) {
        txt_corpus <- Corpus(VectorSource(txt_vctr))
        txt_corpus <- tm_map(txt_corpus, tolower)
        txt_corpus <- tm_map(txt_corpus, PlainTextDocument)
        txt_corpus <- tm_map(txt_corpus, removePunctuation, 
                             preserve_intra_word_dashes=TRUE)
        full_Tf_DTM <- DocumentTermMatrix(txt_corpus, 
                                          control=list(weighting=weightTf))
        print("   Full TermMatrix:"); print(full_Tf_DTM)
        full_Tf_mtrx <- as.matrix(full_Tf_DTM)
        rownames(full_Tf_mtrx) <- rownames(glb_allobs_df) # print undreadable otherwise
        full_Tf_vctr <- colSums(full_Tf_mtrx)
        names(full_Tf_vctr) <- dimnames(full_Tf_DTM)[[2]]
        #grep("year", names(full_Tf_vctr), value=TRUE)
        #which.max(full_Tf_mtrx[, "yearlong"])
        full_Tf_df <- as.data.frame(full_Tf_vctr)
        names(full_Tf_df) <- "Tf.full"
        full_Tf_df$term <- rownames(full_Tf_df)
        #full_Tf_df$freq.full <- colSums(full_Tf_mtrx != 0)
        full_Tf_df <- orderBy(~ -Tf.full, full_Tf_df)
        cmpnd_Tf_df <- full_Tf_df[grep("-", full_Tf_df$term, value=TRUE) ,]
        
        filter_df <- read.csv("mytxt_compound.csv", comment.char="#", strip.white=TRUE)
        cmpnd_Tf_df$filter <- FALSE
        for (row_ix in 1:nrow(filter_df))
            cmpnd_Tf_df[!cmpnd_Tf_df$filter, "filter"] <- 
            grepl(filter_df[row_ix, "rex_str"], 
                  cmpnd_Tf_df[!cmpnd_Tf_df$filter, "term"], ignore.case=TRUE)
        cmpnd_Tf_df <- subset(cmpnd_Tf_df, !filter)
        # Bug in tm_map(txt_corpus, removePunctuation, preserve_intra_word_dashes=TRUE) ???
        #   "net-a-porter" gets converted to "net-aporter"
        #grep("net-a-porter", txt_vctr, ignore.case=TRUE, value=TRUE)
        #grep("maser-laser", txt_vctr, ignore.case=TRUE, value=TRUE)
        #txt_corpus[[which(grepl("net-a-porter", txt_vctr, ignore.case=TRUE))]]
        #grep("\\b(across|longer)-(\\w)", cmpnd_Tf_df$term, ignore.case=TRUE, value=TRUE)
        #grep("(\\w)-(affected|term)\\b", cmpnd_Tf_df$term, ignore.case=TRUE, value=TRUE)
        
        print(sprintf("nrow(cmpnd_Tf_df): %d", nrow(cmpnd_Tf_df)))
        myprint_df(cmpnd_Tf_df)
    }

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "process.text_reporting_compound_terms"), major.inc=FALSE)
    
    for (txt_var in glb_txt_vars) {
        print(sprintf("Remaining compound terms in %s: ", txt_var))        
        txt_vctr <- glb_txt_lst[[txt_var]]                        
#         find_cmpnd_wrds(txt_vctr)
        #grep("thirty-five", txt_vctr, ignore.case=TRUE, value=TRUE)
        #rex_str <- glb_txt_map_df[grepl("hirty", glb_txt_map_df$rex_str), "rex_str"]
    }

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "build.corpus"), major.inc=TRUE)
    
    glb_corpus_lst <- list()
    print(sprintf("Building glb_corpus_lst..."))
    glb_corpus_lst <- foreach(txt_var=glb_txt_vars) %dopar% {   
#     for (txt_var in glb_txt_vars) {
        txt_corpus <- Corpus(VectorSource(glb_txt_lst[[txt_var]]))
        txt_corpus <- tm_map(txt_corpus, tolower) #nuppr
        txt_corpus <- tm_map(txt_corpus, PlainTextDocument)
        txt_corpus <- tm_map(txt_corpus, removePunctuation) #npnct<chr_ix>
#         txt-corpus <- tm_map(txt_corpus, content_transformer(function(x, pattern) gsub(pattern, "", x))   

        # Not to be run in production
        inspect_terms <- function() {
            full_Tf_DTM <- DocumentTermMatrix(txt_corpus, 
                                              control=list(weighting=weightTf))
            print("   Full TermMatrix:"); print(full_Tf_DTM)
            full_Tf_mtrx <- as.matrix(full_Tf_DTM)
            rownames(full_Tf_mtrx) <- rownames(glb_allobs_df) # print undreadable otherwise
            full_Tf_vctr <- colSums(full_Tf_mtrx)
            names(full_Tf_vctr) <- dimnames(full_Tf_DTM)[[2]]
            #grep("year", names(full_Tf_vctr), value=TRUE)
            #which.max(full_Tf_mtrx[, "yearlong"])
            full_Tf_df <- as.data.frame(full_Tf_vctr)
            names(full_Tf_df) <- "Tf.full"
            full_Tf_df$term <- rownames(full_Tf_df)
            #full_Tf_df$freq.full <- colSums(full_Tf_mtrx != 0)
            full_Tf_df <- orderBy(~ -Tf.full +term, full_Tf_df)
            print(myplot_histogram(full_Tf_df, "Tf.full"))
            myprint_df(full_Tf_df)
            #txt_corpus[[which(grepl("zun", txt_vctr, ignore.case=TRUE))]]
            digit_terms_df <- subset(full_Tf_df, grepl("[[:digit:]]", term))
            myprint_df(digit_terms_df)
            return(full_Tf_df)
        }    
        #print("RemovePunct:"); remove_punct_Tf_df <- inspect_terms()

        txt_corpus <- tm_map(txt_corpus, removeWords, 
                             c(glb_append_stop_words[[txt_var]], 
                               stopwords("english"))) #nstopwrds
        #print("StoppedWords:"); stopped_words_Tf_df <- inspect_terms()
        txt_corpus <- tm_map(txt_corpus, stemDocument) #Features for lost information: Difference/ratio in density of full_TfIdf_DTM ???
        #txt_corpus <- tm_map(txt_corpus, content_transformer(stemDocument))        
        #print("StemmedWords:"); stemmed_words_Tf_df <- inspect_terms()
        #stemmed_stopped_Tf_df <- merge(stemmed_words_Tf_df, stopped_words_Tf_df, by="term", all=TRUE, suffixes=c(".stem", ".stop"))
        #myprint_df(stemmed_stopped_Tf_df)
        #print(subset(stemmed_stopped_Tf_df, grepl("compan", term)))
        #glb_corpus_lst[[txt_var]] <- txt_corpus
    }
    names(glb_corpus_lst) <- glb_txt_vars
        
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "extract.DTM"), major.inc=TRUE)

    glb_full_DTM_lst <- list(); glb_sprs_DTM_lst <- list();
    for (txt_var in glb_txt_vars) {
        print(sprintf("Extracting TfIDf terms for %s...", txt_var))        
        txt_corpus <- glb_corpus_lst[[txt_var]]
        
#         full_Tf_DTM <- DocumentTermMatrix(txt_corpus, 
#                                           control=list(weighting=weightTf))
        full_TfIdf_DTM <- DocumentTermMatrix(txt_corpus, 
                                          control=list(weighting=weightTfIdf))
        sprs_TfIdf_DTM <- removeSparseTerms(full_TfIdf_DTM, 
                                            glb_sprs_thresholds[txt_var])
        
#         glb_full_DTM_lst[[txt_var]] <- full_Tf_DTM
#         glb_sprs_DTM_lst[[txt_var]] <- sprs_Tf_DTM
        glb_full_DTM_lst[[txt_var]] <- full_TfIdf_DTM
        glb_sprs_DTM_lst[[txt_var]] <- sprs_TfIdf_DTM
    }

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
            paste0("extract.features_", "report.DTM"), major.inc=TRUE)
    
    for (txt_var in glb_txt_vars) {
        print(sprintf("Reporting TfIDf terms for %s...", txt_var))        
        full_TfIdf_DTM <- glb_full_DTM_lst[[txt_var]]
        sprs_TfIdf_DTM <- glb_sprs_DTM_lst[[txt_var]]        

        print("   Full TermMatrix:"); print(full_TfIdf_DTM)
        full_TfIdf_mtrx <- as.matrix(full_TfIdf_DTM)
        rownames(full_TfIdf_mtrx) <- rownames(glb_allobs_df) # print undreadable otherwise
        full_TfIdf_vctr <- colSums(full_TfIdf_mtrx)
        names(full_TfIdf_vctr) <- dimnames(full_TfIdf_DTM)[[2]]
        #grep("scene", names(full_TfIdf_vctr), value=TRUE)
        #which.max(full_TfIdf_mtrx[, "yearlong"])
        full_TfIdf_df <- as.data.frame(full_TfIdf_vctr)
        names(full_TfIdf_df) <- "TfIdf.full"
        full_TfIdf_df$term <- rownames(full_TfIdf_df)
        full_TfIdf_df$freq.full <- colSums(full_TfIdf_mtrx != 0)
        full_TfIdf_df <- orderBy(~ -TfIdf.full, full_TfIdf_df)

        print("   Sparse TermMatrix:"); print(sprs_TfIdf_DTM)
        sprs_TfIdf_vctr <- colSums(as.matrix(sprs_TfIdf_DTM))
        names(sprs_TfIdf_vctr) <- dimnames(sprs_TfIdf_DTM)[[2]]
        sprs_TfIdf_df <- as.data.frame(sprs_TfIdf_vctr)
        names(sprs_TfIdf_df) <- "TfIdf.sprs"
        sprs_TfIdf_df$term <- rownames(sprs_TfIdf_df)
        sprs_TfIdf_df$freq.sprs <- colSums(as.matrix(sprs_TfIdf_DTM) != 0)        
        sprs_TfIdf_df <- orderBy(~ -TfIdf.sprs, sprs_TfIdf_df)
        
        terms_TfIdf_df <- merge(full_TfIdf_df, sprs_TfIdf_df, all.x=TRUE)
        terms_TfIdf_df$in.sprs <- !is.na(terms_TfIdf_df$freq.sprs)
        plt_TfIdf_df <- subset(terms_TfIdf_df, 
                               TfIdf.full >= min(terms_TfIdf_df$TfIdf.sprs, na.rm=TRUE))
        plt_TfIdf_df$label <- ""
        plt_TfIdf_df[is.na(plt_TfIdf_df$TfIdf.sprs), "label"] <- 
            plt_TfIdf_df[is.na(plt_TfIdf_df$TfIdf.sprs), "term"]
        glb_important_terms[[txt_var]] <- union(glb_important_terms[[txt_var]],
            plt_TfIdf_df[is.na(plt_TfIdf_df$TfIdf.sprs), "term"])
        print(myplot_scatter(plt_TfIdf_df, "freq.full", "TfIdf.full", 
                             colorcol_name="in.sprs") + 
                  geom_text(aes(label=label), color="Black", size=3.5))
        
        melt_TfIdf_df <- orderBy(~ -value, melt(terms_TfIdf_df, id.var="term"))
        print(ggplot(melt_TfIdf_df, aes(value, color=variable)) + stat_ecdf() + 
                  geom_hline(yintercept=glb_sprs_thresholds[txt_var], 
                             linetype = "dotted"))
        
        melt_TfIdf_df <- orderBy(~ -value, 
                        melt(subset(terms_TfIdf_df, !is.na(TfIdf.sprs)), id.var="term"))
        print(myplot_hbar(melt_TfIdf_df, "term", "value", 
                          colorcol_name="variable"))
        
        melt_TfIdf_df <- orderBy(~ -value, 
                        melt(subset(terms_TfIdf_df, is.na(TfIdf.sprs)), id.var="term"))
        print(myplot_hbar(head(melt_TfIdf_df, 10), "term", "value", 
                          colorcol_name="variable"))
    }

#     sav_full_DTM_lst <- glb_full_DTM_lst
#     sav_sprs_DTM_lst <- glb_sprs_DTM_lst
#     print(identical(sav_glb_corpus_lst, glb_corpus_lst))
#     print(all.equal(length(sav_glb_corpus_lst), length(glb_corpus_lst)))
#     print(all.equal(names(sav_glb_corpus_lst), names(glb_corpus_lst)))
#     print(all.equal(sav_glb_corpus_lst[["Headline"]], glb_corpus_lst[["Headline"]]))

#     print(identical(sav_full_DTM_lst, glb_full_DTM_lst))
#     print(identical(sav_sprs_DTM_lst, glb_sprs_DTM_lst))
        
    rm(full_TfIdf_mtrx, full_TfIdf_df, melt_TfIdf_df, terms_TfIdf_df)

    # Create txt features
    if ((length(glb_txt_vars) > 1) &&
        (length(unique(pfxs <- sapply(glb_txt_vars, 
                    function(txt) toupper(substr(txt, 1, 1))))) < length(glb_txt_vars)))
            stop("Prefixes for corpus freq terms not unique: ", pfxs)
    
    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
                            paste0("extract.features_", "bind.DTM"), 
                                         major.inc=TRUE)
    for (txt_var in glb_txt_vars) {
        print(sprintf("Binding DTM for %s...", txt_var))
        txt_var_pfx <- toupper(substr(txt_var, 1, 1))
        txt_X_df <- as.data.frame(as.matrix(glb_sprs_DTM_lst[[txt_var]]))
        colnames(txt_X_df) <- paste(txt_var_pfx, ".T.",
                                    make.names(colnames(txt_X_df)), sep="")
        rownames(txt_X_df) <- rownames(glb_allobs_df) # warning otherwise
#         plt_X_df <- cbind(txt_X_df, glb_allobs_df[, c(glb_id_var, glb_rsp_var)])
#         print(myplot_box(df=plt_X_df, ycol_names="H.T.today", xcol_name=glb_rsp_var))

#         log_X_df <- log(1 + txt_X_df)
#         colnames(log_X_df) <- paste(colnames(txt_X_df), ".log", sep="")
#         plt_X_df <- cbind(log_X_df, glb_allobs_df[, c(glb_id_var, glb_rsp_var)])
#         print(myplot_box(df=plt_X_df, ycol_names="H.T.today.log", xcol_name=glb_rsp_var))
        glb_allobs_df <- cbind(glb_allobs_df, txt_X_df) # TfIdf is normalized
        #glb_allobs_df <- cbind(glb_allobs_df, log_X_df) # if using non-normalized metrics 
    }
    #identical(chk_entity_df, glb_allobs_df)
    #chk_entity_df <- glb_allobs_df

    extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, 
                            paste0("extract.features_", "bind.DXM"), 
                                         major.inc=TRUE)

#sav_allobs_df <- glb_allobs_df
    glb_punct_vctr <- c("!", "\"", "#", "\\$", "%", "&", "'", 
                        "\\(|\\)",# "\\(", "\\)", 
                        "\\*", "\\+", ",", "-", "\\.", "/", ":", ";", 
                        "<|>", # "<", 
                        "=", 
                        # ">", 
                        "\\?", "@", "\\[", "\\\\", "\\]", "^", "_", "`", 
                        "\\{", "\\|", "\\}", "~")
    txt_X_df <- glb_allobs_df[, c(glb_id_var, ".rnorm"), FALSE]
    txt_X_df <- foreach(txt_var=glb_txt_vars, .combine=cbind) %dopar% {   
    #for (txt_var in glb_txt_vars) {
        print(sprintf("Binding DXM for %s...", txt_var))
        txt_var_pfx <- toupper(substr(txt_var, 1, 1))        
        #txt_X_df <- glb_allobs_df[, c(glb_id_var, ".rnorm"), FALSE]
        
        txt_full_DTM_mtrx <- as.matrix(glb_full_DTM_lst[[txt_var]])
        rownames(txt_full_DTM_mtrx) <- rownames(glb_allobs_df) # print undreadable otherwise
        #print(txt_full_DTM_mtrx[txt_full_DTM_mtrx[, "ebola"] != 0, "ebola"])
        
        # Create <txt_var>.T.<term> for glb_important_terms
        for (term in glb_important_terms[[txt_var]])
            txt_X_df[, paste0(txt_var_pfx, ".T.", make.names(term))] <- 
                txt_full_DTM_mtrx[, term]
                
        # Create <txt_var>.nwrds.log & .nwrds.unq.log
        txt_X_df[, paste0(txt_var_pfx, ".nwrds.log")] <- 
            log(1 + mycount_pattern_occ("\\w+", glb_txt_lst[[txt_var]]))
        txt_X_df[, paste0(txt_var_pfx, ".nwrds.unq.log")] <- 
            log(1 + rowSums(txt_full_DTM_mtrx != 0))
        txt_X_df[, paste0(txt_var_pfx, ".sum.TfIdf")] <- 
            rowSums(txt_full_DTM_mtrx) 
        txt_X_df[, paste0(txt_var_pfx, ".ratio.sum.TfIdf.nwrds")] <- 
            txt_X_df[, paste0(txt_var_pfx, ".sum.TfIdf")] / 
            (exp(txt_X_df[, paste0(txt_var_pfx, ".nwrds.log")]) - 1)
        txt_X_df[is.nan(txt_X_df[, paste0(txt_var_pfx, ".ratio.sum.TfIdf.nwrds")]),
                 paste0(txt_var_pfx, ".ratio.sum.TfIdf.nwrds")] <- 0

        # Create <txt_var>.nchrs.log
        txt_X_df[, paste0(txt_var_pfx, ".nchrs.log")] <- 
            log(1 + mycount_pattern_occ(".", glb_allobs_df[, txt_var]))
        txt_X_df[, paste0(txt_var_pfx, ".nuppr.log")] <- 
            log(1 + mycount_pattern_occ("[[:upper:]]", glb_allobs_df[, txt_var]))
        txt_X_df[, paste0(txt_var_pfx, ".ndgts.log")] <- 
            log(1 + mycount_pattern_occ("[[:digit:]]", glb_allobs_df[, txt_var]))

        # Create <txt_var>.npnct?.log
        # would this be faster if it's iterated over each row instead of 
        #   each created column ???
        for (punct_ix in 1:length(glb_punct_vctr)) { 
#             smp0 <- " "
#             smp1 <- "! \" # $ % & ' ( ) * + , - . / : ; < = > ? @ [ \ ] ^ _ ` { | } ~"
#             smp2 <- paste(smp1, smp1, sep=" ")
#             print(sprintf("Testing %s pattern:", glb_punct_vctr[punct_ix])) 
#             results <- mycount_pattern_occ(glb_punct_vctr[punct_ix], c(smp0, smp1, smp2))
#             names(results) <- NULL; print(results)
            txt_X_df[, 
                paste0(txt_var_pfx, ".npnct", sprintf("%02d", punct_ix), ".log")] <-
                log(1 + mycount_pattern_occ(glb_punct_vctr[punct_ix], 
                                            glb_allobs_df[, txt_var]))
        }
#         print(head(glb_allobs_df[glb_allobs_df[, "A.npnct23.log"] > 0, 
#                                     c("UniqueID", "Popular", "Abstract", "A.npnct23.log")]))    
        
        # Create <txt_var>.nstopwrds.log & <txt_var>ratio.nstopwrds.nwrds
        stop_words_rex_str <- paste0("\\b(", paste0(c(glb_append_stop_words[[txt_var]], 
                                       stopwords("english")), collapse="|"),
                                     ")\\b")
        txt_X_df[, paste0(txt_var_pfx, ".nstopwrds", ".log")] <-
            log(1 + mycount_pattern_occ(stop_words_rex_str, glb_txt_lst[[txt_var]]))
        txt_X_df[, paste0(txt_var_pfx, ".ratio.nstopwrds.nwrds")] <-
            exp(txt_X_df[, paste0(txt_var_pfx, ".nstopwrds", ".log")] - 
                txt_X_df[, paste0(txt_var_pfx, ".nwrds", ".log")])

        # Create <txt_var>.P.http
        txt_X_df[, paste(txt_var_pfx, ".P.http", sep="")] <- 
            as.integer(0 + mycount_pattern_occ("http", glb_allobs_df[, txt_var]))    
    
        txt_X_df <- subset(txt_X_df, select=-.rnorm)
        txt_X_df <- txt_X_df[, -grep(glb_id_var, names(txt_X_df), fixed=TRUE), FALSE]
        #glb_allobs_df <- cbind(glb_allobs_df, txt_X_df)
    }
    glb_allobs_df <- cbind(glb_allobs_df, txt_X_df)
    #myplot_box(glb_allobs_df, "A.sum.TfIdf", glb_rsp_var)

    # Generate summaries
#     print(summary(glb_allobs_df))
#     print(sapply(names(glb_allobs_df), function(col) sum(is.na(glb_allobs_df[, col]))))
#     print(summary(glb_trnobs_df))
#     print(sapply(names(glb_trnobs_df), function(col) sum(is.na(glb_trnobs_df[, col]))))
#     print(summary(glb_newobs_df))
#     print(sapply(names(glb_newobs_df), function(col) sum(is.na(glb_newobs_df[, col]))))

    glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
                                          glb_txt_vars)
    rm(log_X_df, txt_X_df)
}

# print(sapply(names(glb_trnobs_df), function(col) sum(is.na(glb_trnobs_df[, col]))))
# print(sapply(names(glb_newobs_df), function(col) sum(is.na(glb_newobs_df[, col]))))

# print(myplot_scatter(glb_trnobs_df, "<col1_name>", "<col2_name>", smooth=TRUE))

rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr, 
   glb_full_DTM_lst, glb_sprs_DTM_lst, txt_corpus, txt_vctr)
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'corpus_lst' not found
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'full_TfIdf_DTM' not found
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'full_TfIdf_vctr' not found
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'glb_full_DTM_lst' not found
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'glb_sprs_DTM_lst' not found
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'txt_corpus' not found
```

```
## Warning in rm(corpus_lst, full_TfIdf_DTM, full_TfIdf_vctr,
## glb_full_DTM_lst, : object 'txt_vctr' not found
```

```r
extract.features_chunk_df <- myadd_chunk(extract.features_chunk_df, "extract.features_end", 
                                     major.inc=TRUE)
```

```
##                                 label step_major step_minor     bgn
## 2 extract.features_factorize.str.vars          2          0 108.047
## 3                extract.features_end          3          0 108.204
##       end elapsed
## 2 108.203   0.157
## 3      NA      NA
```

```r
myplt_chunk(extract.features_chunk_df)
```

```
##                                 label step_major step_minor     bgn
## 1                extract.features_bgn          1          0  10.230
## 2 extract.features_factorize.str.vars          2          0 108.047
##       end elapsed duration
## 1 108.046  97.816   97.816
## 2 108.203   0.157    0.156
## [1] "Total Elapsed Time: 108.203 secs"
```

![](Google_Flu_template2_files/figure-html/extract.features-13.png) 

```r
# if (glb_save_envir)
#     save(glb_feats_df, 
#          glb_allobs_df, #glb_trnobs_df, glb_fitobs_df, glb_OOBobs_df, glb_newobs_df,
#          file=paste0(glb_out_pfx, "extract_features_dsk.RData"))
# load(paste0(glb_out_pfx, "extract_features_dsk.RData"))

replay.petrisim(pn=glb_analytics_pn, 
    replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
        "data.training.all","data.new")), flip_coord=TRUE)
```

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0
```

![](Google_Flu_template2_files/figure-html/extract.features-14.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "cluster.data", major.inc=TRUE)
```

```
##              label step_major step_minor     bgn     end elapsed
## 5 extract.features          3          0  10.224 109.764   99.54
## 6     cluster.data          4          0 109.765      NA      NA
```

### Step `4.0: cluster data`

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "manage.missing.data", major.inc=FALSE)
```

```
##                 label step_major step_minor     bgn     end elapsed
## 6        cluster.data          4          0 109.765 127.266  17.501
## 7 manage.missing.data          4          1 127.267      NA      NA
```

```r
# print(sapply(names(glb_trnobs_df), function(col) sum(is.na(glb_trnobs_df[, col]))))
# print(sapply(names(glb_newobs_df), function(col) sum(is.na(glb_newobs_df[, col]))))
# glb_trnobs_df <- na.omit(glb_trnobs_df)
# glb_newobs_df <- na.omit(glb_newobs_df)
# df[is.na(df)] <- 0

mycheck_problem_data(glb_allobs_df)
```

```
## [1] "numeric data missing in : "
##     ILI.2.lag ILI.2.lag.log 
##             2             2 
## [1] "numeric data w/ 0s in : "
##  Week.bgn.wkday.fctr   Week.bgn.hour.fctr Week.bgn.minute.fctr 
##                  469                  469                  469 
## Week.bgn.second.fctr  Week.end.wkday.fctr   Week.end.hour.fctr 
##                  469                  469                  469 
## Week.end.minute.fctr Week.end.second.fctr   Week.bgn.last1.log 
##                  469                  469                    1 
##   Week.bgn.last2.log  Week.bgn.last10.log Week.bgn.last100.log 
##                    2                   10                  100 
##   Week.end.last1.log   Week.end.last2.log  Week.end.last10.log 
##                    1                    2                   10 
## Week.end.last100.log 
##                  100 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##     Week Week.bgn Week.end 
##        0        0        0
```

```r
# glb_allobs_df <- na.omit(glb_allobs_df)

# Not refactored into mydsutils.R since glb_*_df might be reassigned
glb_impute_missing_data <- function() {
    
    require(mice)
    set.seed(glb_mice_complete.seed)
    inp_impent_df <- glb_allobs_df[, setdiff(names(glb_allobs_df), 
                                union(glb_exclude_vars_as_features, glb_rsp_var))]
    print("Summary before imputation: ")
    print(summary(inp_impent_df))
    out_impent_df <- complete(mice(inp_impent_df))
    print(summary(out_impent_df))
    
    # complete(mice()) changes attributes of factors even though values don't change
    ret_vars <- sapply(names(out_impent_df), 
                       function(col) ifelse(!identical(out_impent_df[, col], inp_impent_df[, col]), 
                                            col, ""))
    ret_vars <- ret_vars[ret_vars != ""]
    return(out_impent_df[, ret_vars])
}

if (glb_impute_na_data && 
    (length(myfind_numerics_missing(glb_allobs_df)) > 0) &&
    (ncol(nonna_df <- glb_impute_missing_data()) > 0)) {
    for (col in names(nonna_df)) {
        glb_allobs_df[, paste0(col, ".nonNA")] <- nonna_df[, col]
        glb_exclude_vars_as_features <- c(glb_exclude_vars_as_features, col)        
    }
}    
```

```
##     ILI.2.lag ILI.2.lag.log 
##             2             2
```

```
## Loading required package: mice
## Loading required package: Rcpp
## mice 2.22 2014-06-10
```

```
## [1] "Summary before imputation: "
##     Queries            .rnorm         ILI.2.lag.log     Week.bgn.year.fctr
##  Min.   :0.04117   Min.   :-2.63119   Min.   :-0.6272   2006   : 53       
##  1st Qu.:0.17928   1st Qu.:-0.64215   1st Qu.:-0.0770   2004   : 52       
##  Median :0.28685   Median : 0.03471   Median : 0.2403   2005   : 52       
##  Mean   :0.29970   Mean   :-0.01675   Mean   : 0.3494   2007   : 52       
##  3rd Qu.:0.39177   3rd Qu.: 0.67690   3rd Qu.: 0.7023   2008   : 52       
##  Max.   :1.00000   Max.   : 2.59377   Max.   : 2.0306   2009   : 52       
##                                       NA's   :2         (Other):156       
##  Week.bgn.month.fctr Week.bgn.date.fctr Week.bgn.wkday.fctr Week.bgn.wkend
##  01     : 41         (0.97,7]:108       0:469               Min.   :1     
##  05     : 41         (7,13]  : 94                           1st Qu.:1     
##  07     : 41         (13,19] : 92                           Median :1     
##  10     : 41         (19,25] : 93                           Mean   :1     
##  08     : 40         (25,31] : 82                           3rd Qu.:1     
##  04     : 39                                                Max.   :1     
##  (Other):226                                                              
##  Week.bgn.hour.fctr Week.bgn.minute.fctr Week.bgn.second.fctr
##  Min.   :0          Min.   :0            Min.   :0           
##  1st Qu.:0          1st Qu.:0            1st Qu.:0           
##  Median :0          Median :0            Median :0           
##  Mean   :0          Mean   :0            Mean   :0           
##  3rd Qu.:0          3rd Qu.:0            3rd Qu.:0           
##  Max.   :0          Max.   :0            Max.   :0           
##                                                              
##  Week.end.year.fctr Week.end.month.fctr Week.end.date.fctr
##  2006   : 53        01     : 41         (0.97,7]:108      
##  2004   : 52        05     : 41         (7,13]  : 94      
##  2005   : 52        07     : 41         (13,19] : 92      
##  2007   : 52        10     : 41         (19,25] : 93      
##  2008   : 52        08     : 40         (25,31] : 82      
##  2009   : 52        04     : 39                           
##  (Other):156        (Other):226                           
##  Week.end.wkday.fctr Week.end.wkend Week.end.hour.fctr
##  0:469               Min.   :1      Min.   :0         
##                      1st Qu.:1      1st Qu.:0         
##                      Median :1      Median :0         
##                      Mean   :1      Mean   :0         
##                      3rd Qu.:1      3rd Qu.:0         
##                      Max.   :1      Max.   :0         
##                                                       
##  Week.end.minute.fctr Week.end.second.fctr Week.bgn.last1.log
##  Min.   :0            Min.   :0            Min.   : 0.00     
##  1st Qu.:0            1st Qu.:0            1st Qu.:13.31     
##  Median :0            Median :0            Median :13.31     
##  Mean   :0            Mean   :0            Mean   :13.28     
##  3rd Qu.:0            3rd Qu.:0            3rd Qu.:13.31     
##  Max.   :0            Max.   :0            Max.   :13.32     
##                                                              
##  Week.bgn.last2.log Week.bgn.last10.log Week.bgn.last100.log
##  Min.   : 0.00      Min.   : 0.00       Min.   : 0.00       
##  1st Qu.:14.01      1st Qu.:15.62       1st Qu.:17.92       
##  Median :14.01      Median :15.62       Median :17.92       
##  Mean   :13.95      Mean   :15.28       Mean   :14.10       
##  3rd Qu.:14.01      3rd Qu.:15.62       3rd Qu.:17.92       
##  Max.   :14.01      Max.   :15.62       Max.   :17.92       
##                                                             
##  Week.end.last1.log Week.end.last2.log Week.end.last10.log
##  Min.   : 0.00      Min.   : 0.00      Min.   : 0.00      
##  1st Qu.:13.31      1st Qu.:14.01      1st Qu.:15.62      
##  Median :13.31      Median :14.01      Median :15.62      
##  Mean   :13.28      Mean   :13.95      Mean   :15.28      
##  3rd Qu.:13.31      3rd Qu.:14.01      3rd Qu.:15.62      
##  Max.   :13.32      Max.   :14.01      Max.   :15.62      
##                                                           
##  Week.end.last100.log
##  Min.   : 0.00       
##  1st Qu.:17.92       
##  Median :17.92       
##  Mean   :14.10       
##  3rd Qu.:17.92       
##  Max.   :17.92       
##                      
## 
##  iter imp variable
##   1   1  ILI.2.lag.log
##   1   2  ILI.2.lag.log
##   1   3  ILI.2.lag.log
##   1   4  ILI.2.lag.log
##   1   5  ILI.2.lag.log
##   2   1  ILI.2.lag.log
##   2   2  ILI.2.lag.log
##   2   3  ILI.2.lag.log
##   2   4  ILI.2.lag.log
##   2   5  ILI.2.lag.log
##   3   1  ILI.2.lag.log
##   3   2  ILI.2.lag.log
##   3   3  ILI.2.lag.log
##   3   4  ILI.2.lag.log
##   3   5  ILI.2.lag.log
##   4   1  ILI.2.lag.log
##   4   2  ILI.2.lag.log
##   4   3  ILI.2.lag.log
##   4   4  ILI.2.lag.log
##   4   5  ILI.2.lag.log
##   5   1  ILI.2.lag.log
##   5   2  ILI.2.lag.log
##   5   3  ILI.2.lag.log
##   5   4  ILI.2.lag.log
##   5   5  ILI.2.lag.log
##     Queries            .rnorm         ILI.2.lag.log     
##  Min.   :0.04117   Min.   :-2.63119   Min.   :-0.62719  
##  1st Qu.:0.17928   1st Qu.:-0.64215   1st Qu.:-0.07625  
##  Median :0.28685   Median : 0.03471   Median : 0.24476  
##  Mean   :0.29970   Mean   :-0.01675   Mean   : 0.34952  
##  3rd Qu.:0.39177   3rd Qu.: 0.67690   3rd Qu.: 0.70151  
##  Max.   :1.00000   Max.   : 2.59377   Max.   : 2.03063  
##                                                         
##  Week.bgn.year.fctr Week.bgn.month.fctr Week.bgn.date.fctr
##  2006   : 53        01     : 41         (0.97,7]:108      
##  2004   : 52        05     : 41         (7,13]  : 94      
##  2005   : 52        07     : 41         (13,19] : 92      
##  2007   : 52        10     : 41         (19,25] : 93      
##  2008   : 52        08     : 40         (25,31] : 82      
##  2009   : 52        04     : 39                           
##  (Other):156        (Other):226                           
##  Week.bgn.wkday.fctr Week.bgn.wkend Week.bgn.hour.fctr
##  0:469               Min.   :1      Min.   :0         
##                      1st Qu.:1      1st Qu.:0         
##                      Median :1      Median :0         
##                      Mean   :1      Mean   :0         
##                      3rd Qu.:1      3rd Qu.:0         
##                      Max.   :1      Max.   :0         
##                                                       
##  Week.bgn.minute.fctr Week.bgn.second.fctr Week.end.year.fctr
##  Min.   :0            Min.   :0            2006   : 53       
##  1st Qu.:0            1st Qu.:0            2004   : 52       
##  Median :0            Median :0            2005   : 52       
##  Mean   :0            Mean   :0            2007   : 52       
##  3rd Qu.:0            3rd Qu.:0            2008   : 52       
##  Max.   :0            Max.   :0            2009   : 52       
##                                            (Other):156       
##  Week.end.month.fctr Week.end.date.fctr Week.end.wkday.fctr Week.end.wkend
##  01     : 41         (0.97,7]:108       0:469               Min.   :1     
##  05     : 41         (7,13]  : 94                           1st Qu.:1     
##  07     : 41         (13,19] : 92                           Median :1     
##  10     : 41         (19,25] : 93                           Mean   :1     
##  08     : 40         (25,31] : 82                           3rd Qu.:1     
##  04     : 39                                                Max.   :1     
##  (Other):226                                                              
##  Week.end.hour.fctr Week.end.minute.fctr Week.end.second.fctr
##  Min.   :0          Min.   :0            Min.   :0           
##  1st Qu.:0          1st Qu.:0            1st Qu.:0           
##  Median :0          Median :0            Median :0           
##  Mean   :0          Mean   :0            Mean   :0           
##  3rd Qu.:0          3rd Qu.:0            3rd Qu.:0           
##  Max.   :0          Max.   :0            Max.   :0           
##                                                              
##  Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
##  Min.   : 0.00      Min.   : 0.00      Min.   : 0.00      
##  1st Qu.:13.31      1st Qu.:14.01      1st Qu.:15.62      
##  Median :13.31      Median :14.01      Median :15.62      
##  Mean   :13.28      Mean   :13.95      Mean   :15.28      
##  3rd Qu.:13.31      3rd Qu.:14.01      3rd Qu.:15.62      
##  Max.   :13.32      Max.   :14.01      Max.   :15.62      
##                                                           
##  Week.bgn.last100.log Week.end.last1.log Week.end.last2.log
##  Min.   : 0.00        Min.   : 0.00      Min.   : 0.00     
##  1st Qu.:17.92        1st Qu.:13.31      1st Qu.:14.01     
##  Median :17.92        Median :13.31      Median :14.01     
##  Mean   :14.10        Mean   :13.28      Mean   :13.95     
##  3rd Qu.:17.92        3rd Qu.:13.31      3rd Qu.:14.01     
##  Max.   :17.92        Max.   :13.32      Max.   :14.01     
##                                                            
##  Week.end.last10.log Week.end.last100.log
##  Min.   : 0.00       Min.   : 0.00       
##  1st Qu.:15.62       1st Qu.:17.92       
##  Median :15.62       Median :17.92       
##  Mean   :15.28       Mean   :14.10       
##  3rd Qu.:15.62       3rd Qu.:17.92       
##  Max.   :15.62       Max.   :17.92       
## 
```

```r
mycheck_problem_data(glb_allobs_df, terminate = TRUE)
```

```
## [1] "numeric data missing in : "
##     ILI.2.lag ILI.2.lag.log 
##             2             2 
## [1] "numeric data w/ 0s in : "
##  Week.bgn.wkday.fctr   Week.bgn.hour.fctr Week.bgn.minute.fctr 
##                  469                  469                  469 
## Week.bgn.second.fctr  Week.end.wkday.fctr   Week.end.hour.fctr 
##                  469                  469                  469 
## Week.end.minute.fctr Week.end.second.fctr   Week.bgn.last1.log 
##                  469                  469                    1 
##   Week.bgn.last2.log  Week.bgn.last10.log Week.bgn.last100.log 
##                    2                   10                  100 
##   Week.end.last1.log   Week.end.last2.log  Week.end.last10.log 
##                    1                    2                   10 
## Week.end.last100.log 
##                  100 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##     Week Week.bgn Week.end 
##        0        0        0
```

## Step `4.1: manage missing data`

```r
if (glb_cluster) {
    require(proxy)
    #require(hash)
    require(dynamicTreeCut)

#     glb_hash <- hash(key=unique(glb_allobs_df$myCategory), 
#                      values=1:length(unique(glb_allobs_df$myCategory)))
#     glb_hash_lst <- hash(key=unique(glb_allobs_df$myCategory), 
#                      values=1:length(unique(glb_allobs_df$myCategory)))
#stophere; sav_allobs_df <- glb_allobs_df; 
    print("Clustering features: ")
    print(cluster_vars <- grep("[HSA]\\.[PT]\\.", names(glb_allobs_df), value=TRUE))
    #print(cluster_vars <- grep("[HSA]\\.", names(glb_allobs_df), value=TRUE))
    glb_allobs_df$.clusterid <- 1    
    #print(max(table(glb_allobs_df$myCategory.fctr) / 20))
    for (myCategory in c("##", "Business#Business Day#Dealbook", "OpEd#Opinion#", 
                         "Styles#U.S.#", "Business#Technology#", "Science#Health#",
                         "Culture#Arts#")) {
        ctgry_allobs_df <- glb_allobs_df[glb_allobs_df$myCategory == myCategory, ]
        
        dstns_dist <- dist(ctgry_allobs_df[, cluster_vars], method = "cosine")
        dstns_mtrx <- as.matrix(dstns_dist)
        print(sprintf("max distance(%0.4f) pair:", max(dstns_mtrx)))
        row_ix <- ceiling(which.max(dstns_mtrx) / ncol(dstns_mtrx))
        col_ix <- which.max(dstns_mtrx[row_ix, ])
        print(ctgry_allobs_df[c(row_ix, col_ix), 
            c("UniqueID", "Popular", "myCategory", "Headline", cluster_vars)])
    
        min_dstns_mtrx <- dstns_mtrx
        diag(min_dstns_mtrx) <- 1
        print(sprintf("min distance(%0.4f) pair:", min(min_dstns_mtrx)))
        row_ix <- ceiling(which.min(min_dstns_mtrx) / ncol(min_dstns_mtrx))
        col_ix <- which.min(min_dstns_mtrx[row_ix, ])
        print(ctgry_allobs_df[c(row_ix, col_ix), 
            c("UniqueID", "Popular", "myCategory", "Headline", cluster_vars)])                          
    
        clusters <- hclust(dstns_dist, method = "ward.D2")
        #plot(clusters, labels=NULL, hang=-1)
        myplclust(clusters, lab.col=unclass(ctgry_allobs_df[, glb_rsp_var]))
        
        #clusterGroups = cutree(clusters, k=7)
        clusterGroups <- cutreeDynamic(clusters, minClusterSize=20, method="tree", deepSplit=0)
        # Unassigned groups are labeled 0; the largest group has label 1
        table(clusterGroups, ctgry_allobs_df[, glb_rsp_var], useNA="ifany")   
        #print(ctgry_allobs_df[which(clusterGroups == 1), c("UniqueID", "Popular", "Headline")])
        #print(ctgry_allobs_df[(clusterGroups == 1) & !is.na(ctgry_allobs_df$Popular) & (ctgry_allobs_df$Popular==1), c("UniqueID", "Popular", "Headline")])
        clusterGroups[clusterGroups == 0] <- 1
        table(clusterGroups, ctgry_allobs_df[, glb_rsp_var], useNA="ifany")        
        #summary(factor(clusterGroups))
#         clusterGroups <- clusterGroups + 
#                 100 * # has to be > max(table(glb_allobs_df$myCategory.fctr) / minClusterSize=20)
#                             which(levels(glb_allobs_df$myCategory.fctr) == myCategory)
#         table(clusterGroups, ctgry_allobs_df[, glb_rsp_var], useNA="ifany")        
    
        # add to glb_allobs_df - then split the data again
        glb_allobs_df[glb_allobs_df$myCategory==myCategory,]$.clusterid <- clusterGroups
        #print(unique(glb_allobs_df$.clusterid))
        #print(glb_feats_df[glb_feats_df$id == ".clusterid.fctr", ])
    }
    
    ctgry_xtab_df <- orderBy(reformulate(c("-", ".n")),
                              mycreate_sqlxtab_df(glb_allobs_df,
        c("myCategory", ".clusterid", glb_rsp_var)))
    ctgry_cast_df <- orderBy(~ -Y -NA, dcast(ctgry_xtab_df, 
                           myCategory + .clusterid ~ 
                               Popular.fctr, sum, value.var=".n"))
    print(ctgry_cast_df)
    #print(orderBy(~ myCategory -Y -NA, ctgry_cast_df))
    # write.table(ctgry_cast_df, paste0(glb_out_pfx, "ctgry_clst.csv"), 
    #             row.names=FALSE)
    
    print(ctgry_sum_tbl <- table(glb_allobs_df$myCategory, glb_allobs_df$.clusterid, 
                                 glb_allobs_df[, glb_rsp_var], 
                                 useNA="ifany"))
#     dsp_obs(.clusterid=1, myCategory="OpEd#Opinion#", 
#             cols=c("UniqueID", "Popular", "myCategory", ".clusterid", "Headline"),
#             all=TRUE)
    
    glb_allobs_df$.clusterid.fctr <- as.factor(glb_allobs_df$.clusterid)
    glb_exclude_vars_as_features <- c(glb_exclude_vars_as_features, 
                                      ".clusterid")
    glb_interaction_only_features["myCategory.fctr"] <- c(".clusterid.fctr")
    glb_exclude_vars_as_features <- c(glb_exclude_vars_as_features, 
                                      cluster_vars)
}

# Re-partition
glb_trnobs_df <- subset(glb_allobs_df, .src == "Train")
glb_newobs_df <- subset(glb_allobs_df, .src == "Test")

glb_chunks_df <- myadd_chunk(glb_chunks_df, "select.features", major.inc=TRUE)
```

```
##                 label step_major step_minor     bgn     end elapsed
## 7 manage.missing.data          4          1 127.267 128.257    0.99
## 8     select.features          5          0 128.258      NA      NA
```

## Step `5.0: select features`

```r
print(glb_feats_df <- myselect_features(entity_df=glb_trnobs_df, 
                       exclude_vars_as_features=glb_exclude_vars_as_features, 
                       rsp_var=glb_rsp_var))
```

```
## Warning in cor(data.matrix(entity_df[, sel_feats]), y =
## as.numeric(entity_df[, : the standard deviation is zero
```

```
##                                                  id         cor.y
## ILI                                             ILI  0.9451682372
## ILI.2.lag.log                         ILI.2.lag.log  0.9214355113
## ILI.2.lag.log.nonNA             ILI.2.lag.log.nonNA  0.9204139777
## ILI.2.lag                                 ILI.2.lag  0.8591322165
## Queries                                     Queries  0.8420332864
## Week.bgn.last100.log           Week.bgn.last100.log  0.2000557734
## Week.end.last100.log           Week.end.last100.log  0.2000557734
## Week.bgn.year.fctr               Week.bgn.year.fctr  0.1987250820
## Week.end.year.fctr               Week.end.year.fctr  0.1987250820
## Week.bgn.year.fctr.nonNA   Week.bgn.year.fctr.nonNA  0.1987250820
## Week.bgn.month.fctr             Week.bgn.month.fctr -0.1979117680
## Week.end.month.fctr             Week.end.month.fctr -0.1979117680
## Week.bgn.month.fctr.nonNA Week.bgn.month.fctr.nonNA -0.1979117680
## Week.bgn.POSIX                       Week.bgn.POSIX  0.1730840769
## Week.end.POSIX                       Week.end.POSIX  0.1730840769
## Week.bgn.zoo                           Week.bgn.zoo  0.1730840769
## Week.end.zoo                           Week.end.zoo  0.1730840769
## Week.bgn.last2.log               Week.bgn.last2.log -0.0488875907
## Week.end.last2.log               Week.end.last2.log -0.0488875907
## Week.bgn.last1.log               Week.bgn.last1.log -0.0473819651
## Week.end.last1.log               Week.end.last1.log -0.0473819651
## .rnorm                                       .rnorm  0.0400349990
## Week.bgn.date.fctr               Week.bgn.date.fctr  0.0088237363
## Week.end.date.fctr               Week.end.date.fctr  0.0088237363
## Week.bgn.date.fctr.nonNA   Week.bgn.date.fctr.nonNA  0.0088237363
## Week.bgn.last10.log             Week.bgn.last10.log  0.0006098826
## Week.end.last10.log             Week.end.last10.log  0.0006098826
## Week.bgn.wkend                       Week.bgn.wkend            NA
## Week.bgn.hour.fctr               Week.bgn.hour.fctr            NA
## Week.bgn.minute.fctr           Week.bgn.minute.fctr            NA
## Week.bgn.second.fctr           Week.bgn.second.fctr            NA
## Week.end.wkend                       Week.end.wkend            NA
## Week.end.hour.fctr               Week.end.hour.fctr            NA
## Week.end.minute.fctr           Week.end.minute.fctr            NA
## Week.end.second.fctr           Week.end.second.fctr            NA
## Week.bgn.wkday.fctr             Week.bgn.wkday.fctr            NA
## Week.end.wkday.fctr             Week.end.wkday.fctr            NA
##                           exclude.as.feat    cor.y.abs
## ILI                                     1 0.9451682372
## ILI.2.lag.log                           1 0.9214355113
## ILI.2.lag.log.nonNA                     0 0.9204139777
## ILI.2.lag                               1 0.8591322165
## Queries                                 0 0.8420332864
## Week.bgn.last100.log                    0 0.2000557734
## Week.end.last100.log                    0 0.2000557734
## Week.bgn.year.fctr                      1 0.1987250820
## Week.end.year.fctr                      0 0.1987250820
## Week.bgn.year.fctr.nonNA                0 0.1987250820
## Week.bgn.month.fctr                     1 0.1979117680
## Week.end.month.fctr                     0 0.1979117680
## Week.bgn.month.fctr.nonNA               0 0.1979117680
## Week.bgn.POSIX                          1 0.1730840769
## Week.end.POSIX                          1 0.1730840769
## Week.bgn.zoo                            1 0.1730840769
## Week.end.zoo                            1 0.1730840769
## Week.bgn.last2.log                      0 0.0488875907
## Week.end.last2.log                      0 0.0488875907
## Week.bgn.last1.log                      0 0.0473819651
## Week.end.last1.log                      0 0.0473819651
## .rnorm                                  0 0.0400349990
## Week.bgn.date.fctr                      1 0.0088237363
## Week.end.date.fctr                      0 0.0088237363
## Week.bgn.date.fctr.nonNA                0 0.0088237363
## Week.bgn.last10.log                     0 0.0006098826
## Week.end.last10.log                     0 0.0006098826
## Week.bgn.wkend                          0           NA
## Week.bgn.hour.fctr                      0           NA
## Week.bgn.minute.fctr                    0           NA
## Week.bgn.second.fctr                    0           NA
## Week.end.wkend                          0           NA
## Week.end.hour.fctr                      0           NA
## Week.end.minute.fctr                    0           NA
## Week.end.second.fctr                    0           NA
## Week.bgn.wkday.fctr                     0           NA
## Week.end.wkday.fctr                     0           NA
```

```r
# sav_feats_df <- glb_feats_df; glb_feats_df <- sav_feats_df
print(glb_feats_df <- orderBy(~-cor.y, 
          myfind_cor_features(feats_df=glb_feats_df, obs_df=glb_trnobs_df, 
                              rsp_var=glb_rsp_var)))
```

```
## Loading required package: reshape2
```

```
## [1] "cor(Week.bgn.last1.log, Week.end.last1.log)=1.0000"
## [1] "cor(ILI.log, Week.bgn.last1.log)=-0.0474"
## [1] "cor(ILI.log, Week.end.last1.log)=-0.0474"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.end.last1.log as highly correlated with
## Week.bgn.last1.log
```

```
## [1] "cor(Week.bgn.last100.log, Week.end.last100.log)=1.0000"
## [1] "cor(ILI.log, Week.bgn.last100.log)=0.2001"
## [1] "cor(ILI.log, Week.end.last100.log)=0.2001"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.end.last100.log as highly correlated with
## Week.bgn.last100.log
```

```
## [1] "cor(Week.bgn.last2.log, Week.end.last2.log)=1.0000"
## [1] "cor(ILI.log, Week.bgn.last2.log)=-0.0489"
## [1] "cor(ILI.log, Week.end.last2.log)=-0.0489"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.end.last2.log as highly correlated with
## Week.bgn.last2.log
```

```
## [1] "cor(Week.bgn.month.fctr.nonNA, Week.end.month.fctr)=1.0000"
## [1] "cor(ILI.log, Week.bgn.month.fctr.nonNA)=-0.1979"
## [1] "cor(ILI.log, Week.end.month.fctr)=-0.1979"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.end.month.fctr as highly correlated with
## Week.bgn.month.fctr.nonNA
```

```
## [1] "cor(Week.bgn.year.fctr.nonNA, Week.end.year.fctr)=1.0000"
## [1] "cor(ILI.log, Week.bgn.year.fctr.nonNA)=0.1987"
## [1] "cor(ILI.log, Week.end.year.fctr)=0.1987"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.end.year.fctr as highly correlated with
## Week.bgn.year.fctr.nonNA
```

```
## [1] "cor(ILI.2.lag.log.nonNA, Queries)=0.7424"
## [1] "cor(ILI.log, ILI.2.lag.log.nonNA)=0.9204"
## [1] "cor(ILI.log, Queries)=0.8420"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Queries as highly correlated with ILI.
## 2.lag.log.nonNA
```

```
## [1] "cor(Week.bgn.last100.log, Week.bgn.year.fctr.nonNA)=0.7399"
## [1] "cor(ILI.log, Week.bgn.last100.log)=0.2001"
## [1] "cor(ILI.log, Week.bgn.year.fctr.nonNA)=0.1987"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.bgn.year.fctr.nonNA as highly correlated
## with Week.bgn.last100.log
```

```
## [1] "cor(Week.bgn.last1.log, Week.bgn.last2.log)=0.7063"
## [1] "cor(ILI.log, Week.bgn.last1.log)=-0.0474"
## [1] "cor(ILI.log, Week.bgn.last2.log)=-0.0489"
```

```
## Warning in myfind_cor_features(feats_df = glb_feats_df, obs_df =
## glb_trnobs_df, : Identified Week.bgn.last1.log as highly correlated with
## Week.bgn.last2.log
```

```
##                           id         cor.y exclude.as.feat    cor.y.abs
## 2                        ILI  0.9451682372               1 0.9451682372
## 4              ILI.2.lag.log  0.9214355113               1 0.9214355113
## 5        ILI.2.lag.log.nonNA  0.9204139777               0 0.9204139777
## 3                  ILI.2.lag  0.8591322165               1 0.8591322165
## 6                    Queries  0.8420332864               0 0.8420332864
## 12      Week.bgn.last100.log  0.2000557734               0 0.2000557734
## 28      Week.end.last100.log  0.2000557734               0 0.2000557734
## 21        Week.bgn.year.fctr  0.1987250820               1 0.1987250820
## 22  Week.bgn.year.fctr.nonNA  0.1987250820               0 0.1987250820
## 36        Week.end.year.fctr  0.1987250820               0 0.1987250820
## 17            Week.bgn.POSIX  0.1730840769               1 0.1730840769
## 23              Week.bgn.zoo  0.1730840769               1 0.1730840769
## 32            Week.end.POSIX  0.1730840769               1 0.1730840769
## 37              Week.end.zoo  0.1730840769               1 0.1730840769
## 1                     .rnorm  0.0400349990               0 0.0400349990
## 7         Week.bgn.date.fctr  0.0088237363               1 0.0088237363
## 8   Week.bgn.date.fctr.nonNA  0.0088237363               0 0.0088237363
## 24        Week.end.date.fctr  0.0088237363               0 0.0088237363
## 11       Week.bgn.last10.log  0.0006098826               0 0.0006098826
## 27       Week.end.last10.log  0.0006098826               0 0.0006098826
## 10        Week.bgn.last1.log -0.0473819651               0 0.0473819651
## 26        Week.end.last1.log -0.0473819651               0 0.0473819651
## 13        Week.bgn.last2.log -0.0488875907               0 0.0488875907
## 29        Week.end.last2.log -0.0488875907               0 0.0488875907
## 15       Week.bgn.month.fctr -0.1979117680               1 0.1979117680
## 16 Week.bgn.month.fctr.nonNA -0.1979117680               0 0.1979117680
## 31       Week.end.month.fctr -0.1979117680               0 0.1979117680
## 9         Week.bgn.hour.fctr            NA               0           NA
## 14      Week.bgn.minute.fctr            NA               0           NA
## 18      Week.bgn.second.fctr            NA               0           NA
## 19       Week.bgn.wkday.fctr            NA               0           NA
## 20            Week.bgn.wkend            NA               0           NA
## 25        Week.end.hour.fctr            NA               0           NA
## 30      Week.end.minute.fctr            NA               0           NA
## 33      Week.end.second.fctr            NA               0           NA
## 34       Week.end.wkday.fctr            NA               0           NA
## 35            Week.end.wkend            NA               0           NA
##                   cor.high.X freqRatio percentUnique zeroVar   nzv
## 2                       <NA>  1.000000   100.0000000   FALSE FALSE
## 4                       <NA>  1.000000    99.5203837   FALSE FALSE
## 5                       <NA>  2.000000    99.7601918   FALSE FALSE
## 3                       <NA>  1.000000    99.5203837   FALSE FALSE
## 6        ILI.2.lag.log.nonNA  1.250000    62.5899281   FALSE FALSE
## 12                      <NA>  2.690000     0.9592326   FALSE FALSE
## 28      Week.bgn.last100.log  2.690000     0.9592326   FALSE FALSE
## 21                      <NA>  1.019231     1.9184652   FALSE FALSE
## 22      Week.bgn.last100.log  1.019231     1.9184652   FALSE FALSE
## 36  Week.bgn.year.fctr.nonNA  1.019231     1.9184652   FALSE FALSE
## 17                      <NA>  1.000000   100.0000000   FALSE FALSE
## 23                      <NA>  1.000000   100.0000000   FALSE FALSE
## 32                      <NA>  1.000000   100.0000000   FALSE FALSE
## 37                      <NA>  1.000000   100.0000000   FALSE FALSE
## 1                       <NA>  1.000000   100.0000000   FALSE FALSE
## 7                       <NA>  1.156627     1.1990408   FALSE FALSE
## 8                       <NA>  1.156627     1.1990408   FALSE FALSE
## 24                      <NA>  1.156627     1.1990408   FALSE FALSE
## 11                      <NA>  3.125000     0.9592326   FALSE FALSE
## 27                      <NA>  3.125000     0.9592326   FALSE FALSE
## 10        Week.bgn.last2.log 50.000000     0.9592326   FALSE  TRUE
## 26        Week.bgn.last1.log 50.000000     0.9592326   FALSE  TRUE
## 13                      <NA> 23.937500     0.9592326   FALSE  TRUE
## 29        Week.bgn.last2.log 23.937500     0.9592326   FALSE  TRUE
## 15                      <NA>  1.000000     2.8776978   FALSE FALSE
## 16                      <NA>  1.000000     2.8776978   FALSE FALSE
## 31 Week.bgn.month.fctr.nonNA  1.000000     2.8776978   FALSE FALSE
## 9                       <NA>  0.000000     0.2398082    TRUE  TRUE
## 14                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 18                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 19                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 20                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 25                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 30                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 33                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 34                      <NA>  0.000000     0.2398082    TRUE  TRUE
## 35                      <NA>  0.000000     0.2398082    TRUE  TRUE
##    myNearZV is.cor.y.abs.low
## 2     FALSE            FALSE
## 4     FALSE            FALSE
## 5     FALSE            FALSE
## 3     FALSE            FALSE
## 6     FALSE            FALSE
## 12    FALSE            FALSE
## 28    FALSE            FALSE
## 21    FALSE            FALSE
## 22    FALSE            FALSE
## 36    FALSE            FALSE
## 17    FALSE            FALSE
## 23    FALSE            FALSE
## 32    FALSE            FALSE
## 37    FALSE            FALSE
## 1     FALSE            FALSE
## 7     FALSE             TRUE
## 8     FALSE             TRUE
## 24    FALSE             TRUE
## 11    FALSE             TRUE
## 27    FALSE             TRUE
## 10    FALSE            FALSE
## 26    FALSE            FALSE
## 13    FALSE            FALSE
## 29    FALSE            FALSE
## 15    FALSE            FALSE
## 16    FALSE            FALSE
## 31    FALSE            FALSE
## 9      TRUE               NA
## 14     TRUE               NA
## 18     TRUE               NA
## 19     TRUE               NA
## 20     TRUE               NA
## 25     TRUE               NA
## 30     TRUE               NA
## 33     TRUE               NA
## 34     TRUE               NA
## 35     TRUE               NA
```

```r
#subset(glb_feats_df, id %in% c("A.nuppr.log", "S.nuppr.log"))
print(myplot_scatter(glb_feats_df, "percentUnique", "freqRatio", 
                     colorcol_name="myNearZV", jitter=TRUE) + 
          geom_point(aes(shape=nzv)) + xlim(-5, 25))
```

```
## Warning in myplot_scatter(glb_feats_df, "percentUnique", "freqRatio",
## colorcol_name = "myNearZV", : converting myNearZV to class:factor
```

```
## Warning in loop_apply(n, do.ply): Removed 10 rows containing missing values
## (geom_point).
```

```
## Warning in loop_apply(n, do.ply): Removed 10 rows containing missing values
## (geom_point).
```

```
## Warning in loop_apply(n, do.ply): Removed 10 rows containing missing values
## (geom_point).
```

![](Google_Flu_template2_files/figure-html/select.features-1.png) 

```r
print(subset(glb_feats_df, myNearZV))
```

```
##                      id cor.y exclude.as.feat cor.y.abs cor.high.X
## 9    Week.bgn.hour.fctr    NA               0        NA       <NA>
## 14 Week.bgn.minute.fctr    NA               0        NA       <NA>
## 18 Week.bgn.second.fctr    NA               0        NA       <NA>
## 19  Week.bgn.wkday.fctr    NA               0        NA       <NA>
## 20       Week.bgn.wkend    NA               0        NA       <NA>
## 25   Week.end.hour.fctr    NA               0        NA       <NA>
## 30 Week.end.minute.fctr    NA               0        NA       <NA>
## 33 Week.end.second.fctr    NA               0        NA       <NA>
## 34  Week.end.wkday.fctr    NA               0        NA       <NA>
## 35       Week.end.wkend    NA               0        NA       <NA>
##    freqRatio percentUnique zeroVar  nzv myNearZV is.cor.y.abs.low
## 9          0     0.2398082    TRUE TRUE     TRUE               NA
## 14         0     0.2398082    TRUE TRUE     TRUE               NA
## 18         0     0.2398082    TRUE TRUE     TRUE               NA
## 19         0     0.2398082    TRUE TRUE     TRUE               NA
## 20         0     0.2398082    TRUE TRUE     TRUE               NA
## 25         0     0.2398082    TRUE TRUE     TRUE               NA
## 30         0     0.2398082    TRUE TRUE     TRUE               NA
## 33         0     0.2398082    TRUE TRUE     TRUE               NA
## 34         0     0.2398082    TRUE TRUE     TRUE               NA
## 35         0     0.2398082    TRUE TRUE     TRUE               NA
```

```r
glb_allobs_df <- glb_allobs_df[, setdiff(names(glb_allobs_df), 
                                         subset(glb_feats_df, myNearZV)$id)]

if (!is.null(glb_interaction_only_features))
    glb_feats_df[glb_feats_df$id %in% glb_interaction_only_features, "interaction.feat"] <-
        names(glb_interaction_only_features) else
    glb_feats_df$interaction.feat <- NA        

mycheck_problem_data(glb_allobs_df, terminate = TRUE)
```

```
## [1] "numeric data missing in : "
##     ILI.2.lag ILI.2.lag.log 
##             2             2 
## [1] "numeric data w/ 0s in : "
##   Week.bgn.last1.log   Week.bgn.last2.log  Week.bgn.last10.log 
##                    1                    2                   10 
## Week.bgn.last100.log   Week.end.last1.log   Week.end.last2.log 
##                  100                    1                    2 
##  Week.end.last10.log Week.end.last100.log 
##                   10                  100 
## [1] "numeric data w/ Infs in : "
## named integer(0)
## [1] "numeric data w/ NaNs in : "
## named integer(0)
## [1] "string data missing in : "
##     Week Week.bgn Week.end 
##        0        0        0
```

```r
# glb_allobs_df %>% filter(is.na(Married.fctr)) %>% tbl_df()
# glb_allobs_df %>% count(Married.fctr)
# levels(glb_allobs_df$Married.fctr)

glb_chunks_df <- myadd_chunk(glb_chunks_df, "partition.data.training", major.inc=TRUE)
```

```
##                     label step_major step_minor     bgn     end elapsed
## 8         select.features          5          0 128.258 128.999   0.742
## 9 partition.data.training          6          0 129.000      NA      NA
```

## Step `6.0: partition data training`

```r
if (all(is.na(glb_newobs_df[, glb_rsp_var]))) {
    require(caTools)
    
    set.seed(glb_split_sample.seed)
    split <- sample.split(glb_trnobs_df[, glb_rsp_var_raw], 
        SplitRatio=1 - (nrow(glb_newobs_df) * 1.1 / nrow(glb_trnobs_df)))
    glb_fitobs_df <- glb_trnobs_df[split, ] 
    glb_OOBobs_df <- glb_trnobs_df[!split ,]    
} else {
    print(sprintf("Newdata contains non-NA data for %s; setting OOB to Newdata", 
                  glb_rsp_var))
    glb_fitobs_df <- glb_trnobs_df; glb_OOBobs_df <- glb_newobs_df
}
```

```
## [1] "Newdata contains non-NA data for ILI.log; setting OOB to Newdata"
```

```r
if (!is.null(glb_max_fitobs) && (nrow(glb_fitobs_df) > glb_max_fitobs)) {
    warning("glb_fitobs_df restricted to glb_max_fitobs: ", 
            format(glb_max_fitobs, big.mark=","))
    org_fitobs_df <- glb_fitobs_df
    glb_fitobs_df <- 
        org_fitobs_df[split <- sample.split(org_fitobs_df[, glb_rsp_var_raw], 
                                            SplitRatio=glb_max_fitobs), ]
    org_fitobs_df <- NULL
}

glb_allobs_df$.lcn <- ""
glb_allobs_df[glb_allobs_df[, glb_id_var] %in% 
              glb_fitobs_df[, glb_id_var], ".lcn"] <- "Fit"
glb_allobs_df[glb_allobs_df[, glb_id_var] %in% 
              glb_OOBobs_df[, glb_id_var], ".lcn"] <- "OOB"

dsp_class_dstrb <- function(obs_df, location_var, partition_var) {
    xtab_df <- mycreate_xtab_df(obs_df, c(location_var, partition_var))
    rownames(xtab_df) <- xtab_df[, location_var]
    xtab_df <- xtab_df[, -grepl(location_var, names(xtab_df))]
    print(xtab_df)
    print(xtab_df / rowSums(xtab_df, na.rm=TRUE))    
}    

# Ensure proper splits by glb_rsp_var_raw & user-specified feature for OOB vs. new
if (!is.null(glb_category_vars)) {
    if (glb_is_classification)
        dsp_class_dstrb(glb_allobs_df, ".lcn", glb_rsp_var_raw)
    newobs_ctgry_df <- mycreate_sqlxtab_df(subset(glb_allobs_df, .src == "Test"), 
                                           glb_category_vars)
    OOBobs_ctgry_df <- mycreate_sqlxtab_df(subset(glb_allobs_df, .lcn == "OOB"), 
                                           glb_category_vars)
    glb_ctgry_df <- merge(newobs_ctgry_df, OOBobs_ctgry_df, by=glb_category_vars
                          , all=TRUE, suffixes=c(".Tst", ".OOB"))
    glb_ctgry_df$.freqRatio.Tst <- glb_ctgry_df$.n.Tst / sum(glb_ctgry_df$.n.Tst, na.rm=TRUE)
    glb_ctgry_df$.freqRatio.OOB <- glb_ctgry_df$.n.OOB / sum(glb_ctgry_df$.n.OOB, na.rm=TRUE)
    print(orderBy(~-.freqRatio.Tst-.freqRatio.OOB, glb_ctgry_df))
}

# Run this line by line
print("glb_feats_df:");   print(dim(glb_feats_df))
```

```
## [1] "glb_feats_df:"
```

```
## [1] 37 12
```

```r
sav_feats_df <- glb_feats_df
glb_feats_df <- sav_feats_df

glb_feats_df[, "rsp_var_raw"] <- FALSE
glb_feats_df[glb_feats_df$id == glb_rsp_var_raw, "rsp_var_raw"] <- TRUE 
glb_feats_df$exclude.as.feat <- (glb_feats_df$exclude.as.feat == 1)
if (!is.null(glb_id_var) && glb_id_var != ".rownames")
    glb_feats_df[glb_feats_df$id %in% glb_id_var, "id_var"] <- TRUE 
add_feats_df <- data.frame(id=glb_rsp_var, exclude.as.feat=TRUE, rsp_var=TRUE)
row.names(add_feats_df) <- add_feats_df$id; print(add_feats_df)
```

```
##              id exclude.as.feat rsp_var
## ILI.log ILI.log            TRUE    TRUE
```

```r
glb_feats_df <- myrbind_df(glb_feats_df, add_feats_df)
if (glb_id_var != ".rownames")
    print(subset(glb_feats_df, rsp_var_raw | rsp_var | id_var)) else
    print(subset(glb_feats_df, rsp_var_raw | rsp_var))    
```

```
##              id     cor.y exclude.as.feat cor.y.abs cor.high.X freqRatio
## 2           ILI 0.9451682            TRUE 0.9451682       <NA>         1
## ILI.log ILI.log        NA            TRUE        NA       <NA>        NA
##         percentUnique zeroVar   nzv myNearZV is.cor.y.abs.low
## 2                 100   FALSE FALSE    FALSE            FALSE
## ILI.log            NA      NA    NA       NA               NA
##         interaction.feat rsp_var_raw id_var rsp_var
## 2                     NA        TRUE     NA      NA
## ILI.log               NA          NA     NA    TRUE
```

```r
print("glb_feats_df vs. glb_allobs_df: "); 
```

```
## [1] "glb_feats_df vs. glb_allobs_df: "
```

```r
print(setdiff(glb_feats_df$id, names(glb_allobs_df)))
```

```
##  [1] "Week.bgn.hour.fctr"   "Week.bgn.minute.fctr" "Week.bgn.second.fctr"
##  [4] "Week.bgn.wkday.fctr"  "Week.bgn.wkend"       "Week.end.hour.fctr"  
##  [7] "Week.end.minute.fctr" "Week.end.second.fctr" "Week.end.wkday.fctr" 
## [10] "Week.end.wkend"
```

```r
print("glb_allobs_df vs. glb_feats_df: "); 
```

```
## [1] "glb_allobs_df vs. glb_feats_df: "
```

```r
# Ensure these are only chr vars
print(setdiff(setdiff(names(glb_allobs_df), glb_feats_df$id), 
                myfind_chr_cols_df(glb_allobs_df)))
```

```
## character(0)
```

```r
#print(setdiff(setdiff(names(glb_allobs_df), glb_exclude_vars_as_features), 
#                glb_feats_df$id))

print("glb_allobs_df: "); print(dim(glb_allobs_df))
```

```
## [1] "glb_allobs_df: "
```

```
## [1] 469  33
```

```r
print("glb_trnobs_df: "); print(dim(glb_trnobs_df))
```

```
## [1] "glb_trnobs_df: "
```

```
## [1] 417  42
```

```r
print("glb_fitobs_df: "); print(dim(glb_fitobs_df))
```

```
## [1] "glb_fitobs_df: "
```

```
## [1] 417  42
```

```r
print("glb_OOBobs_df: "); print(dim(glb_OOBobs_df))
```

```
## [1] "glb_OOBobs_df: "
```

```
## [1] 52 42
```

```r
print("glb_newobs_df: "); print(dim(glb_newobs_df))
```

```
## [1] "glb_newobs_df: "
```

```
## [1] 52 42
```

```r
# # Does not handle NULL or length(glb_id_var) > 1
# glb_allobs_df$.src.trn <- 0
# glb_allobs_df[glb_allobs_df[, glb_id_var] %in% glb_trnobs_df[, glb_id_var], 
#                 ".src.trn"] <- 1 
# glb_allobs_df$.src.fit <- 0
# glb_allobs_df[glb_allobs_df[, glb_id_var] %in% glb_fitobs_df[, glb_id_var], 
#                 ".src.fit"] <- 1 
# glb_allobs_df$.src.OOB <- 0
# glb_allobs_df[glb_allobs_df[, glb_id_var] %in% glb_OOBobs_df[, glb_id_var], 
#                 ".src.OOB"] <- 1 
# glb_allobs_df$.src.new <- 0
# glb_allobs_df[glb_allobs_df[, glb_id_var] %in% glb_newobs_df[, glb_id_var], 
#                 ".src.new"] <- 1 
# #print(unique(glb_allobs_df[, ".src.trn"]))
# write_cols <- c(glb_feats_df$id, 
#                 ".src.trn", ".src.fit", ".src.OOB", ".src.new")
# glb_allobs_df <- glb_allobs_df[, write_cols]
# 
# tmp_feats_df <- glb_feats_df
# tmp_entity_df <- glb_allobs_df

if (glb_save_envir)
    save(glb_feats_df, 
         glb_allobs_df, #glb_trnobs_df, glb_fitobs_df, glb_OOBobs_df, glb_newobs_df,
         file=paste0(glb_out_pfx, "blddfs_dsk.RData"))
# load(paste0(glb_out_pfx, "blddfs_dsk.RData"))

# if (!all.equal(tmp_feats_df, glb_feats_df))
#     stop("glb_feats_df r/w not working")
# if (!all.equal(tmp_entity_df, glb_allobs_df))
#     stop("glb_allobs_df r/w not working")

rm(split)
```

```
## Warning in rm(split): object 'split' not found
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=TRUE)
```

```
##                      label step_major step_minor     bgn     end elapsed
## 9  partition.data.training          6          0 129.000 129.347   0.347
## 10              fit.models          7          0 129.347      NA      NA
```

## Step `7.0: fit models`

```r
# load(paste0(glb_out_pfx, "dsk.RData"))
# keep_cols <- setdiff(names(glb_allobs_df), 
#                      grep("^.src", names(glb_allobs_df), value=TRUE))
# glb_trnobs_df <- glb_allobs_df[glb_allobs_df$.src.trn == 1, keep_cols]
# glb_fitobs_df <- glb_allobs_df[glb_allobs_df$.src.fit == 1, keep_cols]
# glb_OOBobs_df <- glb_allobs_df[glb_allobs_df$.src.OOB == 1, keep_cols]
# glb_newobs_df <- glb_allobs_df[glb_allobs_df$.src.new == 1, keep_cols]
# 
# glb_models_lst <- list(); glb_models_df <- data.frame()
# 
if (glb_is_classification && glb_is_binomial && 
        (length(unique(glb_fitobs_df[, glb_rsp_var])) < 2))
    stop("glb_fitobs_df$", glb_rsp_var, ": contains less than 2 unique values: ",
         paste0(unique(glb_fitobs_df[, glb_rsp_var]), collapse=", "))

max_cor_y_x_vars <- orderBy(~ -cor.y.abs, 
        subset(glb_feats_df, (exclude.as.feat == 0) & !is.cor.y.abs.low & 
                                is.na(cor.high.X)))[1:2, "id"]
# while(length(max_cor_y_x_vars) < 2) {
#     max_cor_y_x_vars <- c(max_cor_y_x_vars, orderBy(~ -cor.y.abs, 
#             subset(glb_feats_df, (exclude.as.feat == 0) & !is.cor.y.abs.low))[3, "id"])    
# }
if (!is.null(glb_Baseline_mdl_var)) {
    if ((max_cor_y_x_vars[1] != glb_Baseline_mdl_var) & 
        (glb_feats_df[max_cor_y_x_vars[1], "cor.y.abs"] > 
         glb_feats_df[glb_Baseline_mdl_var, "cor.y.abs"]))
        stop(max_cor_y_x_vars[1], " has a lower correlation with ", glb_rsp_var, 
             " than the Baseline var: ", glb_Baseline_mdl_var)
}

glb_model_type <- ifelse(glb_is_regression, "regression", "classification")
    
# Baseline
if (!is.null(glb_Baseline_mdl_var)) 
    ret_lst <- myfit_mdl_fn(model_id="Baseline", model_method="mybaseln_classfr",
                            indep_vars_vctr=glb_Baseline_mdl_var,
                            rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                            fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df)

# Most Frequent Outcome "MFO" model: mean(y) for regression
#   Not using caret's nullModel since model stats not avl
#   Cannot use rpart for multinomial classification since it predicts non-MFO
ret_lst <- myfit_mdl(model_id="MFO", 
                     model_method=ifelse(glb_is_regression, "lm", "myMFO_classfr"), 
                     model_type=glb_model_type,
                        indep_vars_vctr=".rnorm",
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df)
```

```
## [1] "fitting model: MFO.lm"
## [1] "    indep_vars: .rnorm"
## Fitting parameter = none on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_0-1.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-2.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-3.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-4.png) 

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.9780 -0.4535 -0.1217  0.3793  1.6924 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.34737    0.02717  12.784   <2e-16 ***
## .rnorm       0.02347    0.02875   0.816    0.415    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.5548 on 415 degrees of freedom
## Multiple R-squared:  0.001603,	Adjusted R-squared:  -0.000803 
## F-statistic: 0.6662 on 1 and 415 DF,  p-value: 0.4148
## 
## [1] "    calling mypredict_mdl for fit:"
## [1] "    calling mypredict_mdl for OOB:"
##   model_id model_method  feats max.nTuningRuns min.elapsedtime.everything
## 1   MFO.lm           lm .rnorm               0                      0.566
##   min.elapsedtime.final max.R.sq.fit min.RMSE.fit max.R.sq.OOB
## 1                 0.002  0.001602801    0.5534848  -0.01575461
##   min.RMSE.OOB max.Adj.R.sq.fit
## 1    0.4080698    -0.0008029752
```

```r
if (glb_is_classification)
    # "random" model - only for classification; 
    #   none needed for regression since it is same as MFO
    ret_lst <- myfit_mdl(model_id="Random", model_method="myrandom_classfr",
                            model_type=glb_model_type,                         
                            indep_vars_vctr=".rnorm",
                            rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                            fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df)

# Any models that have tuning parameters has "better" results with cross-validation
#   (except rf) & "different" results for different outcome metrics

# Max.cor.Y
#   Check impact of cv
#       rpart is not a good candidate since caret does not optimize cp (only tuning parameter of rpart) well
ret_lst <- myfit_mdl(model_id="Max.cor.Y.cv.0", 
                        model_method="rpart",
                     model_type=glb_model_type,
                        indep_vars_vctr=max_cor_y_x_vars,
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df)
```

```
## [1] "fitting model: Max.cor.Y.cv.0.rpart"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log"
```

```
## Loading required package: rpart
```

```
## Fitting cp = 0.635 on full training set
```

```
## Loading required package: rpart.plot
```

![](Google_Flu_template2_files/figure-html/fit.models_0-5.png) 

```
## Call:
## rpart(formula = .outcome ~ ., control = list(minsplit = 20, minbucket = 7, 
##     cp = 0, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
##     surrogatestyle = 0, maxdepth = 30, xval = 0))
##   n= 417 
## 
##          CP nsplit rel error
## 1 0.6353105      0         1
## 
## Node number 1: 417 observations
##   mean=0.3476702, MSE=0.3068372 
## 
## n= 417 
## 
## node), split, n, deviance, yval
##       * denotes terminal node
## 
## 1) root 417 127.9511 0.3476702 *
## [1] "    calling mypredict_mdl for fit:"
## [1] "    calling mypredict_mdl for OOB:"
##               model_id model_method
## 1 Max.cor.Y.cv.0.rpart        rpart
##                                       feats max.nTuningRuns
## 1 ILI.2.lag.log.nonNA, Week.bgn.last100.log               0
##   min.elapsedtime.everything min.elapsedtime.final max.R.sq.fit
## 1                      0.547                 0.016            0
##   min.RMSE.fit max.R.sq.OOB min.RMSE.OOB
## 1    0.5539289            0    0.4048928
```

```r
ret_lst <- myfit_mdl(model_id="Max.cor.Y.cv.0.cp.0", 
                        model_method="rpart",
                     model_type=glb_model_type,
                        indep_vars_vctr=max_cor_y_x_vars,
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                        n_cv_folds=0, 
            tune_models_df=data.frame(parameter="cp", min=0.0, max=0.0, by=0.1))
```

```
## [1] "fitting model: Max.cor.Y.cv.0.cp.0.rpart"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log"
## Fitting cp = 0 on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_0-6.png) 

```
## Call:
## rpart(formula = .outcome ~ ., control = list(minsplit = 20, minbucket = 7, 
##     cp = 0, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
##     surrogatestyle = 0, maxdepth = 30, xval = 0))
##   n= 417 
## 
##              CP nsplit rel error
## 1  6.353105e-01      0 1.0000000
## 2  9.235980e-02      1 0.3646895
## 3  6.881570e-02      2 0.2723297
## 4  2.714155e-02      3 0.2035140
## 5  9.739677e-03      4 0.1763725
## 6  9.553531e-03      5 0.1666328
## 7  8.328212e-03      6 0.1570793
## 8  4.575803e-03      7 0.1487511
## 9  2.026905e-03      8 0.1441753
## 10 1.920792e-03      9 0.1421484
## 11 1.557624e-03     10 0.1402276
## 12 1.495507e-03     12 0.1371123
## 13 1.379224e-03     14 0.1341213
## 14 1.366310e-03     15 0.1327421
## 15 9.483927e-04     16 0.1313758
## 16 8.321743e-04     17 0.1304274
## 17 8.288942e-04     19 0.1287630
## 18 4.238163e-04     20 0.1279341
## 19 3.274522e-04     21 0.1275103
## 20 3.176373e-04     22 0.1271829
## 21 2.254809e-04     24 0.1265476
## 22 2.223878e-04     25 0.1263221
## 23 1.518224e-04     26 0.1260997
## 24 1.222121e-04     27 0.1259479
## 25 8.532486e-05     28 0.1258257
## 26 0.000000e+00     29 0.1257404
## 
## Variable importance
##  ILI.2.lag.log.nonNA Week.bgn.last100.log 
##                   92                    8 
## 
## Node number 1: 417 observations,    complexity param=0.6353105
##   mean=0.3476702, MSE=0.3068372 
##   left son=2 (238 obs) right son=3 (179 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.3642121    to the left,  improve=0.6353105, (0 missing)
##       Week.bgn.last100.log < 17.91785     to the left,  improve=0.0514081, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91785     to the left,  agree=0.619, adj=0.112, (0 split)
## 
## Node number 2: 238 observations,    complexity param=0.0688157
##   mean=-0.03522948, MSE=0.07240075 
##   left son=4 (135 obs) right son=5 (103 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < -0.005288628 to the left,  improve=0.51098900, (0 missing)
##       Week.bgn.last100.log < 17.91779     to the right, improve=0.01577173, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91779     to the right, agree=0.601, adj=0.078, (0 split)
## 
## Node number 3: 179 observations,    complexity param=0.0923598
##   mean=0.8567771, MSE=0.1644193 
##   left son=6 (141 obs) right son=7 (38 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 1.160636     to the left,  improve=0.40153300, (0 missing)
##       Week.bgn.last100.log < 17.91779     to the right, improve=0.02044471, (0 missing)
## 
## Node number 4: 135 observations,    complexity param=0.009739677
##   mean=-0.2032371, MSE=0.02669874 
##   left son=8 (46 obs) right son=9 (89 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < -0.2557664   to the left,  improve=0.34575160, (0 missing)
##       Week.bgn.last100.log < 8.958882     to the left,  improve=0.01554574, (0 missing)
## 
## Node number 5: 103 observations,    complexity param=0.008328212
##   mean=0.1849747, MSE=0.04681557 
##   left son=10 (66 obs) right son=11 (37 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.1886472    to the left,  improve=0.220987800, (0 missing)
##       Week.bgn.last100.log < 8.958882     to the left,  improve=0.003546195, (0 missing)
## 
## Node number 6: 141 observations,    complexity param=0.02714155
##   mean=0.7233883, MSE=0.098441 
##   left son=12 (84 obs) right son=13 (57 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.7871896    to the left,  improve=0.25019790, (0 missing)
##       Week.bgn.last100.log < 17.91785     to the left,  improve=0.00108719, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91785     to the left,  agree=0.603, adj=0.018, (0 split)
## 
## Node number 7: 38 observations,    complexity param=0.009553531
##   mean=1.35172, MSE=0.09824558 
##   left son=14 (24 obs) right son=15 (14 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 1.485112     to the left,  improve=0.3274247, (0 missing)
##       Week.bgn.last100.log < 17.91779     to the right, improve=0.1548345, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91779     to the right, agree=0.711, adj=0.214, (0 split)
## 
## Node number 8: 46 observations,    complexity param=0.001379224
##   mean=-0.3368793, MSE=0.01754249 
##   left son=16 (28 obs) right son=17 (18 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < -0.3364944   to the left,  improve=0.2186906000, (0 missing)
##       Week.bgn.last100.log < 8.958912     to the left,  improve=0.0001459021, (0 missing)
## 
## Node number 9: 89 observations,    complexity param=0.00136631
##   mean=-0.1341636, MSE=0.0174289 
##   left son=18 (53 obs) right son=19 (36 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < -0.1187063   to the left,  improve=0.11270240, (0 missing)
##       Week.bgn.last100.log < 8.958882     to the left,  improve=0.02457038, (0 missing)
## 
## Node number 10: 66 observations,    complexity param=0.002026905
##   mean=0.108818, MSE=0.02579353 
##   left son=20 (35 obs) right son=21 (31 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.1054685    to the left,  improve=0.152343000, (0 missing)
##       Week.bgn.last100.log < 8.958882     to the left,  improve=0.008691956, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91779     to the left,  agree=0.652, adj=0.258, (0 split)
## 
## Node number 11: 37 observations,    complexity param=0.0009483927
##   mean=0.3208218, MSE=0.05551424 
##   left son=22 (18 obs) right son=23 (19 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.2672436    to the right, improve=0.05907805, (0 missing)
##       Week.bgn.last100.log < 8.958882     to the left,  improve=0.01362272, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 8.958882     to the left,  agree=0.595, adj=0.167, (0 split)
## 
## Node number 12: 84 observations,    complexity param=0.004575803
##   mean=0.5941095, MSE=0.07154071 
##   left son=24 (26 obs) right son=25 (58 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.4780125    to the left,  improve=0.0974268900, (0 missing)
##       Week.bgn.last100.log < 8.958912     to the left,  improve=0.0007104286, (0 missing)
## 
## Node number 13: 57 observations,    complexity param=0.0008321743
##   mean=0.9139046, MSE=0.07715736 
##   left son=26 (45 obs) right son=27 (12 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.8769581    to the right, improve=0.01875118, (0 missing)
##       Week.bgn.last100.log < 17.91785     to the right, improve=0.00163561, (0 missing)
## 
## Node number 14: 24 observations,    complexity param=0.001920792
##   mean=1.214735, MSE=0.06238431 
##   left son=28 (14 obs) right son=29 (10 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 1.367642     to the left,  improve=0.1641489, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 8.958912     to the right, agree=0.625, adj=0.1, (0 split)
## 
## Node number 15: 14 observations
##   mean=1.58655, MSE=0.07240882 
## 
## Node number 16: 28 observations,    complexity param=0.0002223878
##   mean=-0.3865406, MSE=0.01641498 
##   left son=32 (15 obs) right son=33 (13 obs)
##   Primary splits:
##       Week.bgn.last100.log < 8.958912     to the right, improve=0.06190942, (0 missing)
##       ILI.2.lag.log.nonNA  < -0.4027548   to the left,  improve=0.04806694, (0 missing)
##   Surrogate splits:
##       ILI.2.lag.log.nonNA < -0.3889181   to the left,  agree=0.643, adj=0.231, (0 split)
## 
## Node number 17: 18 observations
##   mean=-0.2596284, MSE=0.009492319 
## 
## Node number 18: 53 observations,    complexity param=0.0002254809
##   mean=-0.1706907, MSE=0.01735997 
##   left son=36 (19 obs) right son=37 (34 obs)
##   Primary splits:
##       Week.bgn.last100.log < 8.958882     to the left,  improve=0.03135659, (0 missing)
##       ILI.2.lag.log.nonNA  < -0.2037076   to the right, improve=0.02610092, (0 missing)
##   Surrogate splits:
##       ILI.2.lag.log.nonNA < -0.2068175   to the left,  agree=0.698, adj=0.158, (0 split)
## 
## Node number 19: 36 observations,    complexity param=8.532486e-05
##   mean=-0.08038768, MSE=0.01267426 
##   left son=38 (18 obs) right son=39 (18 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < -0.06643646  to the left,  improve=0.023927350, (0 missing)
##       Week.bgn.last100.log < 17.91779     to the right, improve=0.006995762, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91779     to the left,  agree=0.528, adj=0.056, (0 split)
## 
## Node number 20: 35 observations,    complexity param=0.0004238163
##   mean=0.04982323, MSE=0.01327698 
##   left son=40 (18 obs) right son=41 (17 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.05575614   to the left,  improve=0.11669560, (0 missing)
##       Week.bgn.last100.log < 17.91779     to the right, improve=0.07230438, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 8.958882     to the left,  agree=0.571, adj=0.118, (0 split)
## 
## Node number 21: 31 observations,    complexity param=0.0003176373
##   mean=0.175425, MSE=0.03155917 
##   left son=42 (11 obs) right son=43 (20 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.1561068    to the right, improve=0.03380649, (0 missing)
##       Week.bgn.last100.log < 17.91779     to the right, improve=0.00791479, (0 missing)
## 
## Node number 22: 18 observations
##   mean=0.2619841, MSE=0.04141873 
## 
## Node number 23: 19 observations
##   mean=0.3765628, MSE=0.06248115 
## 
## Node number 24: 26 observations,    complexity param=0.0008288942
##   mean=0.4694161, MSE=0.04158716 
##   left son=48 (8 obs) right son=49 (18 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.4304485    to the right, improve=0.09808679, (0 missing)
## 
## Node number 25: 58 observations,    complexity param=0.001557624
##   mean=0.6500065, MSE=0.0748737 
##   left son=50 (35 obs) right son=51 (23 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.6743843    to the left,  improve=0.03720445, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91785     to the left,  agree=0.672, adj=0.174, (0 split)
## 
## Node number 26: 45 observations,    complexity param=0.0008321743
##   mean=0.8942625, MSE=0.07023539 
##   left son=52 (27 obs) right son=53 (18 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 1.026878     to the left,  improve=0.0412859700, (0 missing)
##       Week.bgn.last100.log < 17.91785     to the right, improve=0.0002826772, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91785     to the left,  agree=0.622, adj=0.056, (0 split)
## 
## Node number 27: 12 observations
##   mean=0.9875624, MSE=0.09624252 
## 
## Node number 28: 14 observations
##   mean=1.12921, MSE=0.04984909 
## 
## Node number 29: 10 observations
##   mean=1.33447, MSE=0.05535686 
## 
## Node number 32: 15 observations
##   mean=-0.4162179, MSE=0.01511771 
## 
## Node number 33: 13 observations
##   mean=-0.3522975, MSE=0.01572299 
## 
## Node number 36: 19 observations
##   mean=-0.2019013, MSE=0.01969646 
## 
## Node number 37: 34 observations,    complexity param=0.0001518224
##   mean=-0.1532495, MSE=0.01520574 
##   left son=74 (9 obs) right son=75 (25 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < -0.1486436   to the right, improve=0.03757451, (0 missing)
## 
## Node number 38: 18 observations
##   mean=-0.09780208, MSE=0.01277056 
## 
## Node number 39: 18 observations
##   mean=-0.06297328, MSE=0.01197143 
## 
## Node number 40: 18 observations
##   mean=0.01157027, MSE=0.01121365 
## 
## Node number 41: 17 observations
##   mean=0.09032636, MSE=0.01227181 
## 
## Node number 42: 11 observations
##   mean=0.1313815, MSE=0.01006617 
## 
## Node number 43: 20 observations,    complexity param=0.0003176373
##   mean=0.1996489, MSE=0.04172662 
##   left son=86 (13 obs) right son=87 (7 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.1363493    to the left,  improve=0.05776892, (0 missing)
## 
## Node number 48: 8 observations
##   mean=0.3736137, MSE=0.02161138 
## 
## Node number 49: 18 observations
##   mean=0.5119949, MSE=0.04457318 
## 
## Node number 50: 35 observations,    complexity param=0.001495507
##   mean=0.6072214, MSE=0.07330728 
##   left son=100 (24 obs) right son=101 (11 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.5104656    to the right, improve=0.07044408, (0 missing)
## 
## Node number 51: 23 observations,    complexity param=0.001557624
##   mean=0.7151142, MSE=0.07023273 
##   left son=102 (13 obs) right son=103 (10 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.7008803    to the right, improve=0.1467373, (0 missing)
## 
## Node number 52: 27 observations,    complexity param=0.0003274522
##   mean=0.8502948, MSE=0.07441528 
##   left son=104 (8 obs) right son=105 (19 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.9822524    to the right, improve=0.02085288, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91785     to the right, agree=0.741, adj=0.125, (0 split)
## 
## Node number 53: 18 observations
##   mean=0.960214, MSE=0.05671621 
## 
## Node number 74: 9 observations
##   mean=-0.1930877, MSE=0.003911875 
## 
## Node number 75: 25 observations,    complexity param=0.0001222121
##   mean=-0.1389078, MSE=0.0184945 
##   left son=150 (8 obs) right son=151 (17 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < -0.2185118   to the left,  improve=0.03382015, (0 missing)
## 
## Node number 86: 13 observations
##   mean=0.1636217, MSE=0.02154578 
## 
## Node number 87: 7 observations
##   mean=0.2665567, MSE=0.07231819 
## 
## Node number 100: 24 observations,    complexity param=0.001495507
##   mean=0.558571, MSE=0.05159449 
##   left son=200 (7 obs) right son=201 (17 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA < 0.5560231    to the left,  improve=0.1630998, (0 missing)
## 
## Node number 101: 11 observations
##   mean=0.7133678, MSE=0.1042495 
## 
## Node number 102: 13 observations
##   mean=0.6260776, MSE=0.06166083 
## 
## Node number 103: 10 observations
##   mean=0.8308617, MSE=0.05767296 
## 
## Node number 104: 8 observations
##   mean=0.7895868, MSE=0.03874549 
## 
## Node number 105: 19 observations
##   mean=0.8758561, MSE=0.08722899 
## 
## Node number 150: 8 observations
##   mean=-0.1753654, MSE=0.009977913 
## 
## Node number 151: 17 observations
##   mean=-0.1217513, MSE=0.02158247 
## 
## Node number 200: 7 observations
##   mean=0.4156146, MSE=0.01711153 
## 
## Node number 201: 17 observations
##   mean=0.6174355, MSE=0.05391327 
## 
## n= 417 
## 
## node), split, n, deviance, yval
##       * denotes terminal node
## 
##   1) root 417 127.95110000  0.34767020  
##     2) ILI.2.lag.log.nonNA< 0.3642121 238  17.23138000 -0.03522948  
##       4) ILI.2.lag.log.nonNA< -0.005288628 135   3.60432900 -0.20323710  
##         8) ILI.2.lag.log.nonNA< -0.2557664 46   0.80695440 -0.33687930  
##          16) ILI.2.lag.log.nonNA< -0.3364944 28   0.45961930 -0.38654060  
##            32) Week.bgn.last100.log>=8.958912 15   0.22676570 -0.41621790 *
##            33) Week.bgn.last100.log< 8.958912 13   0.20439890 -0.35229750 *
##          17) ILI.2.lag.log.nonNA>=-0.3364944 18   0.17086170 -0.25962840 *
##         9) ILI.2.lag.log.nonNA>=-0.2557664 89   1.55117300 -0.13416360  
##          18) ILI.2.lag.log.nonNA< -0.1187063 53   0.92007850 -0.17069070  
##            36) Week.bgn.last100.log< 8.958882 19   0.37423280 -0.20190130 *
##            37) Week.bgn.last100.log>=8.958882 34   0.51699520 -0.15324950  
##              74) ILI.2.lag.log.nonNA>=-0.1486436 9   0.03520688 -0.19308770 *
##              75) ILI.2.lag.log.nonNA< -0.1486436 25   0.46236240 -0.13890780  
##               150) ILI.2.lag.log.nonNA< -0.2185118 8   0.07982330 -0.17536540 *
##               151) ILI.2.lag.log.nonNA>=-0.2185118 17   0.36690200 -0.12175130 *
##          19) ILI.2.lag.log.nonNA>=-0.1187063 36   0.45627320 -0.08038768  
##            38) ILI.2.lag.log.nonNA< -0.06643646 18   0.22987000 -0.09780208 *
##            39) ILI.2.lag.log.nonNA>=-0.06643646 18   0.21548580 -0.06297328 *
##       5) ILI.2.lag.log.nonNA>=-0.005288628 103   4.82200400  0.18497470  
##        10) ILI.2.lag.log.nonNA< 0.1886472 66   1.70237300  0.10881800  
##          20) ILI.2.lag.log.nonNA< 0.1054685 35   0.46469410  0.04982323  
##            40) ILI.2.lag.log.nonNA< 0.05575614 18   0.20184570  0.01157027 *
##            41) ILI.2.lag.log.nonNA>=0.05575614 17   0.20862070  0.09032636 *
##          21) ILI.2.lag.log.nonNA>=0.1054685 31   0.97833440  0.17542500  
##            42) ILI.2.lag.log.nonNA>=0.1561068 11   0.11072790  0.13138150 *
##            43) ILI.2.lag.log.nonNA< 0.1561068 20   0.83453250  0.19964890  
##              86) ILI.2.lag.log.nonNA< 0.1363493 13   0.28009510  0.16362170 *
##              87) ILI.2.lag.log.nonNA>=0.1363493 7   0.50622730  0.26655670 *
##        11) ILI.2.lag.log.nonNA>=0.1886472 37   2.05402700  0.32082180  
##          22) ILI.2.lag.log.nonNA>=0.2672436 18   0.74553720  0.26198410 *
##          23) ILI.2.lag.log.nonNA< 0.2672436 19   1.18714200  0.37656280 *
##     3) ILI.2.lag.log.nonNA>=0.3642121 179  29.43105000  0.85677710  
##       6) ILI.2.lag.log.nonNA< 1.160636 141  13.88018000  0.72338830  
##        12) ILI.2.lag.log.nonNA< 0.7871896 84   6.00942000  0.59410950  
##          24) ILI.2.lag.log.nonNA< 0.4780125 26   1.08126600  0.46941610  
##            48) ILI.2.lag.log.nonNA>=0.4304485 8   0.17289100  0.37361370 *
##            49) ILI.2.lag.log.nonNA< 0.4304485 18   0.80231720  0.51199490 *
##          25) ILI.2.lag.log.nonNA>=0.4780125 58   4.34267400  0.65000650  
##            50) ILI.2.lag.log.nonNA< 0.6743843 35   2.56575500  0.60722140  
##             100) ILI.2.lag.log.nonNA>=0.5104656 24   1.23826800  0.55857100  
##               200) ILI.2.lag.log.nonNA< 0.5560231 7   0.11978070  0.41561460 *
##               201) ILI.2.lag.log.nonNA>=0.5560231 17   0.91652560  0.61743550 *
##             101) ILI.2.lag.log.nonNA< 0.5104656 11   1.14674500  0.71336780 *
##            51) ILI.2.lag.log.nonNA>=0.6743843 23   1.61535300  0.71511420  
##             102) ILI.2.lag.log.nonNA>=0.7008803 13   0.80159070  0.62607760 *
##             103) ILI.2.lag.log.nonNA< 0.7008803 10   0.57672960  0.83086170 *
##        13) ILI.2.lag.log.nonNA>=0.7871896 57   4.39797000  0.91390460  
##          26) ILI.2.lag.log.nonNA>=0.8769581 45   3.16059200  0.89426250  
##            52) ILI.2.lag.log.nonNA< 1.026878 27   2.00921300  0.85029480  
##             104) ILI.2.lag.log.nonNA>=0.9822524 8   0.30996400  0.78958680 *
##             105) ILI.2.lag.log.nonNA< 0.9822524 19   1.65735100  0.87585610 *
##            53) ILI.2.lag.log.nonNA>=1.026878 18   1.02089200  0.96021400 *
##          27) ILI.2.lag.log.nonNA< 0.8769581 12   1.15491000  0.98756240 *
##       7) ILI.2.lag.log.nonNA>=1.160636 38   3.73333200  1.35172000  
##        14) ILI.2.lag.log.nonNA< 1.485112 24   1.49722300  1.21473500  
##          28) ILI.2.lag.log.nonNA< 1.367642 14   0.69788730  1.12921000 *
##          29) ILI.2.lag.log.nonNA>=1.367642 10   0.55356860  1.33447000 *
##        15) ILI.2.lag.log.nonNA>=1.485112 14   1.01372400  1.58655000 *
## [1] "    calling mypredict_mdl for fit:"
## [1] "    calling mypredict_mdl for OOB:"
##                    model_id model_method
## 1 Max.cor.Y.cv.0.cp.0.rpart        rpart
##                                       feats max.nTuningRuns
## 1 ILI.2.lag.log.nonNA, Week.bgn.last100.log               0
##   min.elapsedtime.everything min.elapsedtime.final max.R.sq.fit
## 1                      0.476                 0.013    0.8742596
##   min.RMSE.fit max.R.sq.OOB min.RMSE.OOB
## 1    0.1964226    0.7442082     0.204778
```

```r
if (glb_is_regression || glb_is_binomial) # For multinomials this model will be run next by default
ret_lst <- myfit_mdl(model_id="Max.cor.Y", 
                        model_method="rpart",
                     model_type=glb_model_type,
                        indep_vars_vctr=max_cor_y_x_vars,
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                        n_cv_folds=glb_n_cv_folds, tune_models_df=NULL)
```

```
## [1] "fitting model: Max.cor.Y.rpart"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log"
```

```
## Warning in nominalTrainWorkflow(x = x, y = y, wts = weights, info =
## trainInfo, : There were missing values in resampled performance measures.
```

```
## Aggregating results
## Selecting tuning parameters
## Fitting cp = 0.0688 on full training set
```

```
## Warning in myfit_mdl(model_id = "Max.cor.Y", model_method = "rpart",
## model_type = glb_model_type, : model's bestTune found at an extreme of
## tuneGrid for parameter: cp
```

![](Google_Flu_template2_files/figure-html/fit.models_0-7.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-8.png) 

```
## Call:
## rpart(formula = .outcome ~ ., control = list(minsplit = 20, minbucket = 7, 
##     cp = 0, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
##     surrogatestyle = 0, maxdepth = 30, xval = 0))
##   n= 417 
## 
##          CP nsplit rel error
## 1 0.6353105      0 1.0000000
## 2 0.0923598      1 0.3646895
## 3 0.0688157      2 0.2723297
## 
## Variable importance
##  ILI.2.lag.log.nonNA Week.bgn.last100.log 
##                   91                    9 
## 
## Node number 1: 417 observations,    complexity param=0.6353105
##   mean=0.3476702, MSE=0.3068372 
##   left son=2 (238 obs) right son=3 (179 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 0.3642121 to the left,  improve=0.6353105, (0 missing)
##       Week.bgn.last100.log < 17.91785  to the left,  improve=0.0514081, (0 missing)
##   Surrogate splits:
##       Week.bgn.last100.log < 17.91785  to the left,  agree=0.619, adj=0.112, (0 split)
## 
## Node number 2: 238 observations
##   mean=-0.03522948, MSE=0.07240075 
## 
## Node number 3: 179 observations,    complexity param=0.0923598
##   mean=0.8567771, MSE=0.1644193 
##   left son=6 (141 obs) right son=7 (38 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA  < 1.160636  to the left,  improve=0.40153300, (0 missing)
##       Week.bgn.last100.log < 17.91779  to the right, improve=0.02044471, (0 missing)
## 
## Node number 6: 141 observations
##   mean=0.7233883, MSE=0.098441 
## 
## Node number 7: 38 observations
##   mean=1.35172, MSE=0.09824558 
## 
## n= 417 
## 
## node), split, n, deviance, yval
##       * denotes terminal node
## 
## 1) root 417 127.951100  0.34767020  
##   2) ILI.2.lag.log.nonNA< 0.3642121 238  17.231380 -0.03522948 *
##   3) ILI.2.lag.log.nonNA>=0.3642121 179  29.431050  0.85677710  
##     6) ILI.2.lag.log.nonNA< 1.160636 141  13.880180  0.72338830 *
##     7) ILI.2.lag.log.nonNA>=1.160636 38   3.733332  1.35172000 *
## [1] "    calling mypredict_mdl for fit:"
## [1] "    calling mypredict_mdl for OOB:"
##          model_id model_method                                     feats
## 1 Max.cor.Y.rpart        rpart ILI.2.lag.log.nonNA, Week.bgn.last100.log
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               3                      1.045                 0.014
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Rsquared.fit
## 1    0.7276703     0.305688    0.5715219    0.2650357        0.7000844
##   min.RMSESD.fit max.RsquaredSD.fit
## 1     0.03621604          0.0873017
```

```r
# Used to compare vs. Interactions.High.cor.Y and/or Max.cor.Y.TmSrs
ret_lst <- myfit_mdl(model_id="Max.cor.Y", 
                        model_method=ifelse(glb_is_regression, "lm", 
                                        ifelse(glb_is_binomial, "glm", "rpart")),
                     model_type=glb_model_type,
                        indep_vars_vctr=max_cor_y_x_vars,
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                        n_cv_folds=glb_n_cv_folds, tune_models_df=NULL)
```

```
## [1] "fitting model: Max.cor.Y.lm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_0-9.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-10.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-11.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-12.png) 

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.45573 -0.16134 -0.02313  0.11680  0.82385 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          0.013593   0.021911   0.620    0.535    
## ILI.2.lag.log.nonNA  0.917923   0.019611  46.806   <2e-16 ***
## Week.bgn.last100.log 0.001200   0.001419   0.846    0.398    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.2172 on 414 degrees of freedom
## Multiple R-squared:  0.8474,	Adjusted R-squared:  0.8467 
## F-statistic:  1150 on 2 and 414 DF,  p-value: < 2.2e-16
## 
## [1] "    calling mypredict_mdl for fit:"
## [1] "    calling mypredict_mdl for OOB:"
##       model_id model_method                                     feats
## 1 Max.cor.Y.lm           lm ILI.2.lag.log.nonNA, Week.bgn.last100.log
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      0.894                 0.002
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Adj.R.sq.fit
## 1    0.8474256    0.2172781    0.8008593    0.1806841        0.8466885
##   max.Rsquared.fit min.RMSESD.fit max.RsquaredSD.fit
## 1        0.8456678     0.01629213         0.03979421
```

```r
if (!is.null(glb_date_vars) && 
    (sum(grepl(paste(glb_date_vars, "\\.day\\.minutes\\.poly\\.", sep=""),
               names(glb_allobs_df))) > 0)) {
# ret_lst <- myfit_mdl(model_id="Max.cor.Y.TmSrs.poly1", 
#                         model_method=ifelse(glb_is_regression, "lm", 
#                                         ifelse(glb_is_binomial, "glm", "rpart")),
#                      model_type=glb_model_type,
#                         indep_vars_vctr=c(max_cor_y_x_vars, paste0(glb_date_vars, ".day.minutes")),
#                         rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
#                         fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
#                         n_cv_folds=glb_n_cv_folds, tune_models_df=NULL)
# 
ret_lst <- myfit_mdl(model_id="Max.cor.Y.TmSrs.poly", 
                        model_method=ifelse(glb_is_regression, "lm", 
                                        ifelse(glb_is_binomial, "glm", "rpart")),
                     model_type=glb_model_type,
                        indep_vars_vctr=c(max_cor_y_x_vars, 
            grep(paste(glb_date_vars, "\\.day\\.minutes\\.poly\\.", sep=""),
                        names(glb_allobs_df), value=TRUE)),
                        rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                        n_cv_folds=glb_n_cv_folds, tune_models_df=NULL)
}
```

```
## Warning in grepl(paste(glb_date_vars, "\\.day\\.minutes\\.poly\\.", sep =
## ""), : argument 'pattern' has length > 1 and only the first element will be
## used
```

```r
# Interactions.High.cor.Y
if (length(int_feats <- setdiff(unique(glb_feats_df$cor.high.X), NA)) > 0) {
    # lm & glm handle interaction terms; rpart & rf do not
    if (glb_is_regression || glb_is_binomial) {
        indep_vars_vctr <- 
            c(max_cor_y_x_vars, paste(max_cor_y_x_vars[1], int_feats, sep=":"))            
    } else { indep_vars_vctr <- union(max_cor_y_x_vars, int_feats) }
    
    ret_lst <- myfit_mdl(model_id="Interact.High.cor.Y", 
                            model_method=ifelse(glb_is_regression, "lm", 
                                        ifelse(glb_is_binomial, "glm", "rpart")),
                         model_type=glb_model_type,
                            indep_vars_vctr,
                            glb_rsp_var, glb_rsp_var_out,
                            fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                            n_cv_folds=glb_n_cv_folds, tune_models_df=NULL)                        
}    
```

```
## [1] "fitting model: Interact.High.cor.Y.lm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log, ILI.2.lag.log.nonNA:ILI.2.lag.log.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last100.log, ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last2.log, ILI.2.lag.log.nonNA:Week.bgn.last1.log, ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_0-13.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-14.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-15.png) 

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.49287 -0.12271 -0.00814  0.09731  0.93015 
## 
## Coefficients: (1 not defined because of singularities)
##                                                    Estimate Std. Error
## (Intercept)                                        0.018213   0.021403
## ILI.2.lag.log.nonNA                                2.024306   0.449141
## Week.bgn.last100.log                               0.002023   0.001367
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`         0.017678   0.010910
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`    0.148444   0.100875
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`   -0.211701   0.219974
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`   -0.230738   0.220369
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`   -0.177708   0.221369
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`   -0.260484   0.217276
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`   -0.243517   0.221416
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`   -0.234396   0.220642
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA9`          NA         NA
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`          -0.073575   0.045252
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`          -0.010988   0.057602
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`   0.025315   0.049117
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`  -0.241012   0.050936
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`  -0.483300   0.081461
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`   0.099225   0.098290
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`  -0.230973   0.097117
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`   0.091939   0.119895
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`   0.148733   0.112588
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`   0.083180   0.086770
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`  0.200677   0.078522
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11` -0.067658   0.068132
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`  0.163824   0.066647
##                                                   t value Pr(>|t|)    
## (Intercept)                                         0.851   0.3953    
## ILI.2.lag.log.nonNA                                 4.507 8.68e-06 ***
## Week.bgn.last100.log                                1.480   0.1397    
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`          1.620   0.1059    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`     1.472   0.1419    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`    -0.962   0.3364    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`    -1.047   0.2957    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`    -0.803   0.4226    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`    -1.199   0.2313    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`    -1.100   0.2721    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`    -1.062   0.2887    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA9`        NA       NA    
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`           -1.626   0.1048    
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`           -0.191   0.8488    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`    0.515   0.6066    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`   -4.732 3.11e-06 ***
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`   -5.933 6.54e-09 ***
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`    1.010   0.3133    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`   -2.378   0.0179 *  
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`    0.767   0.4436    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`    1.321   0.1873    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`    0.959   0.3383    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`   2.556   0.0110 *  
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`  -0.993   0.3213    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`   2.458   0.0144 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1906 on 393 degrees of freedom
## Multiple R-squared:  0.8884,	Adjusted R-squared:  0.8818 
## F-statistic:   136 on 23 and 393 DF,  p-value: < 2.2e-16
## 
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/fit.models_0-16.png) 

```
##                 model_id model_method
## 1 Interact.High.cor.Y.lm           lm
##                                                                                                                                                                                                                                                                                                       feats
## 1 ILI.2.lag.log.nonNA, Week.bgn.last100.log, ILI.2.lag.log.nonNA:ILI.2.lag.log.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last100.log, ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last2.log, ILI.2.lag.log.nonNA:Week.bgn.last1.log, ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      0.901                  0.01
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Adj.R.sq.fit
## 1    0.8883614     2.744383    0.8869087    0.1361616        0.8818278
##   max.Rsquared.fit min.RMSESD.fit max.RsquaredSD.fit
## 1        0.3980416       3.981253          0.4354334
```

```r
# Low.cor.X
# if (glb_is_classification && glb_is_binomial)
#     indep_vars_vctr <- subset(glb_feats_df, is.na(cor.high.X) & 
#                                             is.ConditionalX.y & 
#                                             (exclude.as.feat != 1))[, "id"] else
indep_vars_vctr <- subset(glb_feats_df, is.na(cor.high.X) & !myNearZV & 
                              (exclude.as.feat != 1))[, "id"]  
myadjust_interaction_feats <- function(vars_vctr) {
    for (feat in subset(glb_feats_df, !is.na(interaction.feat))$id)
        if (feat %in% vars_vctr)
            vars_vctr <- union(setdiff(vars_vctr, feat), 
                paste0(glb_feats_df[glb_feats_df$id == feat, "interaction.feat"], ":", feat))
    return(vars_vctr)
}
indep_vars_vctr <- myadjust_interaction_feats(indep_vars_vctr)
ret_lst <- myfit_mdl(model_id="Low.cor.X", 
                        model_method=ifelse(glb_is_regression, "lm", 
                                        ifelse(glb_is_binomial, "glm", "rpart")),
                        indep_vars_vctr=indep_vars_vctr,
                        model_type=glb_model_type,                     
                        glb_rsp_var, glb_rsp_var_out,
                        fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                        n_cv_folds=glb_n_cv_folds, tune_models_df=NULL)
```

```
## [1] "fitting model: Low.cor.X.lm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last2.log, Week.bgn.month.fctr.nonNA"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_0-17.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-18.png) ![](Google_Flu_template2_files/figure-html/fit.models_0-19.png) 

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.49176 -0.09908 -0.02105  0.06266  0.97631 
## 
## Coefficients: (5 not defined because of singularities)
##                               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                  0.4238511  0.1263135   3.356 0.000869 ***
## ILI.2.lag.log.nonNA          0.8670850  0.0245238  35.357  < 2e-16 ***
## Week.bgn.last100.log         0.0007402  0.0012222   0.606 0.545139    
## .rnorm                      -0.0042796  0.0093738  -0.457 0.648244    
## Week.bgn.date.fctr.nonNA2   -0.0083818  0.0265841  -0.315 0.752705    
## Week.bgn.date.fctr.nonNA3    0.0236838  0.0268024   0.884 0.377424    
## Week.bgn.date.fctr.nonNA4    0.0293814  0.0267301   1.099 0.272353    
## Week.bgn.date.fctr.nonNA5    0.0388797  0.0275316   1.412 0.158681    
## `Week.end.date.fctr(7,13]`          NA         NA      NA       NA    
## `Week.end.date.fctr(13,19]`         NA         NA      NA       NA    
## `Week.end.date.fctr(19,25]`         NA         NA      NA       NA    
## `Week.end.date.fctr(25,31]`         NA         NA      NA       NA    
## Week.bgn.last10.log          0.0166568  0.0044772   3.720 0.000228 ***
## Week.end.last10.log                 NA         NA      NA       NA    
## Week.bgn.last2.log          -0.0402170  0.0102724  -3.915 0.000106 ***
## Week.bgn.month.fctr.nonNA2   0.0526295  0.0441268   1.193 0.233706    
## Week.bgn.month.fctr.nonNA3  -0.2544180  0.0431501  -5.896 7.98e-09 ***
## Week.bgn.month.fctr.nonNA4  -0.2787947  0.0449952  -6.196 1.45e-09 ***
## Week.bgn.month.fctr.nonNA5  -0.1140584  0.0460397  -2.477 0.013651 *  
## Week.bgn.month.fctr.nonNA6  -0.2857019  0.0477678  -5.981 4.96e-09 ***
## Week.bgn.month.fctr.nonNA7  -0.2289393  0.0499140  -4.587 6.05e-06 ***
## Week.bgn.month.fctr.nonNA8  -0.0875292  0.0506637  -1.728 0.084830 .  
## Week.bgn.month.fctr.nonNA9  -0.0415082  0.0486449  -0.853 0.394014    
## Week.bgn.month.fctr.nonNA10  0.0129504  0.0451944   0.287 0.774607    
## Week.bgn.month.fctr.nonNA11 -0.0208782  0.0442255  -0.472 0.637126    
## Week.bgn.month.fctr.nonNA12  0.0966740  0.0437413   2.210 0.027667 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1772 on 396 degrees of freedom
## Multiple R-squared:  0.9028,	Adjusted R-squared:  0.8979 
## F-statistic: 183.9 on 20 and 396 DF,  p-value: < 2.2e-16
## 
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/fit.models_0-20.png) 

```
##       model_id model_method
## 1 Low.cor.X.lm           lm
##                                                                                                                                                                                      feats
## 1 ILI.2.lag.log.nonNA, Week.bgn.last100.log, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last2.log, Week.bgn.month.fctr.nonNA
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      0.905                  0.01
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Adj.R.sq.fit
## 1    0.9028104    0.1858368    0.8799352    0.1402968        0.8979018
##   max.Rsquared.fit min.RMSESD.fit max.RsquaredSD.fit
## 1        0.8905275     0.02000994         0.02855833
```

```r
rm(ret_lst)

glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor     bgn     end elapsed
## 10 fit.models          7          0 129.347 145.145  15.798
## 11 fit.models          7          1 145.145      NA      NA
```


```r
fit.models_1_chunk_df <- myadd_chunk(NULL, "fit.models_1_bgn")
```

```
##              label step_major step_minor     bgn end elapsed
## 1 fit.models_1_bgn          1          0 147.577  NA      NA
```

```r
# Options:
#   1. rpart & rf manual tuning
#   2. rf without pca (default: with pca)

#stop(here); sav_models_lst <- glb_models_lst; sav_models_df <- glb_models_df
#glb_models_lst <- sav_models_lst; glb_models_df <- sav_models_df

# All X that is not user excluded
# if (glb_is_classification && glb_is_binomial) {
#     model_id_pfx <- "Conditional.X"
# # indep_vars_vctr <- setdiff(names(glb_fitobs_df), union(glb_rsp_var, glb_exclude_vars_as_features))
#     indep_vars_vctr <- subset(glb_feats_df, is.ConditionalX.y & 
#                                             (exclude.as.feat != 1))[, "id"]
# } else {
    model_id_pfx <- "All.X"
    indep_vars_vctr <- subset(glb_feats_df, !myNearZV &
                                            (exclude.as.feat != 1))[, "id"]
# }

indep_vars_vctr <- myadjust_interaction_feats(indep_vars_vctr)

for (method in glb_models_method_vctr) {
    fit.models_1_chunk_df <- myadd_chunk(fit.models_1_chunk_df, 
                                paste0("fit.models_1_", method), major.inc=TRUE)
    if (method %in% c("rpart", "rf")) {
        # rpart:    fubar's the tree
        # rf:       skip the scenario w/ .rnorm for speed
        indep_vars_vctr <- setdiff(indep_vars_vctr, c(".rnorm"))
        model_id <- paste0(model_id_pfx, ".no.rnorm")
    } else model_id <- model_id_pfx
    
    ret_lst <- myfit_mdl(model_id=model_id, model_method=method,
                            indep_vars_vctr=indep_vars_vctr,
                            model_type=glb_model_type,
                            rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                            fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                n_cv_folds=glb_n_cv_folds, tune_models_df=glb_tune_models_df)
    
    # If All.X.glm is less accurate than Low.Cor.X.glm
    #   check NA coefficients & filter appropriate terms in indep_vars_vctr
#     if (method == "glm") {
#         orig_glm <- glb_models_lst[[paste0(model_id, ".", model_method)]]$finalModel
#         orig_glm <- glb_models_lst[["All.X.glm"]]$finalModel; print(summary(orig_glm))
#           vif_orig_glm <- vif(orig_glm); print(vif_orig_glm)
#           print(vif_orig_glm[!is.na(vif_orig_glm) & (vif_orig_glm == Inf)])
#           print(which.max(vif_orig_glm))
#           print(sort(vif_orig_glm[vif_orig_glm >= 1.0e+03], decreasing=TRUE))
#           glb_fitobs_df[c(1143, 3637, 3953, 4105), c("UniqueID", "Popular", "H.P.quandary", "Headline")]
#           glb_feats_df[glb_feats_df$id %in% grep("[HSA]\\.nchrs.log", glb_feats_df$id, value=TRUE) | glb_feats_df$cor.high.X %in%    grep("[HSA]\\.nchrs.log", glb_feats_df$id, value=TRUE), ]
#           glb_feats_df[glb_feats_df$id %in% grep("[HSA]\\.npnct14.log", glb_feats_df$id, value=TRUE) | glb_feats_df$cor.high.X %in%    grep("[HSA]\\.npnct14.log", glb_feats_df$id, value=TRUE), ]
#           glb_feats_df[glb_feats_df$id %in% grep("[HSA]\\.T.scen", glb_feats_df$id, value=TRUE) | glb_feats_df$cor.high.X %in%         grep("[HSA]\\.T.scen", glb_feats_df$id, value=TRUE), ]
#           glb_feats_df[glb_feats_df$id %in% grep("[HSA]\\.P.first", glb_feats_df$id, value=TRUE) | glb_feats_df$cor.high.X %in%         grep("[HSA]\\.P.first", glb_feats_df$id, value=TRUE), ]
#           all.equal(glb_allobs_df$S.nuppr.log, glb_allobs_df$A.nuppr.log)
#           all.equal(glb_allobs_df$S.npnct19.log, glb_allobs_df$A.npnct19.log)
#           all.equal(glb_allobs_df$S.P.year.colon, glb_allobs_df$A.P.year.colon)
#           all.equal(glb_allobs_df$S.T.share, glb_allobs_df$A.T.share)
#           all.equal(glb_allobs_df$H.T.clip, glb_allobs_df$H.P.daily.clip.report)
#           cor(glb_allobs_df$S.T.herald, glb_allobs_df$S.T.tribun)
#           dsp_obs(Abstract.contains="[Dd]iar", cols=("Abstract"), all=TRUE)
#           dsp_obs(Abstract.contains="[Ss]hare", cols=("Abstract"), all=TRUE)
#           subset(glb_feats_df, cor.y.abs <= glb_feats_df[glb_feats_df$id == ".rnorm", "cor.y.abs"])
#         corxx_mtrx <- cor(data.matrix(glb_allobs_df[, setdiff(names(glb_allobs_df), myfind_chr_cols_df(glb_allobs_df))]), use="pairwise.complete.obs"); abs_corxx_mtrx <- abs(corxx_mtrx); diag(abs_corxx_mtrx) <- 0
#           which.max(abs_corxx_mtrx["S.T.tribun", ])
#           abs_corxx_mtrx["A.npnct08.log", "S.npnct08.log"]
#         step_glm <- step(orig_glm)
#     }
    # Since caret does not optimize rpart well
#     if (method == "rpart")
#         ret_lst <- myfit_mdl(model_id=paste0(model_id_pfx, ".cp.0"), model_method=method,
#                                 indep_vars_vctr=indep_vars_vctr,
#                                 model_type=glb_model_type,
#                                 rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
#                                 fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,        
#             n_cv_folds=0, tune_models_df=data.frame(parameter="cp", min=0.0, max=0.0, by=0.1))
}
```

```
##              label step_major step_minor     bgn     end elapsed
## 1 fit.models_1_bgn          1          0 147.577 147.591   0.014
## 2  fit.models_1_lm          2          0 147.592      NA      NA
## [1] "fitting model: All.X.lm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_1-1.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-2.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-3.png) 

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.44670 -0.08553 -0.00315  0.06802  0.68025 
## 
## Coefficients: (28 not defined because of singularities)
##                              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                  0.307809   0.152075   2.024 0.043650 *  
## ILI.2.lag.log.nonNA          0.539252   0.033974  15.872  < 2e-16 ***
## Queries                      1.443489   0.126061  11.451  < 2e-16 ***
## Week.bgn.last100.log         0.007966   0.004655   1.711 0.087854 .  
## Week.end.last100.log               NA         NA      NA       NA    
## Week.bgn.year.fctr.nonNA2    0.006957   0.032882   0.212 0.832542    
## Week.bgn.year.fctr.nonNA3   -0.123304   0.089531  -1.377 0.169239    
## Week.bgn.year.fctr.nonNA4   -0.177282   0.089341  -1.984 0.047926 *  
## Week.bgn.year.fctr.nonNA5   -0.157908   0.089833  -1.758 0.079573 .  
## Week.bgn.year.fctr.nonNA6   -0.122828   0.093380  -1.315 0.189167    
## Week.bgn.year.fctr.nonNA7   -0.269157   0.090353  -2.979 0.003075 ** 
## Week.bgn.year.fctr.nonNA8   -0.304090   0.090975  -3.343 0.000911 ***
## Week.bgn.year.fctr.nonNA9          NA         NA      NA       NA    
## Week.end.year.fctr2005             NA         NA      NA       NA    
## Week.end.year.fctr2006             NA         NA      NA       NA    
## Week.end.year.fctr2007             NA         NA      NA       NA    
## Week.end.year.fctr2008             NA         NA      NA       NA    
## Week.end.year.fctr2009             NA         NA      NA       NA    
## Week.end.year.fctr2010             NA         NA      NA       NA    
## Week.end.year.fctr2011             NA         NA      NA       NA    
## Week.end.year.fctr2012             NA         NA      NA       NA    
## .rnorm                      -0.006227   0.008019  -0.777 0.437916    
## Week.bgn.date.fctr.nonNA2    0.008102   0.022702   0.357 0.721386    
## Week.bgn.date.fctr.nonNA3    0.020264   0.022760   0.890 0.373854    
## Week.bgn.date.fctr.nonNA4    0.032501   0.022695   1.432 0.152934    
## Week.bgn.date.fctr.nonNA5    0.020456   0.023449   0.872 0.383563    
## `Week.end.date.fctr(7,13]`         NA         NA      NA       NA    
## `Week.end.date.fctr(13,19]`        NA         NA      NA       NA    
## `Week.end.date.fctr(19,25]`        NA         NA      NA       NA    
## `Week.end.date.fctr(25,31]`        NA         NA      NA       NA    
## Week.bgn.last10.log          0.014859   0.003981   3.733 0.000218 ***
## Week.end.last10.log                NA         NA      NA       NA    
## Week.bgn.last1.log          -0.015451   0.016052  -0.963 0.336391    
## Week.end.last1.log                 NA         NA      NA       NA    
## Week.bgn.last2.log          -0.035227   0.011587  -3.040 0.002525 ** 
## Week.end.last2.log                 NA         NA      NA       NA    
## Week.bgn.month.fctr.nonNA2   0.123809   0.037870   3.269 0.001175 ** 
## Week.bgn.month.fctr.nonNA3  -0.050204   0.040735  -1.232 0.218520    
## Week.bgn.month.fctr.nonNA4  -0.136562   0.041911  -3.258 0.001219 ** 
## Week.bgn.month.fctr.nonNA5  -0.001851   0.043435  -0.043 0.966027    
## Week.bgn.month.fctr.nonNA6  -0.145283   0.045953  -3.162 0.001693 ** 
## Week.bgn.month.fctr.nonNA7  -0.167452   0.046920  -3.569 0.000404 ***
## Week.bgn.month.fctr.nonNA8  -0.082832   0.047103  -1.759 0.079451 .  
## Week.bgn.month.fctr.nonNA9  -0.113525   0.043241  -2.625 0.008997 ** 
## Week.bgn.month.fctr.nonNA10 -0.124632   0.040198  -3.100 0.002074 ** 
## Week.bgn.month.fctr.nonNA11 -0.086368   0.038130  -2.265 0.024058 *  
## Week.bgn.month.fctr.nonNA12 -0.058992   0.039944  -1.477 0.140526    
## Week.end.month.fctr02              NA         NA      NA       NA    
## Week.end.month.fctr03              NA         NA      NA       NA    
## Week.end.month.fctr04              NA         NA      NA       NA    
## Week.end.month.fctr05              NA         NA      NA       NA    
## Week.end.month.fctr06              NA         NA      NA       NA    
## Week.end.month.fctr07              NA         NA      NA       NA    
## Week.end.month.fctr08              NA         NA      NA       NA    
## Week.end.month.fctr09              NA         NA      NA       NA    
## Week.end.month.fctr10              NA         NA      NA       NA    
## Week.end.month.fctr11              NA         NA      NA       NA    
## Week.end.month.fctr12              NA         NA      NA       NA    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1502 on 387 degrees of freedom
## Multiple R-squared:  0.9318,	Adjusted R-squared:  0.9266 
## F-statistic: 182.2 on 29 and 387 DF,  p-value: < 2.2e-16
## 
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/fit.models_1-4.png) 

```
##   model_id model_method
## 1 All.X.lm           lm
##                                                                                                                                                                                                                                                                                                                                                    feats
## 1 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      0.977                 0.019
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Adj.R.sq.fit
## 1    0.9317513     3.454309     0.327414     0.332058        0.9266371
##   max.Rsquared.fit min.RMSESD.fit max.RsquaredSD.fit
## 1        0.3062242       2.969678          0.5263306
##              label step_major step_minor     bgn     end elapsed
## 2  fit.models_1_lm          2          0 147.592 150.085   2.493
## 3 fit.models_1_glm          3          0 150.086      NA      NA
## [1] "fitting model: All.X.glm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_1-5.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-6.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-7.png) 

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## 
## Call:
## NULL
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -0.44670  -0.08553  -0.00315   0.06802   0.68025  
## 
## Coefficients: (28 not defined because of singularities)
##                              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                  0.307809   0.152075   2.024 0.043650 *  
## ILI.2.lag.log.nonNA          0.539252   0.033974  15.872  < 2e-16 ***
## Queries                      1.443489   0.126061  11.451  < 2e-16 ***
## Week.bgn.last100.log         0.007966   0.004655   1.711 0.087854 .  
## Week.end.last100.log               NA         NA      NA       NA    
## Week.bgn.year.fctr.nonNA2    0.006957   0.032882   0.212 0.832542    
## Week.bgn.year.fctr.nonNA3   -0.123304   0.089531  -1.377 0.169239    
## Week.bgn.year.fctr.nonNA4   -0.177282   0.089341  -1.984 0.047926 *  
## Week.bgn.year.fctr.nonNA5   -0.157908   0.089833  -1.758 0.079573 .  
## Week.bgn.year.fctr.nonNA6   -0.122828   0.093380  -1.315 0.189167    
## Week.bgn.year.fctr.nonNA7   -0.269157   0.090353  -2.979 0.003075 ** 
## Week.bgn.year.fctr.nonNA8   -0.304090   0.090975  -3.343 0.000911 ***
## Week.bgn.year.fctr.nonNA9          NA         NA      NA       NA    
## Week.end.year.fctr2005             NA         NA      NA       NA    
## Week.end.year.fctr2006             NA         NA      NA       NA    
## Week.end.year.fctr2007             NA         NA      NA       NA    
## Week.end.year.fctr2008             NA         NA      NA       NA    
## Week.end.year.fctr2009             NA         NA      NA       NA    
## Week.end.year.fctr2010             NA         NA      NA       NA    
## Week.end.year.fctr2011             NA         NA      NA       NA    
## Week.end.year.fctr2012             NA         NA      NA       NA    
## .rnorm                      -0.006227   0.008019  -0.777 0.437916    
## Week.bgn.date.fctr.nonNA2    0.008102   0.022702   0.357 0.721386    
## Week.bgn.date.fctr.nonNA3    0.020264   0.022760   0.890 0.373854    
## Week.bgn.date.fctr.nonNA4    0.032501   0.022695   1.432 0.152934    
## Week.bgn.date.fctr.nonNA5    0.020456   0.023449   0.872 0.383563    
## `Week.end.date.fctr(7,13]`         NA         NA      NA       NA    
## `Week.end.date.fctr(13,19]`        NA         NA      NA       NA    
## `Week.end.date.fctr(19,25]`        NA         NA      NA       NA    
## `Week.end.date.fctr(25,31]`        NA         NA      NA       NA    
## Week.bgn.last10.log          0.014859   0.003981   3.733 0.000218 ***
## Week.end.last10.log                NA         NA      NA       NA    
## Week.bgn.last1.log          -0.015451   0.016052  -0.963 0.336391    
## Week.end.last1.log                 NA         NA      NA       NA    
## Week.bgn.last2.log          -0.035227   0.011587  -3.040 0.002525 ** 
## Week.end.last2.log                 NA         NA      NA       NA    
## Week.bgn.month.fctr.nonNA2   0.123809   0.037870   3.269 0.001175 ** 
## Week.bgn.month.fctr.nonNA3  -0.050204   0.040735  -1.232 0.218520    
## Week.bgn.month.fctr.nonNA4  -0.136562   0.041911  -3.258 0.001219 ** 
## Week.bgn.month.fctr.nonNA5  -0.001851   0.043435  -0.043 0.966027    
## Week.bgn.month.fctr.nonNA6  -0.145283   0.045953  -3.162 0.001693 ** 
## Week.bgn.month.fctr.nonNA7  -0.167452   0.046920  -3.569 0.000404 ***
## Week.bgn.month.fctr.nonNA8  -0.082832   0.047103  -1.759 0.079451 .  
## Week.bgn.month.fctr.nonNA9  -0.113525   0.043241  -2.625 0.008997 ** 
## Week.bgn.month.fctr.nonNA10 -0.124632   0.040198  -3.100 0.002074 ** 
## Week.bgn.month.fctr.nonNA11 -0.086368   0.038130  -2.265 0.024058 *  
## Week.bgn.month.fctr.nonNA12 -0.058992   0.039944  -1.477 0.140526    
## Week.end.month.fctr02              NA         NA      NA       NA    
## Week.end.month.fctr03              NA         NA      NA       NA    
## Week.end.month.fctr04              NA         NA      NA       NA    
## Week.end.month.fctr05              NA         NA      NA       NA    
## Week.end.month.fctr06              NA         NA      NA       NA    
## Week.end.month.fctr07              NA         NA      NA       NA    
## Week.end.month.fctr08              NA         NA      NA       NA    
## Week.end.month.fctr09              NA         NA      NA       NA    
## Week.end.month.fctr10              NA         NA      NA       NA    
## Week.end.month.fctr11              NA         NA      NA       NA    
## Week.end.month.fctr12              NA         NA      NA       NA    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 0.02256459)
## 
##     Null deviance: 127.9511  on 416  degrees of freedom
## Residual deviance:   8.7325  on 387  degrees of freedom
## AIC: -366.74
## 
## Number of Fisher Scoring iterations: 2
## 
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
##    model_id model_method
## 1 All.X.glm          glm
##                                                                                                                                                                                                                                                                                                                                                    feats
## 1 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                       1.03                 0.083
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB min.aic.fit
## 1    0.9317513     3.454309     0.327414     0.332058   -366.7419
##   max.Rsquared.fit min.RMSESD.fit max.RsquaredSD.fit
## 1        0.3062242       2.969678          0.5263306
##                   label step_major step_minor     bgn     end elapsed
## 3      fit.models_1_glm          3          0 150.086 152.705   2.619
## 4 fit.models_1_bayesglm          4          0 152.705      NA      NA
## [1] "fitting model: All.X.bayesglm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr"
```

```
## Loading required package: arm
## Loading required package: MASS
## 
## Attaching package: 'MASS'
## 
## The following object is masked from 'package:dplyr':
## 
##     select
## 
## Loading required package: Matrix
## Loading required package: lme4
## 
## arm (Version 1.8-5, built: 2015-05-13)
## 
## Working directory is /Users/bbalaji-2012/Documents/Work/Courses/MIT/Analytics_Edge_15_071x/Assignments/HW2_CDC_Google_Flu
```

```
## Warning in nominalTrainWorkflow(x = x, y = y, wts = weights, info =
## trainInfo, : There were missing values in resampled performance measures.
```

```
## Aggregating results
## Fitting final model on full training set
## 
## Call:
## NULL
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -0.44572  -0.08624  -0.00336   0.06790   0.68002  
## 
## Coefficients:
##                              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                  0.308556   0.157887   1.954   0.0514 .  
## ILI.2.lag.log.nonNA          0.539923   0.035243  15.320   <2e-16 ***
## Queries                      1.439142   0.130652  11.015   <2e-16 ***
## Week.bgn.last100.log         0.003961   1.662214   0.002   0.9981    
## Week.end.last100.log         0.003961   1.662214   0.002   0.9981    
## Week.bgn.year.fctr.nonNA2    0.003563   1.662403   0.002   0.9983    
## Week.bgn.year.fctr.nonNA3   -0.061232   1.663553  -0.037   0.9707    
## Week.bgn.year.fctr.nonNA4   -0.088169   1.663984  -0.053   0.9578    
## Week.bgn.year.fctr.nonNA5   -0.078468   1.663818  -0.047   0.9624    
## Week.bgn.year.fctr.nonNA6   -0.060749   1.663622  -0.037   0.9709    
## Week.bgn.year.fctr.nonNA7   -0.133938   1.665100  -0.080   0.9359    
## Week.bgn.year.fctr.nonNA8   -0.151354   1.665648  -0.091   0.9276    
## Week.bgn.year.fctr.nonNA9    0.000000   2.879119   0.000   1.0000    
## Week.end.year.fctr2005       0.003563   1.662403   0.002   0.9983    
## Week.end.year.fctr2006      -0.061232   1.663553  -0.037   0.9707    
## Week.end.year.fctr2007      -0.088169   1.663984  -0.053   0.9578    
## Week.end.year.fctr2008      -0.078468   1.663818  -0.047   0.9624    
## Week.end.year.fctr2009      -0.060749   1.663622  -0.037   0.9709    
## Week.end.year.fctr2010      -0.133938   1.665100  -0.080   0.9359    
## Week.end.year.fctr2011      -0.151354   1.665648  -0.091   0.9276    
## Week.end.year.fctr2012       0.000000   2.879119   0.000   1.0000    
## .rnorm                      -0.006215   0.008325  -0.747   0.4558    
## Week.bgn.date.fctr.nonNA2    0.004029   1.662343   0.002   0.9981    
## Week.bgn.date.fctr.nonNA3    0.010134   1.662352   0.006   0.9951    
## Week.bgn.date.fctr.nonNA4    0.016246   1.662370   0.010   0.9922    
## Week.bgn.date.fctr.nonNA5    0.010252   1.662356   0.006   0.9951    
## `Week.end.date.fctr(7,13]`   0.004029   1.662343   0.002   0.9981    
## `Week.end.date.fctr(13,19]`  0.010134   1.662352   0.006   0.9951    
## `Week.end.date.fctr(19,25]`  0.016246   1.662370   0.010   0.9922    
## `Week.end.date.fctr(25,31]`  0.010252   1.662356   0.006   0.9951    
## Week.bgn.last10.log          0.007432   1.662220   0.004   0.9964    
## Week.end.last10.log          0.007432   1.662220   0.004   0.9964    
## Week.bgn.last1.log          -0.007724   1.662288  -0.005   0.9963    
## Week.end.last1.log          -0.007724   1.662288  -0.005   0.9963    
## Week.bgn.last2.log          -0.017612   1.662278  -0.011   0.9916    
## Week.end.last2.log          -0.017612   1.662278  -0.011   0.9916    
## Week.bgn.month.fctr.nonNA2   0.061852   1.662854   0.037   0.9703    
## Week.bgn.month.fctr.nonNA3  -0.025372   1.662534  -0.015   0.9878    
## Week.bgn.month.fctr.nonNA4  -0.068533   1.662982  -0.041   0.9672    
## Week.bgn.month.fctr.nonNA5  -0.001180   1.662489  -0.001   0.9994    
## Week.bgn.month.fctr.nonNA6  -0.072942   1.663088  -0.044   0.9650    
## Week.bgn.month.fctr.nonNA7  -0.083944   1.663284  -0.050   0.9598    
## Week.bgn.month.fctr.nonNA8  -0.041565   1.662712  -0.025   0.9801    
## Week.bgn.month.fctr.nonNA9  -0.056751   1.662835  -0.034   0.9728    
## Week.bgn.month.fctr.nonNA10 -0.062169   1.662877  -0.037   0.9702    
## Week.bgn.month.fctr.nonNA11 -0.043112   1.662643  -0.026   0.9793    
## Week.bgn.month.fctr.nonNA12 -0.029243   1.662550  -0.018   0.9860    
## Week.end.month.fctr02        0.061852   1.662854   0.037   0.9703    
## Week.end.month.fctr03       -0.025372   1.662534  -0.015   0.9878    
## Week.end.month.fctr04       -0.068533   1.662982  -0.041   0.9672    
## Week.end.month.fctr05       -0.001180   1.662489  -0.001   0.9994    
## Week.end.month.fctr06       -0.072942   1.663088  -0.044   0.9650    
## Week.end.month.fctr07       -0.083944   1.663284  -0.050   0.9598    
## Week.end.month.fctr08       -0.041565   1.662712  -0.025   0.9801    
## Week.end.month.fctr09       -0.056751   1.662835  -0.034   0.9728    
## Week.end.month.fctr10       -0.062169   1.662877  -0.037   0.9702    
## Week.end.month.fctr11       -0.043112   1.662643  -0.026   0.9793    
## Week.end.month.fctr12       -0.029243   1.662550  -0.018   0.9860    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 0.02432458)
## 
##     Null deviance: 127.9511  on 416  degrees of freedom
## Residual deviance:   8.7325  on 359  degrees of freedom
## AIC: -310.74
## 
## Number of Fisher Scoring iterations: 8
## 
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
##         model_id model_method
## 1 All.X.bayesglm     bayesglm
##                                                                                                                                                                                                                                                                                                                                                    feats
## 1 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      1.669                 0.148
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB min.aic.fit
## 1    0.9317511    0.2484433    0.3333059    0.3306004   -310.7405
##   max.Rsquared.fit min.RMSESD.fit max.RsquaredSD.fit
## 1        0.7906456      0.1107067          0.1744511
##                   label step_major step_minor     bgn     end elapsed
## 4 fit.models_1_bayesglm          4          0 152.705 155.339   2.634
## 5    fit.models_1_rpart          5          0 155.340      NA      NA
## [1] "fitting model: All.X.no.rnorm.rpart"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr"
```

```
## Warning in nominalTrainWorkflow(x = x, y = y, wts = weights, info =
## trainInfo, : There were missing values in resampled performance measures.
```

![](Google_Flu_template2_files/figure-html/fit.models_1-8.png) 

```
## Aggregating results
## Selecting tuning parameters
## Fitting cp = 0.0688 on full training set
```

```
## Warning in myfit_mdl(model_id = model_id, model_method = method,
## indep_vars_vctr = indep_vars_vctr, : model's bestTune found at an extreme
## of tuneGrid for parameter: cp
```

![](Google_Flu_template2_files/figure-html/fit.models_1-9.png) 

```
## Call:
## rpart(formula = .outcome ~ ., control = list(minsplit = 20, minbucket = 7, 
##     cp = 0, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, 
##     surrogatestyle = 0, maxdepth = 30, xval = 0))
##   n= 417 
## 
##          CP nsplit rel error
## 1 0.6353105      0 1.0000000
## 2 0.0923598      1 0.3646895
## 3 0.0688157      2 0.2723297
## 
## Variable importance
##         ILI.2.lag.log.nonNA                     Queries 
##                          43                          22 
##         Week.bgn.last10.log         Week.end.last10.log 
##                           9                           9 
##   Week.bgn.year.fctr.nonNA6      Week.end.year.fctr2009 
##                           8                           8 
## Week.bgn.month.fctr.nonNA10       Week.end.month.fctr10 
##                           1                           1 
## 
## Node number 1: 417 observations,    complexity param=0.6353105
##   mean=0.3476702, MSE=0.3068372 
##   left son=2 (238 obs) right son=3 (179 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA        < 0.3642121 to the left,  improve=0.6353105, (0 missing)
##       Queries                    < 0.2848606 to the left,  improve=0.5265077, (0 missing)
##       Week.bgn.year.fctr.nonNA6  < 0.5       to the left,  improve=0.1640573, (0 missing)
##       Week.end.year.fctr2009     < 0.5       to the left,  improve=0.1640573, (0 missing)
##       Week.bgn.month.fctr.nonNA2 < 0.5       to the left,  improve=0.1356085, (0 missing)
##   Surrogate splits:
##       Queries                   < 0.3233732 to the left,  agree=0.808, adj=0.553, (0 split)
##       Week.bgn.last10.log       < 15.61554  to the left,  agree=0.674, adj=0.240, (0 split)
##       Week.end.last10.log       < 15.61554  to the left,  agree=0.674, adj=0.240, (0 split)
##       Week.bgn.year.fctr.nonNA6 < 0.5       to the left,  agree=0.657, adj=0.201, (0 split)
##       Week.end.year.fctr2009    < 0.5       to the left,  agree=0.657, adj=0.201, (0 split)
## 
## Node number 2: 238 observations
##   mean=-0.03522948, MSE=0.07240075 
## 
## Node number 3: 179 observations,    complexity param=0.0923598
##   mean=0.8567771, MSE=0.1644193 
##   left son=6 (141 obs) right son=7 (38 obs)
##   Primary splits:
##       ILI.2.lag.log.nonNA        < 1.160636  to the left,  improve=0.4015330, (0 missing)
##       Queries                    < 0.5053121 to the left,  improve=0.3302168, (0 missing)
##       Week.bgn.last10.log        < 15.61494  to the left,  improve=0.2225811, (0 missing)
##       Week.end.last10.log        < 15.61494  to the left,  improve=0.2225811, (0 missing)
##       Week.bgn.month.fctr.nonNA4 < 0.5       to the right, improve=0.2102918, (0 missing)
##   Surrogate splits:
##       Queries                     < 0.5471448 to the left,  agree=0.838, adj=0.237, (0 split)
##       Week.bgn.month.fctr.nonNA10 < 0.5       to the left,  agree=0.810, adj=0.105, (0 split)
##       Week.end.month.fctr10       < 0.5       to the left,  agree=0.810, adj=0.105, (0 split)
##       Week.bgn.month.fctr.nonNA9  < 0.5       to the left,  agree=0.799, adj=0.053, (0 split)
##       Week.end.month.fctr09       < 0.5       to the left,  agree=0.799, adj=0.053, (0 split)
## 
## Node number 6: 141 observations
##   mean=0.7233883, MSE=0.098441 
## 
## Node number 7: 38 observations
##   mean=1.35172, MSE=0.09824558 
## 
## n= 417 
## 
## node), split, n, deviance, yval
##       * denotes terminal node
## 
## 1) root 417 127.951100  0.34767020  
##   2) ILI.2.lag.log.nonNA< 0.3642121 238  17.231380 -0.03522948 *
##   3) ILI.2.lag.log.nonNA>=0.3642121 179  29.431050  0.85677710  
##     6) ILI.2.lag.log.nonNA< 1.160636 141  13.880180  0.72338830 *
##     7) ILI.2.lag.log.nonNA>=1.160636 38   3.733332  1.35172000 *
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
##               model_id model_method
## 1 All.X.no.rnorm.rpart        rpart
##                                                                                                                                                                                                                                                                                                                                            feats
## 1 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               3                      1.121                  0.08
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Rsquared.fit
## 1    0.7276703     0.305688    0.5715219    0.2650357        0.7000844
##   min.RMSESD.fit max.RsquaredSD.fit
## 1     0.03621604          0.0873017
##                label step_major step_minor     bgn     end elapsed
## 5 fit.models_1_rpart          5          0 155.340 158.458   3.118
## 6    fit.models_1_rf          6          0 158.458      NA      NA
## [1] "fitting model: All.X.no.rnorm.rf"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr"
```

```
## Loading required package: randomForest
## randomForest 4.6-10
## Type rfNews() to see new features/changes/bug fixes.
## 
## Attaching package: 'randomForest'
## 
## The following object is masked from 'package:dplyr':
## 
##     combine
```

![](Google_Flu_template2_files/figure-html/fit.models_1-10.png) 

```
## Aggregating results
## Selecting tuning parameters
## Fitting mtry = 29 on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_1-11.png) 

```
##                 Length Class      Mode     
## call              4    -none-     call     
## type              1    -none-     character
## predicted       417    -none-     numeric  
## mse             500    -none-     numeric  
## rsq             500    -none-     numeric  
## oob.times       417    -none-     numeric  
## importance       56    -none-     numeric  
## importanceSD      0    -none-     NULL     
## localImportance   0    -none-     NULL     
## proximity         0    -none-     NULL     
## ntree             1    -none-     numeric  
## mtry              1    -none-     numeric  
## forest           11    -none-     list     
## coefs             0    -none-     NULL     
## y               417    -none-     numeric  
## test              0    -none-     NULL     
## inbag             0    -none-     NULL     
## xNames           56    -none-     character
## problemType       1    -none-     character
## tuneValue         1    data.frame list     
## obsLevels         1    -none-     logical  
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## [1] "    calling mypredict_mdl for OOB:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.date.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

![](Google_Flu_template2_files/figure-html/fit.models_1-12.png) 

```
##            model_id model_method
## 1 All.X.no.rnorm.rf           rf
##                                                                                                                                                                                                                                                                                                                                            feats
## 1 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               3                      6.443                 1.857
##   max.R.sq.fit min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Rsquared.fit
## 1    0.9869949    0.1574949    0.7833803    0.1885728        0.9189164
##   min.RMSESD.fit max.RsquaredSD.fit
## 1     0.01511178         0.02617253
```

```r
# User specified
#   Ensure at least 2 vars in each regression; else varImp crashes
# sav_models_lst <- glb_models_lst; sav_models_df <- glb_models_df; sav_featsimp_df <- glb_featsimp_df
# glb_models_lst <- sav_models_lst; glb_models_df <- sav_models_df; glm_featsimp_df <- sav_featsimp_df

    # easier to exclude features
#model_id <- "";
# indep_vars_vctr <- head(subset(glb_models_df, grepl("All\\.X\\.", model_id), select=feats), 1)
# indep_vars_vctr <- setdiff(indep_vars_vctr, ".rnorm")

    # easier to include features
model_id <- "Flu.Trend2"; indep_vars_vctr <- c("Queries", "ILI.2.lag.log.nonNA")
for (method in c("lm")) {
    ret_lst <- myfit_mdl(model_id=model_id, model_method=method,
                                indep_vars_vctr=indep_vars_vctr,
                                model_type=glb_model_type,
                                rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
                                fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
                    n_cv_folds=glb_n_cv_folds, tune_models_df=glb_tune_models_df)
    csm_mdl_id <- paste0(model_id, ".", method)
    csm_featsimp_df <- myget_feats_importance(glb_models_lst[[paste0(model_id, ".", method)]]);         print(head(csm_featsimp_df))
}
```

```
## [1] "fitting model: Flu.Trend2.lm"
## [1] "    indep_vars: Queries, ILI.2.lag.log.nonNA"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.models_1-13.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-14.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-15.png) ![](Google_Flu_template2_files/figure-html/fit.models_1-16.png) 

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.52096 -0.11299 -0.01915  0.08230  0.76446 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         -0.23599    0.01978  -11.93   <2e-16 ***
## Queries              1.24372    0.08022   15.50   <2e-16 ***
## ILI.2.lag.log.nonNA  0.65847    0.02283   28.84   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1729 on 414 degrees of freedom
## Multiple R-squared:  0.9033,	Adjusted R-squared:  0.9028 
## F-statistic:  1934 on 2 and 414 DF,  p-value: < 2.2e-16
## 
## [1] "    calling mypredict_mdl for fit:"
## [1] "    calling mypredict_mdl for OOB:"
##        model_id model_method                        feats max.nTuningRuns
## 1 Flu.Trend2.lm           lm Queries, ILI.2.lag.log.nonNA               1
##   min.elapsedtime.everything min.elapsedtime.final max.R.sq.fit
## 1                      0.863                 0.003    0.9033037
##   min.RMSE.fit max.R.sq.OOB min.RMSE.OOB max.Adj.R.sq.fit max.Rsquared.fit
## 1    0.1729956    0.8565912    0.1533303        0.9028365        0.9020136
##   min.RMSESD.fit max.RsquaredSD.fit
## 1      0.0135587         0.02723253
##                     importance
## ILI.2.lag.log.nonNA        100
## Queries                      0
```

```r
# Ntv.1.lm <- lm(reformulate(indep_vars_vctr, glb_rsp_var), glb_trnobs_df); print(summary(Ntv.1.lm))

#print(dsp_models_df <- orderBy(model_sel_frmla, glb_models_df)[, dsp_models_cols])
#csm_featsimp_df[grepl("H.npnct19.log", row.names(csm_featsimp_df)), , FALSE]
#csm_OOBobs_df <- glb_get_predictions(glb_OOBobs_df, mdl_id=csm_mdl_id, rsp_var_out=glb_rsp_var_out, prob_threshold_def=glb_models_df[glb_models_df$model_id == csm_mdl_id, "opt.prob.threshold.OOB"])
#print(sprintf("%s OOB confusion matrix & accuracy: ", csm_mdl_id)); print(t(confusionMatrix(csm_OOBobs_df[, paste0(glb_rsp_var_out, csm_mdl_id)], csm_OOBobs_df[, glb_rsp_var])$table))

#glb_models_df[, "max.Accuracy.OOB", FALSE]
#varImp(glb_models_lst[["Low.cor.X.glm"]])
#orderBy(~ -Overall, varImp(glb_models_lst[["All.X.2.glm"]])$importance)
#orderBy(~ -Overall, varImp(glb_models_lst[["All.X.3.glm"]])$importance)
#glb_feats_df[grepl("npnct28", glb_feats_df$id), ]
#print(sprintf("%s OOB confusion matrix & accuracy: ", glb_sel_mdl_id)); print(t(confusionMatrix(glb_OOBobs_df[, paste0(glb_rsp_var_out, glb_sel_mdl_id)], glb_OOBobs_df[, glb_rsp_var])$table))

    # User specified bivariate models
#     indep_vars_vctr_lst <- list()
#     for (feat in setdiff(names(glb_fitobs_df), 
#                          union(glb_rsp_var, glb_exclude_vars_as_features)))
#         indep_vars_vctr_lst[["feat"]] <- feat

    # User specified combinatorial models
#     indep_vars_vctr_lst <- list()
#     combn_mtrx <- combn(c("<feat1_name>", "<feat2_name>", "<featn_name>"), 
#                           <num_feats_to_choose>)
#     for (combn_ix in 1:ncol(combn_mtrx))
#         #print(combn_mtrx[, combn_ix])
#         indep_vars_vctr_lst[[combn_ix]] <- combn_mtrx[, combn_ix]
    
    # template for myfit_mdl
    #   rf is hard-coded in caret to recognize only Accuracy / Kappa evaluation metrics
    #       only for OOB in trainControl ?
    
#     ret_lst <- myfit_mdl_fn(model_id=paste0(model_id_pfx, ""), model_method=method,
#                             indep_vars_vctr=indep_vars_vctr,
#                             rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out,
#                             fit_df=glb_fitobs_df, OOB_df=glb_OOBobs_df,
#                             n_cv_folds=glb_n_cv_folds, tune_models_df=glb_tune_models_df,
#                             model_loss_mtrx=glb_model_metric_terms,
#                             model_summaryFunction=glb_model_metric_smmry,
#                             model_metric=glb_model_metric,
#                             model_metric_maximize=glb_model_metric_maximize)

# Simplify a model
# fit_df <- glb_fitobs_df; glb_mdl <- step(<complex>_mdl)

# Non-caret models
#     rpart_area_mdl <- rpart(reformulate("Area", response=glb_rsp_var), 
#                                data=glb_fitobs_df, #method="class", 
#                                control=rpart.control(cp=0.12),
#                            parms=list(loss=glb_model_metric_terms))
#     print("rpart_sel_wlm_mdl"); prp(rpart_sel_wlm_mdl)
# 

print(glb_models_df)
```

```
##                                            model_id model_method
## MFO.lm                                       MFO.lm           lm
## Max.cor.Y.cv.0.rpart           Max.cor.Y.cv.0.rpart        rpart
## Max.cor.Y.cv.0.cp.0.rpart Max.cor.Y.cv.0.cp.0.rpart        rpart
## Max.cor.Y.rpart                     Max.cor.Y.rpart        rpart
## Max.cor.Y.lm                           Max.cor.Y.lm           lm
## Interact.High.cor.Y.lm       Interact.High.cor.Y.lm           lm
## Low.cor.X.lm                           Low.cor.X.lm           lm
## All.X.lm                                   All.X.lm           lm
## All.X.glm                                 All.X.glm          glm
## All.X.bayesglm                       All.X.bayesglm     bayesglm
## All.X.no.rnorm.rpart           All.X.no.rnorm.rpart        rpart
## All.X.no.rnorm.rf                 All.X.no.rnorm.rf           rf
## Flu.Trend2.lm                         Flu.Trend2.lm           lm
##                                                                                                                                                                                                                                                                                                                                                                            feats
## MFO.lm                                                                                                                                                                                                                                                                                                                                                                    .rnorm
## Max.cor.Y.cv.0.rpart                                                                                                                                                                                                                                                                                                                   ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Max.cor.Y.cv.0.cp.0.rpart                                                                                                                                                                                                                                                                                                              ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Max.cor.Y.rpart                                                                                                                                                                                                                                                                                                                        ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Max.cor.Y.lm                                                                                                                                                                                                                                                                                                                           ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Interact.High.cor.Y.lm                                                 ILI.2.lag.log.nonNA, Week.bgn.last100.log, ILI.2.lag.log.nonNA:ILI.2.lag.log.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last100.log, ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last2.log, ILI.2.lag.log.nonNA:Week.bgn.last1.log, ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA
## Low.cor.X.lm                                                                                                                                                                            ILI.2.lag.log.nonNA, Week.bgn.last100.log, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last2.log, Week.bgn.month.fctr.nonNA
## All.X.lm                  ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.glm                 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.bayesglm            ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.no.rnorm.rpart              ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.no.rnorm.rf                 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## Flu.Trend2.lm                                                                                                                                                                                                                                                                                                                                       Queries, ILI.2.lag.log.nonNA
##                           max.nTuningRuns min.elapsedtime.everything
## MFO.lm                                  0                      0.566
## Max.cor.Y.cv.0.rpart                    0                      0.547
## Max.cor.Y.cv.0.cp.0.rpart               0                      0.476
## Max.cor.Y.rpart                         3                      1.045
## Max.cor.Y.lm                            1                      0.894
## Interact.High.cor.Y.lm                  1                      0.901
## Low.cor.X.lm                            1                      0.905
## All.X.lm                                1                      0.977
## All.X.glm                               1                      1.030
## All.X.bayesglm                          1                      1.669
## All.X.no.rnorm.rpart                    3                      1.121
## All.X.no.rnorm.rf                       3                      6.443
## Flu.Trend2.lm                           1                      0.863
##                           min.elapsedtime.final max.R.sq.fit min.RMSE.fit
## MFO.lm                                    0.002  0.001602801    0.5534848
## Max.cor.Y.cv.0.rpart                      0.016  0.000000000    0.5539289
## Max.cor.Y.cv.0.cp.0.rpart                 0.013  0.874259639    0.1964226
## Max.cor.Y.rpart                           0.014  0.727670264    0.3056880
## Max.cor.Y.lm                              0.002  0.847425566    0.2172781
## Interact.High.cor.Y.lm                    0.010  0.888361359    2.7443835
## Low.cor.X.lm                              0.010  0.902810393    0.1858368
## All.X.lm                                  0.019  0.931751316    3.4543092
## All.X.glm                                 0.083  0.931751316    3.4543092
## All.X.bayesglm                            0.148  0.931751086    0.2484433
## All.X.no.rnorm.rpart                      0.080  0.727670264    0.3056880
## All.X.no.rnorm.rf                         1.857  0.986994930    0.1574949
## Flu.Trend2.lm                             0.003  0.903303681    0.1729956
##                           max.R.sq.OOB min.RMSE.OOB max.Adj.R.sq.fit
## MFO.lm                     -0.01575461    0.4080698    -0.0008029752
## Max.cor.Y.cv.0.rpart        0.00000000    0.4048928               NA
## Max.cor.Y.cv.0.cp.0.rpart   0.74420816    0.2047780               NA
## Max.cor.Y.rpart             0.57152193    0.2650357               NA
## Max.cor.Y.lm                0.80085932    0.1806841     0.8466884911
## Interact.High.cor.Y.lm      0.88690870    0.1361616     0.8818277999
## Low.cor.X.lm                0.87993525    0.1402968     0.8979018265
## All.X.lm                    0.32741401    0.3320580     0.9266370734
## All.X.glm                   0.32741401    0.3320580               NA
## All.X.bayesglm              0.33330588    0.3306004               NA
## All.X.no.rnorm.rpart        0.57152193    0.2650357               NA
## All.X.no.rnorm.rf           0.78338031    0.1885728               NA
## Flu.Trend2.lm               0.85659116    0.1533303     0.9028365494
##                           max.Rsquared.fit min.RMSESD.fit
## MFO.lm                                  NA             NA
## Max.cor.Y.cv.0.rpart                    NA             NA
## Max.cor.Y.cv.0.cp.0.rpart               NA             NA
## Max.cor.Y.rpart                  0.7000844     0.03621604
## Max.cor.Y.lm                     0.8456678     0.01629213
## Interact.High.cor.Y.lm           0.3980416     3.98125325
## Low.cor.X.lm                     0.8905275     0.02000994
## All.X.lm                         0.3062242     2.96967781
## All.X.glm                        0.3062242     2.96967781
## All.X.bayesglm                   0.7906456     0.11070673
## All.X.no.rnorm.rpart             0.7000844     0.03621604
## All.X.no.rnorm.rf                0.9189164     0.01511178
## Flu.Trend2.lm                    0.9020136     0.01355870
##                           max.RsquaredSD.fit min.aic.fit
## MFO.lm                                    NA          NA
## Max.cor.Y.cv.0.rpart                      NA          NA
## Max.cor.Y.cv.0.cp.0.rpart                 NA          NA
## Max.cor.Y.rpart                   0.08730170          NA
## Max.cor.Y.lm                      0.03979421          NA
## Interact.High.cor.Y.lm            0.43543342          NA
## Low.cor.X.lm                      0.02855833          NA
## All.X.lm                          0.52633057          NA
## All.X.glm                         0.52633057   -366.7419
## All.X.bayesglm                    0.17445114   -310.7405
## All.X.no.rnorm.rpart              0.08730170          NA
## All.X.no.rnorm.rf                 0.02617253          NA
## Flu.Trend2.lm                     0.02723253          NA
```

```r
rm(ret_lst)
fit.models_1_chunk_df <- myadd_chunk(fit.models_1_chunk_df, "fit.models_1_end", 
                                     major.inc=TRUE)
```

```
##              label step_major step_minor     bgn     end elapsed
## 6  fit.models_1_rf          6          0 158.458 168.966  10.508
## 7 fit.models_1_end          7          0 168.966      NA      NA
```

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor     bgn     end elapsed
## 11 fit.models          7          1 145.145 168.972  23.828
## 12 fit.models          7          2 168.973      NA      NA
```


```r
if (!is.null(glb_model_metric_smmry)) {
    stats_df <- glb_models_df[, "model_id", FALSE]

    stats_mdl_df <- data.frame()
    for (model_id in stats_df$model_id) {
        stats_mdl_df <- rbind(stats_mdl_df, 
            mypredict_mdl(glb_models_lst[[model_id]], glb_fitobs_df, glb_rsp_var, 
                          glb_rsp_var_out, model_id, "fit",
        						glb_model_metric_smmry, glb_model_metric, 
        						glb_model_metric_maximize, ret_type="stats"))
    }
    stats_df <- merge(stats_df, stats_mdl_df, all.x=TRUE)
    
    stats_mdl_df <- data.frame()
    for (model_id in stats_df$model_id) {
        stats_mdl_df <- rbind(stats_mdl_df, 
            mypredict_mdl(glb_models_lst[[model_id]], glb_OOBobs_df, glb_rsp_var, 
                          glb_rsp_var_out, model_id, "OOB",
            					glb_model_metric_smmry, glb_model_metric, 
        						glb_model_metric_maximize, ret_type="stats"))
    }
    stats_df <- merge(stats_df, stats_mdl_df, all.x=TRUE)
    
    print("Merging following data into glb_models_df:")
    print(stats_mrg_df <- stats_df[, c(1, grep(glb_model_metric, names(stats_df)))])
    print(tmp_models_df <- orderBy(~model_id, glb_models_df[, c("model_id",
                                    grep(glb_model_metric, names(stats_df), value=TRUE))]))

    tmp2_models_df <- glb_models_df[, c("model_id", setdiff(names(glb_models_df),
                                    grep(glb_model_metric, names(stats_df), value=TRUE)))]
    tmp3_models_df <- merge(tmp2_models_df, stats_mrg_df, all.x=TRUE, sort=FALSE)
    print(tmp3_models_df)
    print(names(tmp3_models_df))
    print(glb_models_df <- subset(tmp3_models_df, select=-model_id.1))
}

plt_models_df <- glb_models_df[, -grep("SD|Upper|Lower", names(glb_models_df))]
for (var in grep("^min.", names(plt_models_df), value=TRUE)) {
    plt_models_df[, sub("min.", "inv.", var)] <- 
        #ifelse(all(is.na(tmp <- plt_models_df[, var])), NA, 1.0 / tmp)
        1.0 / plt_models_df[, var]
    plt_models_df <- plt_models_df[ , -grep(var, names(plt_models_df))]
}
print(plt_models_df)
```

```
##                                            model_id model_method
## MFO.lm                                       MFO.lm           lm
## Max.cor.Y.cv.0.rpart           Max.cor.Y.cv.0.rpart        rpart
## Max.cor.Y.cv.0.cp.0.rpart Max.cor.Y.cv.0.cp.0.rpart        rpart
## Max.cor.Y.rpart                     Max.cor.Y.rpart        rpart
## Max.cor.Y.lm                           Max.cor.Y.lm           lm
## Interact.High.cor.Y.lm       Interact.High.cor.Y.lm           lm
## Low.cor.X.lm                           Low.cor.X.lm           lm
## All.X.lm                                   All.X.lm           lm
## All.X.glm                                 All.X.glm          glm
## All.X.bayesglm                       All.X.bayesglm     bayesglm
## All.X.no.rnorm.rpart           All.X.no.rnorm.rpart        rpart
## All.X.no.rnorm.rf                 All.X.no.rnorm.rf           rf
## Flu.Trend2.lm                         Flu.Trend2.lm           lm
##                                                                                                                                                                                                                                                                                                                                                                            feats
## MFO.lm                                                                                                                                                                                                                                                                                                                                                                    .rnorm
## Max.cor.Y.cv.0.rpart                                                                                                                                                                                                                                                                                                                   ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Max.cor.Y.cv.0.cp.0.rpart                                                                                                                                                                                                                                                                                                              ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Max.cor.Y.rpart                                                                                                                                                                                                                                                                                                                        ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Max.cor.Y.lm                                                                                                                                                                                                                                                                                                                           ILI.2.lag.log.nonNA, Week.bgn.last100.log
## Interact.High.cor.Y.lm                                                 ILI.2.lag.log.nonNA, Week.bgn.last100.log, ILI.2.lag.log.nonNA:ILI.2.lag.log.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last100.log, ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last2.log, ILI.2.lag.log.nonNA:Week.bgn.last1.log, ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA
## Low.cor.X.lm                                                                                                                                                                            ILI.2.lag.log.nonNA, Week.bgn.last100.log, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last2.log, Week.bgn.month.fctr.nonNA
## All.X.lm                  ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.glm                 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.bayesglm            ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, .rnorm, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.no.rnorm.rpart              ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## All.X.no.rnorm.rf                 ILI.2.lag.log.nonNA, Queries, Week.bgn.last100.log, Week.end.last100.log, Week.bgn.year.fctr.nonNA, Week.end.year.fctr, Week.bgn.date.fctr.nonNA, Week.end.date.fctr, Week.bgn.last10.log, Week.end.last10.log, Week.bgn.last1.log, Week.end.last1.log, Week.bgn.last2.log, Week.end.last2.log, Week.bgn.month.fctr.nonNA, Week.end.month.fctr
## Flu.Trend2.lm                                                                                                                                                                                                                                                                                                                                       Queries, ILI.2.lag.log.nonNA
##                           max.nTuningRuns max.R.sq.fit max.R.sq.OOB
## MFO.lm                                  0  0.001602801  -0.01575461
## Max.cor.Y.cv.0.rpart                    0  0.000000000   0.00000000
## Max.cor.Y.cv.0.cp.0.rpart               0  0.874259639   0.74420816
## Max.cor.Y.rpart                         3  0.727670264   0.57152193
## Max.cor.Y.lm                            1  0.847425566   0.80085932
## Interact.High.cor.Y.lm                  1  0.888361359   0.88690870
## Low.cor.X.lm                            1  0.902810393   0.87993525
## All.X.lm                                1  0.931751316   0.32741401
## All.X.glm                               1  0.931751316   0.32741401
## All.X.bayesglm                          1  0.931751086   0.33330588
## All.X.no.rnorm.rpart                    3  0.727670264   0.57152193
## All.X.no.rnorm.rf                       3  0.986994930   0.78338031
## Flu.Trend2.lm                           1  0.903303681   0.85659116
##                           max.Adj.R.sq.fit max.Rsquared.fit
## MFO.lm                       -0.0008029752               NA
## Max.cor.Y.cv.0.rpart                    NA               NA
## Max.cor.Y.cv.0.cp.0.rpart               NA               NA
## Max.cor.Y.rpart                         NA        0.7000844
## Max.cor.Y.lm                  0.8466884911        0.8456678
## Interact.High.cor.Y.lm        0.8818277999        0.3980416
## Low.cor.X.lm                  0.8979018265        0.8905275
## All.X.lm                      0.9266370734        0.3062242
## All.X.glm                               NA        0.3062242
## All.X.bayesglm                          NA        0.7906456
## All.X.no.rnorm.rpart                    NA        0.7000844
## All.X.no.rnorm.rf                       NA        0.9189164
## Flu.Trend2.lm                 0.9028365494        0.9020136
##                           inv.elapsedtime.everything inv.elapsedtime.final
## MFO.lm                                     1.7667845            500.000000
## Max.cor.Y.cv.0.rpart                       1.8281536             62.500000
## Max.cor.Y.cv.0.cp.0.rpart                  2.1008403             76.923077
## Max.cor.Y.rpart                            0.9569378             71.428571
## Max.cor.Y.lm                               1.1185682            500.000000
## Interact.High.cor.Y.lm                     1.1098779            100.000000
## Low.cor.X.lm                               1.1049724            100.000000
## All.X.lm                                   1.0235415             52.631579
## All.X.glm                                  0.9708738             12.048193
## All.X.bayesglm                             0.5991612              6.756757
## All.X.no.rnorm.rpart                       0.8920607             12.500000
## All.X.no.rnorm.rf                          0.1552072              0.538503
## Flu.Trend2.lm                              1.1587486            333.333333
##                           inv.RMSE.fit inv.RMSE.OOB  inv.aic.fit
## MFO.lm                       1.8067344     2.450561           NA
## Max.cor.Y.cv.0.rpart         1.8052859     2.469790           NA
## Max.cor.Y.cv.0.cp.0.rpart    5.0910651     4.883336           NA
## Max.cor.Y.rpart              3.2713097     3.773077           NA
## Max.cor.Y.lm                 4.6023969     5.534520           NA
## Interact.High.cor.Y.lm       0.3643806     7.344215           NA
## Low.cor.X.lm                 5.3810642     7.127746           NA
## All.X.lm                     0.2894935     3.011522           NA
## All.X.glm                    0.2894935     3.011522 -0.002726713
## All.X.bayesglm               4.0250626     3.024800 -0.003218119
## All.X.no.rnorm.rpart         3.2713097     3.773077           NA
## All.X.no.rnorm.rf            6.3494112     5.302993           NA
## Flu.Trend2.lm                5.7804924     6.521868           NA
```

```r
print(myplot_radar(radar_inp_df=plt_models_df))
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 13. Consider specifying shapes manually. if you must have them.
```

```
## Warning in loop_apply(n, do.ply): Removed 7 rows containing missing values
## (geom_path).
```

```
## Warning in loop_apply(n, do.ply): Removed 78 rows containing missing values
## (geom_point).
```

```
## Warning in loop_apply(n, do.ply): Removed 21 rows containing missing values
## (geom_text).
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 13. Consider specifying shapes manually. if you must have them.
```

![](Google_Flu_template2_files/figure-html/fit.models_2-1.png) 

```r
# print(myplot_radar(radar_inp_df=subset(plt_models_df, 
#         !(model_id %in% grep("random|MFO", plt_models_df$model_id, value=TRUE)))))

# Compute CI for <metric>SD
glb_models_df <- mutate(glb_models_df, 
                max.df = ifelse(max.nTuningRuns > 1, max.nTuningRuns - 1, NA),
                min.sd2ci.scaler = ifelse(is.na(max.df), NA, qt(0.975, max.df)))
for (var in grep("SD", names(glb_models_df), value=TRUE)) {
    # Does CI alredy exist ?
    var_components <- unlist(strsplit(var, "SD"))
    varActul <- paste0(var_components[1],          var_components[2])
    varUpper <- paste0(var_components[1], "Upper", var_components[2])
    varLower <- paste0(var_components[1], "Lower", var_components[2])
    if (varUpper %in% names(glb_models_df)) {
        warning(varUpper, " already exists in glb_models_df")
        # Assuming Lower also exists
        next
    }    
    print(sprintf("var:%s", var))
    # CI is dependent on sample size in t distribution; df=n-1
    glb_models_df[, varUpper] <- glb_models_df[, varActul] + 
        glb_models_df[, "min.sd2ci.scaler"] * glb_models_df[, var]
    glb_models_df[, varLower] <- glb_models_df[, varActul] - 
        glb_models_df[, "min.sd2ci.scaler"] * glb_models_df[, var]
}
```

```
## [1] "var:min.RMSESD.fit"
## [1] "var:max.RsquaredSD.fit"
```

```r
# Plot metrics with CI
plt_models_df <- glb_models_df[, "model_id", FALSE]
pltCI_models_df <- glb_models_df[, "model_id", FALSE]
for (var in grep("Upper", names(glb_models_df), value=TRUE)) {
    var_components <- unlist(strsplit(var, "Upper"))
    col_name <- unlist(paste(var_components, collapse=""))
    plt_models_df[, col_name] <- glb_models_df[, col_name]
    for (name in paste0(var_components[1], c("Upper", "Lower"), var_components[2]))
        pltCI_models_df[, name] <- glb_models_df[, name]
}

build_statsCI_data <- function(plt_models_df) {
    mltd_models_df <- melt(plt_models_df, id.vars="model_id")
    mltd_models_df$data <- sapply(1:nrow(mltd_models_df), 
        function(row_ix) tail(unlist(strsplit(as.character(
            mltd_models_df[row_ix, "variable"]), "[.]")), 1))
    mltd_models_df$label <- sapply(1:nrow(mltd_models_df), 
        function(row_ix) head(unlist(strsplit(as.character(
            mltd_models_df[row_ix, "variable"]), 
            paste0(".", mltd_models_df[row_ix, "data"]))), 1))
    #print(mltd_models_df)
    
    return(mltd_models_df)
}
mltd_models_df <- build_statsCI_data(plt_models_df)

mltdCI_models_df <- melt(pltCI_models_df, id.vars="model_id")
for (row_ix in 1:nrow(mltdCI_models_df)) {
    for (type in c("Upper", "Lower")) {
        if (length(var_components <- unlist(strsplit(
                as.character(mltdCI_models_df[row_ix, "variable"]), type))) > 1) {
            #print(sprintf("row_ix:%d; type:%s; ", row_ix, type))
            mltdCI_models_df[row_ix, "label"] <- var_components[1]
            mltdCI_models_df[row_ix, "data"] <- 
                unlist(strsplit(var_components[2], "[.]"))[2]
            mltdCI_models_df[row_ix, "type"] <- type
            break
        }
    }    
}
wideCI_models_df <- reshape(subset(mltdCI_models_df, select=-variable), 
                            timevar="type", 
        idvar=setdiff(names(mltdCI_models_df), c("type", "value", "variable")), 
                            direction="wide")
#print(wideCI_models_df)
mrgdCI_models_df <- merge(wideCI_models_df, mltd_models_df, all.x=TRUE)
#print(mrgdCI_models_df)

# Merge stats back in if CIs don't exist
goback_vars <- c()
for (var in unique(mltd_models_df$label)) {
    for (type in unique(mltd_models_df$data)) {
        var_type <- paste0(var, ".", type)
        # if this data is already present, next
        if (var_type %in% unique(paste(mltd_models_df$label, mltd_models_df$data,
                                       sep=".")))
            next
        #print(sprintf("var_type:%s", var_type))
        goback_vars <- c(goback_vars, var_type)
    }
}

if (length(goback_vars) > 0) {
    mltd_goback_df <- build_statsCI_data(glb_models_df[, c("model_id", goback_vars)])
    mltd_models_df <- rbind(mltd_models_df, mltd_goback_df)
}

mltd_models_df <- merge(mltd_models_df, glb_models_df[, c("model_id", "model_method")], 
                        all.x=TRUE)

png(paste0(glb_out_pfx, "models_bar.png"), width=480*3, height=480*2)
print(gp <- myplot_bar(mltd_models_df, "model_id", "value", colorcol_name="model_method") + 
        geom_errorbar(data=mrgdCI_models_df, 
            mapping=aes(x=model_id, ymax=value.Upper, ymin=value.Lower), width=0.5) + 
          facet_grid(label ~ data, scales="free") + 
          theme(axis.text.x = element_text(angle = 90,vjust = 0.5)))
```

```
## Warning in loop_apply(n, do.ply): Removed 3 rows containing missing values
## (position_stack).
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
print(gp)
```

```
## Warning in loop_apply(n, do.ply): Removed 3 rows containing missing values
## (position_stack).
```

![](Google_Flu_template2_files/figure-html/fit.models_2-2.png) 

```r
# used for console inspection
model_evl_terms <- c(NULL)
for (metric in glb_model_evl_criteria)
    model_evl_terms <- c(model_evl_terms, 
                         ifelse(length(grep("max", metric)) > 0, "-", "+"), metric)
if (glb_is_classification && glb_is_binomial)
    model_evl_terms <- c(model_evl_terms, "-", "opt.prob.threshold.OOB")
model_sel_frmla <- as.formula(paste(c("~ ", model_evl_terms), collapse=" "))
dsp_models_cols <- c("model_id", glb_model_evl_criteria) 
if (glb_is_classification && glb_is_binomial) 
    dsp_models_cols <- c(dsp_models_cols, "opt.prob.threshold.OOB")
print(dsp_models_df <- orderBy(model_sel_frmla, glb_models_df)[, dsp_models_cols])
```

```
##                     model_id min.RMSE.OOB max.R.sq.OOB max.Adj.R.sq.fit
## 6     Interact.High.cor.Y.lm    0.1361616   0.88690870     0.8818277999
## 7               Low.cor.X.lm    0.1402968   0.87993525     0.8979018265
## 13             Flu.Trend2.lm    0.1533303   0.85659116     0.9028365494
## 5               Max.cor.Y.lm    0.1806841   0.80085932     0.8466884911
## 12         All.X.no.rnorm.rf    0.1885728   0.78338031               NA
## 3  Max.cor.Y.cv.0.cp.0.rpart    0.2047780   0.74420816               NA
## 4            Max.cor.Y.rpart    0.2650357   0.57152193               NA
## 11      All.X.no.rnorm.rpart    0.2650357   0.57152193               NA
## 10            All.X.bayesglm    0.3306004   0.33330588               NA
## 8                   All.X.lm    0.3320580   0.32741401     0.9266370734
## 9                  All.X.glm    0.3320580   0.32741401               NA
## 2       Max.cor.Y.cv.0.rpart    0.4048928   0.00000000               NA
## 1                     MFO.lm    0.4080698  -0.01575461    -0.0008029752
```

```r
print(myplot_radar(radar_inp_df=dsp_models_df))
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 13. Consider specifying shapes manually. if you must have them.
```

```
## Warning in loop_apply(n, do.ply): Removed 5 rows containing missing values
## (geom_path).
```

```
## Warning in loop_apply(n, do.ply): Removed 29 rows containing missing values
## (geom_point).
```

```
## Warning in loop_apply(n, do.ply): Removed 7 rows containing missing values
## (geom_text).
```

```
## Warning in RColorBrewer::brewer.pal(n, pal): n too large, allowed maximum for palette Set1 is 9
## Returning the palette you asked for with that many colors
```

```
## Warning: The shape palette can deal with a maximum of 6 discrete values
## because more than 6 becomes difficult to discriminate; you have
## 13. Consider specifying shapes manually. if you must have them.
```

![](Google_Flu_template2_files/figure-html/fit.models_2-3.png) 

```r
print("Metrics used for model selection:"); print(model_sel_frmla)
```

```
## [1] "Metrics used for model selection:"
```

```
## ~+min.RMSE.OOB - max.R.sq.OOB - max.Adj.R.sq.fit
```

```r
print(sprintf("Best model id: %s", dsp_models_df[1, "model_id"]))
```

```
## [1] "Best model id: Interact.High.cor.Y.lm"
```

```r
if (is.null(glb_sel_mdl_id)) { 
    glb_sel_mdl_id <- dsp_models_df[1, "model_id"]
    if (glb_sel_mdl_id == "Interact.High.cor.Y.glm") {
        warning("glb_sel_mdl_id: Interact.High.cor.Y.glm; myextract_mdl_feats does not currently support interaction terms")
        glb_sel_mdl_id <- dsp_models_df[2, "model_id"]
    }
} else 
    print(sprintf("User specified selection: %s", glb_sel_mdl_id))   
    
myprint_mdl(glb_sel_mdl <- glb_models_lst[[glb_sel_mdl_id]])
```

![](Google_Flu_template2_files/figure-html/fit.models_2-4.png) ![](Google_Flu_template2_files/figure-html/fit.models_2-5.png) ![](Google_Flu_template2_files/figure-html/fit.models_2-6.png) 

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

![](Google_Flu_template2_files/figure-html/fit.models_2-7.png) 

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.49287 -0.12271 -0.00814  0.09731  0.93015 
## 
## Coefficients: (1 not defined because of singularities)
##                                                    Estimate Std. Error
## (Intercept)                                        0.018213   0.021403
## ILI.2.lag.log.nonNA                                2.024306   0.449141
## Week.bgn.last100.log                               0.002023   0.001367
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`         0.017678   0.010910
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`    0.148444   0.100875
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`   -0.211701   0.219974
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`   -0.230738   0.220369
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`   -0.177708   0.221369
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`   -0.260484   0.217276
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`   -0.243517   0.221416
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`   -0.234396   0.220642
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA9`          NA         NA
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`          -0.073575   0.045252
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`          -0.010988   0.057602
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`   0.025315   0.049117
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`  -0.241012   0.050936
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`  -0.483300   0.081461
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`   0.099225   0.098290
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`  -0.230973   0.097117
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`   0.091939   0.119895
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`   0.148733   0.112588
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`   0.083180   0.086770
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`  0.200677   0.078522
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11` -0.067658   0.068132
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`  0.163824   0.066647
##                                                   t value Pr(>|t|)    
## (Intercept)                                         0.851   0.3953    
## ILI.2.lag.log.nonNA                                 4.507 8.68e-06 ***
## Week.bgn.last100.log                                1.480   0.1397    
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`          1.620   0.1059    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`     1.472   0.1419    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`    -0.962   0.3364    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`    -1.047   0.2957    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`    -0.803   0.4226    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`    -1.199   0.2313    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`    -1.100   0.2721    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`    -1.062   0.2887    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA9`        NA       NA    
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`           -1.626   0.1048    
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`           -0.191   0.8488    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`    0.515   0.6066    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`   -4.732 3.11e-06 ***
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`   -5.933 6.54e-09 ***
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`    1.010   0.3133    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`   -2.378   0.0179 *  
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`    0.767   0.4436    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`    1.321   0.1873    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`    0.959   0.3383    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`   2.556   0.0110 *  
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`  -0.993   0.3213    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`   2.458   0.0144 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1906 on 393 degrees of freedom
## Multiple R-squared:  0.8884,	Adjusted R-squared:  0.8818 
## F-statistic:   136 on 23 and 393 DF,  p-value: < 2.2e-16
```

```
## [1] TRUE
```

```r
# From here to save(), this should all be in one function
#   these are executed in the same seq twice more:
#       fit.data.training & predict.data.new chunks
glb_get_predictions <- function(df, mdl_id, rsp_var_out, prob_threshold_def=NULL) {
    mdl <- glb_models_lst[[mdl_id]]
    rsp_var_out <- paste0(rsp_var_out, mdl_id)

    if (glb_is_regression) {
        df[, rsp_var_out] <- predict(mdl, newdata=df, type="raw")
        print(myplot_scatter(df, glb_rsp_var, rsp_var_out, smooth=TRUE))
        df[, paste0(rsp_var_out, ".err")] <- 
            abs(df[, rsp_var_out] - df[, glb_rsp_var])
        print(head(orderBy(reformulate(c("-", paste0(rsp_var_out, ".err"))), 
                           df)))                             
    }

    if (glb_is_classification && glb_is_binomial) {
        prob_threshold <- glb_models_df[glb_models_df$model_id == mdl_id, 
                                        "opt.prob.threshold.OOB"]
        if (is.null(prob_threshold) || is.na(prob_threshold)) {
            warning("Using default probability threshold: ", prob_threshold_def)
            if (is.null(prob_threshold <- prob_threshold_def))
                stop("Default probability threshold is NULL")
        }
        
        df[, paste0(rsp_var_out, ".prob")] <- 
            predict(mdl, newdata=df, type="prob")[, 2]
        df[, rsp_var_out] <- 
        		factor(levels(df[, glb_rsp_var])[
    				(df[, paste0(rsp_var_out, ".prob")] >=
    					prob_threshold) * 1 + 1], levels(df[, glb_rsp_var]))
    
        # prediction stats already reported by myfit_mdl ???
    }    
    
    if (glb_is_classification && !glb_is_binomial) {
        df[, rsp_var_out] <- predict(mdl, newdata=df, type="raw")
        df[, paste0(rsp_var_out, ".prob")] <- 
            predict(mdl, newdata=df, type="prob")
    }

    return(df)
}    
glb_OOBobs_df <- glb_get_predictions(df=glb_OOBobs_df, mdl_id=glb_sel_mdl_id, 
                                     rsp_var_out=glb_rsp_var_out)
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/fit.models_2-8.png) 

```
##                        Week       ILI   Queries .src     ILI.log
## 419 2012-01-08 - 2012-01-14 1.5434005 0.4993360 Test  0.43398812
## 430 2012-03-25 - 2012-03-31 1.7423860 0.3652058 Test  0.55525545
## 464 2012-11-18 - 2012-11-24 2.3046254 0.5112882 Test  0.83491815
## 445 2012-07-08 - 2012-07-14 0.9281519 0.2656042 Test -0.07455986
## 420 2012-01-15 - 2012-01-21 1.6476154 0.5006640 Test  0.49932902
## 441 2012-06-10 - 2012-06-16 1.0861211 0.2509960 Test  0.08261272
##         .rnorm   Week.bgn   Week.end ILI.2.lag ILI.2.lag.log
## 419 -1.3166927 2012-01-08 2012-01-14  2.124130    0.75336227
## 430 -0.2262124 2012-03-25 2012-03-31  2.293422    0.83004483
## 464 -0.4158885 2012-11-18 2012-11-24  1.610915    0.47680233
## 445  0.6476146 2012-07-08 2012-07-14  1.078713    0.07576853
## 420  0.8742224 2012-01-15 2012-01-21  1.766707    0.56911730
## 441 -1.6211072 2012-06-10 2012-06-16  1.299020    0.26160987
##     Week.bgn.POSIX Week.bgn.year.fctr Week.bgn.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
## 441     2012-06-10               2012                  06
##     Week.bgn.date.fctr Week.bgn.wkday.fctr Week.bgn.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
## 441             (7,13]                   0              1
##     Week.bgn.hour.fctr Week.bgn.minute.fctr Week.bgn.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
## 441                  0                    0                    0
##     Week.end.POSIX Week.end.year.fctr Week.end.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
## 441     2012-06-10               2012                  06
##     Week.end.date.fctr Week.end.wkday.fctr Week.end.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
## 441             (7,13]                   0              1
##     Week.end.hour.fctr Week.end.minute.fctr Week.end.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
## 441                  0                    0                    0
##     Week.bgn.zoo Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
## 419   1325394000           13.31265           14.00580            15.61583
## 430   1325998800           13.31265           14.00282            15.61464
## 464   1326603600           13.31265           14.00877            15.61583
## 445   1327208400           13.31265           14.00580            15.61524
## 420   1327813200           13.31265           14.00580            15.61583
## 441   1328418000           13.31265           14.00580            15.61524
##     Week.bgn.last100.log Week.end.zoo Week.end.last1.log
## 419             17.91782   1325394000           13.31265
## 430             17.91782   1325998800           13.31265
## 464             17.91782   1326603600           13.31265
## 445             17.91782   1327208400           13.31265
## 420             17.91782   1327813200           13.31265
## 441             17.91782   1328418000           13.31265
##     Week.end.last2.log Week.end.last10.log Week.end.last100.log
## 419           14.00580            15.61583             17.91782
## 430           14.00282            15.61464             17.91782
## 464           14.00877            15.61583             17.91782
## 445           14.00580            15.61524             17.91782
## 420           14.00580            15.61583             17.91782
## 441           14.00580            15.61524             17.91782
##     ILI.2.lag.log.nonNA Week.bgn.year.fctr.nonNA Week.bgn.month.fctr.nonNA
## 419          0.75336227                     2012                        01
## 430          0.83004483                     2012                        03
## 464          0.47680233                     2012                        11
## 445          0.07576853                     2012                        07
## 420          0.56911730                     2012                        01
## 441          0.26160987                     2012                        06
##     Week.bgn.date.fctr.nonNA ILI.log.predict.Interact.High.cor.Y.lm
## 419                   (7,13]                              0.9316064
## 430                  (19,25]                              0.8210195
## 464                  (13,19]                              0.5772443
## 445                   (7,13]                              0.1496508
## 420                  (13,19]                              0.7170902
## 441                   (7,13]                              0.2986348
##     ILI.log.predict.Interact.High.cor.Y.lm.err
## 419                                  0.4976183
## 430                                  0.2657641
## 464                                  0.2576738
## 445                                  0.2242107
## 420                                  0.2177612
## 441                                  0.2160221
```

```r
predct_accurate_var_name <- paste0(glb_rsp_var_out, glb_sel_mdl_id, ".accurate")
glb_OOBobs_df[, predct_accurate_var_name] <-
                    (glb_OOBobs_df[, glb_rsp_var] == 
                     glb_OOBobs_df[, paste0(glb_rsp_var_out, glb_sel_mdl_id)])

#stop(here"); #sav_models_lst <- glb_models_lst; sav_models_df <- glb_models_df
glb_featsimp_df <- 
    myget_feats_importance(mdl=glb_sel_mdl, featsimp_df=NULL)
glb_featsimp_df[, paste0(glb_sel_mdl_id, ".importance")] <- glb_featsimp_df$importance
print(glb_featsimp_df)
```

```
##                                                   importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`  100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`   79.081237
## ILI.2.lag.log.nonNA                                75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`  41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`  39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`   38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`           24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`         24.897758
## Week.bgn.last100.log                               22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`    22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`   19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`    17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`    15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`    15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`    14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`   14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`  13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`    13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`   13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`    10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`   10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`    5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`            0.000000
##                                                   Interact.High.cor.Y.lm.importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`                         100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`                          79.081237
## ILI.2.lag.log.nonNA                                                       75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`                         41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`                         39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`                          38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`                                  24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`                                24.897758
## Week.bgn.last100.log                                                      22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`                           22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`                          19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`                           17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`                           15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`                           15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`                           14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`                          14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`                         13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`                           13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`                          13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`                           10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`                          10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`                           5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`                                   0.000000
```

```r
# Used again in fit.data.training & predict.data.new chunks
glb_analytics_diag_plots <- function(obs_df, mdl_id, prob_threshold=NULL) {
    featsimp_df <- glb_featsimp_df
    featsimp_df$feat <- gsub("`(.*?)`", "\\1", row.names(featsimp_df))    
    featsimp_df$feat.interact <- gsub("(.*?):(.*)", "\\2", featsimp_df$feat)
    featsimp_df$feat <- gsub("(.*?):(.*)", "\\1", featsimp_df$feat)    
    featsimp_df$feat.interact <- ifelse(featsimp_df$feat.interact == featsimp_df$feat, 
                                        NA, featsimp_df$feat.interact)
    featsimp_df$feat <- gsub("(.*?)\\.fctr(.*)", "\\1\\.fctr", featsimp_df$feat)
    featsimp_df$feat.interact <- gsub("(.*?)\\.fctr(.*)", "\\1\\.fctr", featsimp_df$feat.interact) 
    featsimp_df <- orderBy(~ -importance.max, summaryBy(importance ~ feat + feat.interact, 
                                                        data=featsimp_df, FUN=max))    
    #rex_str=":(.*)"; txt_vctr=tail(featsimp_df$feat); ret_lst <- regexec(rex_str, txt_vctr); ret_lst <- regmatches(txt_vctr, ret_lst); ret_vctr <- sapply(1:length(ret_lst), function(pos_ix) ifelse(length(ret_lst[[pos_ix]]) > 0, ret_lst[[pos_ix]], "")); print(ret_vctr <- ret_vctr[ret_vctr != ""])    
    if (nrow(featsimp_df) > 5) {
        warning("Limiting important feature scatter plots to 5 out of ", nrow(featsimp_df))
        featsimp_df <- head(featsimp_df, 5)
    }
    
#     if (!all(is.na(featsimp_df$feat.interact)))
#         stop("not implemented yet")
    rsp_var_out <- paste0(glb_rsp_var_out, mdl_id)
    for (var in featsimp_df$feat) {
        plot_df <- melt(obs_df, id.vars=var, 
                        measure.vars=c(glb_rsp_var, rsp_var_out))

#         if (var == "<feat_name>") print(myplot_scatter(plot_df, var, "value", 
#                                              facet_colcol_name="variable") + 
#                       geom_vline(xintercept=<divider_val>, linetype="dotted")) else     
            print(myplot_scatter(plot_df, var, "value", colorcol_name="variable",
                                 facet_colcol_name="variable", jitter=TRUE) + 
                      guides(color=FALSE))
    }
    
    if (glb_is_regression) {
        if (nrow(featsimp_df) == 0)
            warning("No important features in glb_fin_mdl") else
            print(myplot_prediction_regression(df=obs_df, 
                        feat_x=ifelse(nrow(featsimp_df) > 1, featsimp_df$feat[2],
                                      ".rownames"), 
                                               feat_y=featsimp_df$feat[1],
                        rsp_var=glb_rsp_var, rsp_var_out=rsp_var_out,
                        id_vars=glb_id_var)
    #               + facet_wrap(reformulate(featsimp_df$feat[2])) # if [1 or 2] is a factor
    #               + geom_point(aes_string(color="<col_name>.fctr")) #  to color the plot
                  )
    }    
    
    if (glb_is_classification) {
        if (nrow(featsimp_df) == 0)
            warning("No features in selected model are statistically important")
        else print(myplot_prediction_classification(df=obs_df, 
                feat_x=ifelse(nrow(featsimp_df) > 1, featsimp_df$feat[2], 
                              ".rownames"),
                                               feat_y=featsimp_df$feat[1],
                     rsp_var=glb_rsp_var, 
                     rsp_var_out=rsp_var_out, 
                     id_vars=glb_id_var,
                    prob_threshold=prob_threshold)
#               + geom_hline(yintercept=<divider_val>, linetype = "dotted")
                )
    }    
}
if (glb_is_classification && glb_is_binomial)
    glb_analytics_diag_plots(obs_df=glb_OOBobs_df, mdl_id=glb_sel_mdl_id, 
            prob_threshold=glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, 
                                         "opt.prob.threshold.OOB"]) else
    glb_analytics_diag_plots(obs_df=glb_OOBobs_df, mdl_id=glb_sel_mdl_id)                  
```

```
## Warning in glb_analytics_diag_plots(obs_df = glb_OOBobs_df, mdl_id =
## glb_sel_mdl_id): Limiting important feature scatter plots to 5 out of 7
```

![](Google_Flu_template2_files/figure-html/fit.models_2-9.png) ![](Google_Flu_template2_files/figure-html/fit.models_2-10.png) ![](Google_Flu_template2_files/figure-html/fit.models_2-11.png) ![](Google_Flu_template2_files/figure-html/fit.models_2-12.png) ![](Google_Flu_template2_files/figure-html/fit.models_2-13.png) 

```
##                        Week       ILI   Queries .src     ILI.log
## 419 2012-01-08 - 2012-01-14 1.5434005 0.4993360 Test  0.43398812
## 430 2012-03-25 - 2012-03-31 1.7423860 0.3652058 Test  0.55525545
## 464 2012-11-18 - 2012-11-24 2.3046254 0.5112882 Test  0.83491815
## 445 2012-07-08 - 2012-07-14 0.9281519 0.2656042 Test -0.07455986
## 420 2012-01-15 - 2012-01-21 1.6476154 0.5006640 Test  0.49932902
##         .rnorm   Week.bgn   Week.end ILI.2.lag ILI.2.lag.log
## 419 -1.3166927 2012-01-08 2012-01-14  2.124130    0.75336227
## 430 -0.2262124 2012-03-25 2012-03-31  2.293422    0.83004483
## 464 -0.4158885 2012-11-18 2012-11-24  1.610915    0.47680233
## 445  0.6476146 2012-07-08 2012-07-14  1.078713    0.07576853
## 420  0.8742224 2012-01-15 2012-01-21  1.766707    0.56911730
##     Week.bgn.POSIX Week.bgn.year.fctr Week.bgn.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
##     Week.bgn.date.fctr Week.bgn.wkday.fctr Week.bgn.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
##     Week.bgn.hour.fctr Week.bgn.minute.fctr Week.bgn.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
##     Week.end.POSIX Week.end.year.fctr Week.end.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
##     Week.end.date.fctr Week.end.wkday.fctr Week.end.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
##     Week.end.hour.fctr Week.end.minute.fctr Week.end.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
##     Week.bgn.zoo Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
## 419   1325394000           13.31265           14.00580            15.61583
## 430   1325998800           13.31265           14.00282            15.61464
## 464   1326603600           13.31265           14.00877            15.61583
## 445   1327208400           13.31265           14.00580            15.61524
## 420   1327813200           13.31265           14.00580            15.61583
##     Week.bgn.last100.log Week.end.zoo Week.end.last1.log
## 419             17.91782   1325394000           13.31265
## 430             17.91782   1325998800           13.31265
## 464             17.91782   1326603600           13.31265
## 445             17.91782   1327208400           13.31265
## 420             17.91782   1327813200           13.31265
##     Week.end.last2.log Week.end.last10.log Week.end.last100.log
## 419           14.00580            15.61583             17.91782
## 430           14.00282            15.61464             17.91782
## 464           14.00877            15.61583             17.91782
## 445           14.00580            15.61524             17.91782
## 420           14.00580            15.61583             17.91782
##     ILI.2.lag.log.nonNA Week.bgn.year.fctr.nonNA Week.bgn.month.fctr.nonNA
## 419          0.75336227                     2012                        01
## 430          0.83004483                     2012                        03
## 464          0.47680233                     2012                        11
## 445          0.07576853                     2012                        07
## 420          0.56911730                     2012                        01
##     Week.bgn.date.fctr.nonNA ILI.log.predict.Interact.High.cor.Y.lm
## 419                   (7,13]                              0.9316064
## 430                  (19,25]                              0.8210195
## 464                  (13,19]                              0.5772443
## 445                   (7,13]                              0.1496508
## 420                  (13,19]                              0.7170902
##     ILI.log.predict.Interact.High.cor.Y.lm.err
## 419                                  0.4976183
## 430                                  0.2657641
## 464                                  0.2576738
## 445                                  0.2242107
## 420                                  0.2177612
##     ILI.log.predict.Interact.High.cor.Y.lm.accurate
## 419                                           FALSE
## 430                                           FALSE
## 464                                           FALSE
## 445                                           FALSE
## 420                                           FALSE
##                      .label
## 419 2012-01-08 - 2012-01-14
## 430 2012-03-25 - 2012-03-31
## 464 2012-11-18 - 2012-11-24
## 445 2012-07-08 - 2012-07-14
## 420 2012-01-15 - 2012-01-21
```

![](Google_Flu_template2_files/figure-html/fit.models_2-14.png) 

```r
# gather predictions from models better than MFO.*
#mdl_id <- "Conditional.X.rf"
#mdl_id <- "Conditional.X.cp.0.rpart"
#mdl_id <- "Conditional.X.rpart"
# glb_OOBobs_df <- glb_get_predictions(df=glb_OOBobs_df, mdl_id,
#                                      glb_rsp_var_out)
# print(t(confusionMatrix(glb_OOBobs_df[, paste0(glb_rsp_var_out, mdl_id)], 
#                         glb_OOBobs_df[, glb_rsp_var])$table))
# FN_OOB_ids <- c(4721, 4020, 693, 92)
# print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                     grep(glb_rsp_var, names(glb_OOBobs_df), value=TRUE)])
# print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                     glb_feats_df$id[1:5]])
# print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                     glb_txt_vars])
write.csv(glb_OOBobs_df[, c(glb_id_var, 
                grep(glb_rsp_var, names(glb_OOBobs_df), fixed=TRUE, value=TRUE))], 
    paste0(gsub(".", "_", paste0(glb_out_pfx, glb_sel_mdl_id), fixed=TRUE), 
           "_OOBobs.csv"), row.names=FALSE)

# print(glb_allobs_df[glb_allobs_df$UniqueID %in% FN_OOB_ids, 
#                     glb_txt_vars])
# dsp_tbl(Headline.contains="[Ee]bola")
# sum(sel_obs(Headline.contains="[Ee]bola"))
# ftable(xtabs(Popular ~ NewsDesk.fctr, data=glb_allobs_df[sel_obs(Headline.contains="[Ee]bola") ,]))
# xtabs(NewsDesk ~ Popular, #Popular ~ NewsDesk.fctr, 
#       data=glb_allobs_df[sel_obs(Headline.contains="[Ee]bola") ,],
#       exclude=NULL)
# print(mycreate_xtab_df(df=glb_allobs_df[sel_obs(Headline.contains="[Ee]bola") ,], c("Popular", "NewsDesk", "SectionName", "SubsectionName")))
# print(mycreate_tbl_df(df=glb_allobs_df[sel_obs(Headline.contains="[Ee]bola") ,], c("Popular", "NewsDesk", "SectionName", "SubsectionName")))
# print(mycreate_tbl_df(df=glb_allobs_df[sel_obs(Headline.contains="[Ee]bola") ,], c("Popular")))
# print(mycreate_tbl_df(df=glb_allobs_df[sel_obs(Headline.contains="[Ee]bola") ,], 
#                       tbl_col_names=c("Popular", "NewsDesk")))

# write.csv(glb_chunks_df, paste0(glb_out_pfx, tail(glb_chunks_df, 1)$label, "_",
#                                 tail(glb_chunks_df, 1)$step_minor,  "_chunks1.csv"),
#           row.names=FALSE)

glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.models", major.inc=FALSE)
```

```
##         label step_major step_minor     bgn     end elapsed
## 12 fit.models          7          2 168.973 180.142  11.169
## 13 fit.models          7          3 180.142      NA      NA
```


```r
print(setdiff(names(glb_trnobs_df), names(glb_allobs_df)))
```

```
##  [1] "Week.bgn.wkday.fctr"  "Week.bgn.wkend"       "Week.bgn.hour.fctr"  
##  [4] "Week.bgn.minute.fctr" "Week.bgn.second.fctr" "Week.end.wkday.fctr" 
##  [7] "Week.end.wkend"       "Week.end.hour.fctr"   "Week.end.minute.fctr"
## [10] "Week.end.second.fctr"
```

```r
print(setdiff(names(glb_fitobs_df), names(glb_allobs_df)))
```

```
##  [1] "Week.bgn.wkday.fctr"  "Week.bgn.wkend"       "Week.bgn.hour.fctr"  
##  [4] "Week.bgn.minute.fctr" "Week.bgn.second.fctr" "Week.end.wkday.fctr" 
##  [7] "Week.end.wkend"       "Week.end.hour.fctr"   "Week.end.minute.fctr"
## [10] "Week.end.second.fctr"
```

```r
print(setdiff(names(glb_OOBobs_df), names(glb_allobs_df)))
```

```
##  [1] "Week.bgn.wkday.fctr"                            
##  [2] "Week.bgn.wkend"                                 
##  [3] "Week.bgn.hour.fctr"                             
##  [4] "Week.bgn.minute.fctr"                           
##  [5] "Week.bgn.second.fctr"                           
##  [6] "Week.end.wkday.fctr"                            
##  [7] "Week.end.wkend"                                 
##  [8] "Week.end.hour.fctr"                             
##  [9] "Week.end.minute.fctr"                           
## [10] "Week.end.second.fctr"                           
## [11] "ILI.log.predict.Interact.High.cor.Y.lm"         
## [12] "ILI.log.predict.Interact.High.cor.Y.lm.err"     
## [13] "ILI.log.predict.Interact.High.cor.Y.lm.accurate"
```

```r
for (col in setdiff(names(glb_OOBobs_df), names(glb_allobs_df)))
    # Merge or cbind ?
    glb_allobs_df[glb_allobs_df$.lcn == "OOB", col] <- glb_OOBobs_df[, col]
    
print(setdiff(names(glb_newobs_df), names(glb_allobs_df)))
```

```
## character(0)
```

```r
if (glb_save_envir)
    save(glb_feats_df, 
         glb_allobs_df, #glb_trnobs_df, glb_fitobs_df, glb_OOBobs_df, glb_newobs_df,
         glb_models_df, dsp_models_df, glb_models_lst, glb_sel_mdl, glb_sel_mdl_id,
         glb_model_type,
        file=paste0(glb_out_pfx, "selmdl_dsk.RData"))
#load(paste0(glb_out_pfx, "selmdl_dsk.RData"))

rm(ret_lst)
```

```
## Warning in rm(ret_lst): object 'ret_lst' not found
```

```r
replay.petrisim(pn=glb_analytics_pn, 
    replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
        "model.selected")), flip_coord=TRUE)
```

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0 
## 2.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction data.new.prediction 	firing:  model.selected 
## 3.0000 	 3 	 0 2 1 0
```

![](Google_Flu_template2_files/figure-html/fit.models_3-1.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=TRUE)
```

```
##                label step_major step_minor     bgn     end elapsed
## 13        fit.models          7          3 180.142 184.107   3.965
## 14 fit.data.training          8          0 184.108      NA      NA
```

## Step `8.0: fit data training`

```r
#load(paste0(glb_inp_pfx, "dsk.RData"))

# To create specific models
# glb_fin_mdl_id <- NULL; glb_fin_mdl <- NULL; 
# glb_sel_mdl_id <- "Conditional.X.cp.0.rpart"; 
# glb_sel_mdl <- glb_models_lst[[glb_sel_mdl_id]]; print(glb_sel_mdl)
    
if (!is.null(glb_fin_mdl_id) && (glb_fin_mdl_id %in% names(glb_models_lst))) {
    warning("Final model same as user selected model")
    glb_fin_mdl <- glb_sel_mdl
} else {    
#     print(mdl_feats_df <- myextract_mdl_feats(sel_mdl=glb_sel_mdl, 
#                                               entity_df=glb_fitobs_df))
    
    if ((model_method <- glb_sel_mdl$method) == "custom")
        # get actual method from the model_id
        model_method <- tail(unlist(strsplit(glb_sel_mdl_id, "[.]")), 1)
        
    tune_finmdl_df <- NULL
    if (nrow(glb_sel_mdl$bestTune) > 0) {
        for (param in names(glb_sel_mdl$bestTune)) {
            #print(sprintf("param: %s", param))
            if (glb_sel_mdl$bestTune[1, param] != "none")
                tune_finmdl_df <- rbind(tune_finmdl_df, 
                    data.frame(parameter=param, 
                               min=glb_sel_mdl$bestTune[1, param], 
                               max=glb_sel_mdl$bestTune[1, param], 
                               by=1)) # by val does not matter
        }
    } 
    
    # Sync with parameters in mydsutils.R
    require(gdata)
    ret_lst <- myfit_mdl(model_id="Final", model_method=model_method,
        indep_vars_vctr=trim(unlist(strsplit(glb_models_df[glb_models_df$model_id == glb_sel_mdl_id,
                                                    "feats"], "[,]"))), 
                         model_type=glb_model_type,
                            rsp_var=glb_rsp_var, rsp_var_out=glb_rsp_var_out, 
                            fit_df=glb_trnobs_df, OOB_df=NULL,
                            n_cv_folds=glb_n_cv_folds, tune_models_df=tune_finmdl_df,
                         # Automate from here
                         #  Issues if glb_sel_mdl$method == "rf" b/c trainControl is "oob"; not "cv"
                            model_loss_mtrx=glb_model_metric_terms,
                            model_summaryFunction=glb_sel_mdl$control$summaryFunction,
                            model_metric=glb_sel_mdl$metric,
                            model_metric_maximize=glb_sel_mdl$maximize)
    glb_fin_mdl <- glb_models_lst[[length(glb_models_lst)]] 
    glb_fin_mdl_id <- glb_models_df[length(glb_models_lst), "model_id"]
}
```

```
## Loading required package: gdata
## gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.
## 
## gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.
## 
## Attaching package: 'gdata'
## 
## The following object is masked from 'package:randomForest':
## 
##     combine
## 
## The following objects are masked from 'package:dplyr':
## 
##     combine, first, last
## 
## The following object is masked from 'package:stats':
## 
##     nobs
## 
## The following object is masked from 'package:utils':
## 
##     object.size
```

```
## [1] "fitting model: Final.lm"
## [1] "    indep_vars: ILI.2.lag.log.nonNA, Week.bgn.last100.log, ILI.2.lag.log.nonNA:ILI.2.lag.log.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last100.log, ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last2.log, ILI.2.lag.log.nonNA:Week.bgn.last1.log, ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA"
## Aggregating results
## Fitting final model on full training set
```

![](Google_Flu_template2_files/figure-html/fit.data.training_0-1.png) ![](Google_Flu_template2_files/figure-html/fit.data.training_0-2.png) ![](Google_Flu_template2_files/figure-html/fit.data.training_0-3.png) 

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
```

```
## 
## Call:
## lm(formula = .outcome ~ ., data = dat)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.49287 -0.12271 -0.00814  0.09731  0.93015 
## 
## Coefficients: (1 not defined because of singularities)
##                                                    Estimate Std. Error
## (Intercept)                                        0.018213   0.021403
## ILI.2.lag.log.nonNA                                2.024306   0.449141
## Week.bgn.last100.log                               0.002023   0.001367
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`         0.017678   0.010910
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`    0.148444   0.100875
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`   -0.211701   0.219974
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`   -0.230738   0.220369
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`   -0.177708   0.221369
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`   -0.260484   0.217276
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`   -0.243517   0.221416
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`   -0.234396   0.220642
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA9`          NA         NA
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`          -0.073575   0.045252
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`          -0.010988   0.057602
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`   0.025315   0.049117
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`  -0.241012   0.050936
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`  -0.483300   0.081461
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`   0.099225   0.098290
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`  -0.230973   0.097117
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`   0.091939   0.119895
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`   0.148733   0.112588
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`   0.083180   0.086770
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`  0.200677   0.078522
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11` -0.067658   0.068132
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`  0.163824   0.066647
##                                                   t value Pr(>|t|)    
## (Intercept)                                         0.851   0.3953    
## ILI.2.lag.log.nonNA                                 4.507 8.68e-06 ***
## Week.bgn.last100.log                                1.480   0.1397    
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`          1.620   0.1059    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`     1.472   0.1419    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`    -0.962   0.3364    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`    -1.047   0.2957    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`    -0.803   0.4226    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`    -1.199   0.2313    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`    -1.100   0.2721    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`    -1.062   0.2887    
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA9`        NA       NA    
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`           -1.626   0.1048    
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`           -0.191   0.8488    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`    0.515   0.6066    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`   -4.732 3.11e-06 ***
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`   -5.933 6.54e-09 ***
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`    1.010   0.3133    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`   -2.378   0.0179 *  
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`    0.767   0.4436    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`    1.321   0.1873    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`    0.959   0.3383    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`   2.556   0.0110 *  
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`  -0.993   0.3213    
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`   2.458   0.0144 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1906 on 393 degrees of freedom
## Multiple R-squared:  0.8884,	Adjusted R-squared:  0.8818 
## F-statistic:   136 on 23 and 393 DF,  p-value: < 2.2e-16
## 
## [1] "    calling mypredict_mdl for fit:"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/fit.data.training_0-4.png) 

```
##   model_id model_method
## 1 Final.lm           lm
##                                                                                                                                                                                                                                                                                                       feats
## 1 ILI.2.lag.log.nonNA, Week.bgn.last100.log, ILI.2.lag.log.nonNA:ILI.2.lag.log.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last100.log, ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA, ILI.2.lag.log.nonNA:Week.bgn.last2.log, ILI.2.lag.log.nonNA:Week.bgn.last1.log, ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA
##   max.nTuningRuns min.elapsedtime.everything min.elapsedtime.final
## 1               1                      0.922                  0.01
##   max.R.sq.fit min.RMSE.fit max.Adj.R.sq.fit max.Rsquared.fit
## 1    0.8883614     2.744383        0.8818278        0.3980416
##   min.RMSESD.fit max.RsquaredSD.fit
## 1       3.981253          0.4354334
```

```r
rm(ret_lst)
glb_chunks_df <- myadd_chunk(glb_chunks_df, "fit.data.training", major.inc=FALSE)
```

```
##                label step_major step_minor     bgn     end elapsed
## 14 fit.data.training          8          0 184.108 188.576   4.468
## 15 fit.data.training          8          1 188.577      NA      NA
```


```r
glb_trnobs_df <- glb_get_predictions(df=glb_trnobs_df, mdl_id=glb_fin_mdl_id, 
                                     rsp_var_out=glb_rsp_var_out,
    prob_threshold_def=ifelse(glb_is_classification && glb_is_binomial, 
        glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, "opt.prob.threshold.OOB"], NULL))
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/fit.data.training_1-1.png) 

```
##                        Week      ILI   Queries  .src   ILI.log      .rnorm
## 278 2009-04-26 - 2009-05-02 2.981589 0.5553810 Train 1.0924564  0.05420224
## 296 2009-08-30 - 2009-09-05 3.719694 0.4116866 Train 1.3136413 -1.70337580
## 295 2009-08-23 - 2009-08-29 2.471660 0.3466135 Train 0.9048899 -0.83932764
## 279 2009-05-03 - 2009-05-09 2.437022 0.5551129 Train 0.8907770  0.76701532
## 213 2008-01-27 - 2008-02-02 4.433810 0.4143426 Train 1.4892593  1.83754480
## 282 2009-05-24 - 2009-05-30 4.213152 0.2948207 Train 1.4382111  0.58885757
##       Week.bgn   Week.end ILI.2.lag ILI.2.lag.log Week.bgn.POSIX
## 278 2009-04-26 2009-05-02  1.292327     0.2564443     2009-04-26
## 296 2009-08-30 2009-09-05  1.641071     0.4953493     2009-08-30
## 295 2009-08-23 2009-08-29  1.161419     0.1496424     2009-08-23
## 279 2009-05-03 2009-05-09  1.271641     0.2403083     2009-05-03
## 213 2008-01-27 2008-02-02  2.359343     0.8583831     2008-01-27
## 282 2009-05-24 2009-05-30  2.281301     0.8247459     2009-05-24
##     Week.bgn.year.fctr Week.bgn.month.fctr Week.bgn.date.fctr
## 278               2009                  04            (25,31]
## 296               2009                  08            (25,31]
## 295               2009                  08            (19,25]
## 279               2009                  05           (0.97,7]
## 213               2008                  01            (25,31]
## 282               2009                  05            (19,25]
##     Week.bgn.wkday.fctr Week.bgn.wkend Week.bgn.hour.fctr
## 278                   0              1                  0
## 296                   0              1                  0
## 295                   0              1                  0
## 279                   0              1                  0
## 213                   0              1                  0
## 282                   0              1                  0
##     Week.bgn.minute.fctr Week.bgn.second.fctr Week.end.POSIX
## 278                    0                    0     2009-04-26
## 296                    0                    0     2009-08-30
## 295                    0                    0     2009-08-23
## 279                    0                    0     2009-05-03
## 213                    0                    0     2008-01-27
## 282                    0                    0     2009-05-24
##     Week.end.year.fctr Week.end.month.fctr Week.end.date.fctr
## 278               2009                  04            (25,31]
## 296               2009                  08            (25,31]
## 295               2009                  08            (19,25]
## 279               2009                  05           (0.97,7]
## 213               2008                  01            (25,31]
## 282               2009                  05            (19,25]
##     Week.end.wkday.fctr Week.end.wkend Week.end.hour.fctr
## 278                   0              1                  0
## 296                   0              1                  0
## 295                   0              1                  0
## 279                   0              1                  0
## 213                   0              1                  0
## 282                   0              1                  0
##     Week.end.minute.fctr Week.end.second.fctr Week.bgn.zoo
## 278                    0                    0   1073192400
## 296                    0                    0   1073797200
## 295                    0                    0   1074402000
## 279                    0                    0   1075006800
## 213                    0                    0   1075611600
## 282                    0                    0   1076216400
##     Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
## 278           13.31265            14.0058            15.61464
## 296           13.31265            14.0058            15.61524
## 295           13.31265            14.0058            15.61524
## 279           13.31265            14.0058            15.61464
## 213           13.31265            14.0058            15.61524
## 282           13.31265            14.0058            15.61524
##     Week.bgn.last100.log Week.end.zoo Week.end.last1.log
## 278             17.91782   1073192400           13.31265
## 296             17.91782   1073797200           13.31265
## 295             17.91782   1074402000           13.31265
## 279             17.91782   1075006800           13.31265
## 213             17.91782   1075611600           13.31265
## 282             17.91782   1076216400           13.31265
##     Week.end.last2.log Week.end.last10.log Week.end.last100.log
## 278            14.0058            15.61464             17.91782
## 296            14.0058            15.61524             17.91782
## 295            14.0058            15.61524             17.91782
## 279            14.0058            15.61464             17.91782
## 213            14.0058            15.61524             17.91782
## 282            14.0058            15.61524             17.91782
##     ILI.2.lag.log.nonNA Week.bgn.year.fctr.nonNA Week.bgn.month.fctr.nonNA
## 278           0.2564443                     2009                        04
## 296           0.4953493                     2009                        08
## 295           0.1496424                     2009                        08
## 279           0.2403083                     2009                        05
## 213           0.8583831                     2008                        01
## 282           0.8247459                     2009                        05
##     Week.bgn.date.fctr.nonNA ILI.log.predict.Final.lm
## 278                  (25,31]                0.1623063
## 296                  (25,31]                0.5758465
## 295                  (19,25]                0.2119733
## 279                 (0.97,7]                0.2955065
## 213                  (25,31]                0.9013407
## 282                  (19,25]                0.8817212
##     ILI.log.predict.Final.lm.err
## 278                    0.9301501
## 296                    0.7377949
## 295                    0.6929166
## 279                    0.5952705
## 213                    0.5879186
## 282                    0.5564899
```

```r
sav_featsimp_df <- glb_featsimp_df
#glb_feats_df <- sav_feats_df
# glb_feats_df <- mymerge_feats_importance(feats_df=glb_feats_df, sel_mdl=glb_fin_mdl, 
#                                                entity_df=glb_trnobs_df)
glb_featsimp_df <- myget_feats_importance(mdl=glb_fin_mdl, featsimp_df=glb_featsimp_df)
glb_featsimp_df[, paste0(glb_fin_mdl_id, ".importance")] <- glb_featsimp_df$importance
print(glb_featsimp_df)
```

```
##                                                   Interact.High.cor.Y.lm.importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`                         100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`                          79.081237
## ILI.2.lag.log.nonNA                                                       75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`                         41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`                         39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`                          38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`                                  24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`                                24.897758
## Week.bgn.last100.log                                                      22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`                           22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`                          19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`                           17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`                           15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`                           15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`                           14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`                          14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`                         13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`                           13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`                          13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`                           10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`                          10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`                           5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`                                   0.000000
##                                                   importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`  100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`   79.081237
## ILI.2.lag.log.nonNA                                75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`  41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`  39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`   38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`           24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`         24.897758
## Week.bgn.last100.log                               22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`    22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`   19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`    17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`    15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`    15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`    14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`   14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`  13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`    13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`   13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`    10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`   10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`    5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`            0.000000
##                                                   Final.lm.importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`           100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`            79.081237
## ILI.2.lag.log.nonNA                                         75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`           41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`           39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`            38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`                    24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`                  24.897758
## Week.bgn.last100.log                                        22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`             22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`            19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`             17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`             15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`             15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`             14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`            14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`           13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`             13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`            13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`             10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`            10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`             5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`                     0.000000
```

```r
if (glb_is_classification && glb_is_binomial)
    glb_analytics_diag_plots(obs_df=glb_trnobs_df, mdl_id=glb_fin_mdl_id, 
            prob_threshold=glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, 
                                         "opt.prob.threshold.OOB"]) else
    glb_analytics_diag_plots(obs_df=glb_trnobs_df, mdl_id=glb_fin_mdl_id)                  
```

```
## Warning in glb_analytics_diag_plots(obs_df = glb_trnobs_df, mdl_id =
## glb_fin_mdl_id): Limiting important feature scatter plots to 5 out of 7
```

![](Google_Flu_template2_files/figure-html/fit.data.training_1-2.png) ![](Google_Flu_template2_files/figure-html/fit.data.training_1-3.png) ![](Google_Flu_template2_files/figure-html/fit.data.training_1-4.png) ![](Google_Flu_template2_files/figure-html/fit.data.training_1-5.png) ![](Google_Flu_template2_files/figure-html/fit.data.training_1-6.png) 

```
##                        Week      ILI   Queries  .src   ILI.log      .rnorm
## 278 2009-04-26 - 2009-05-02 2.981589 0.5553810 Train 1.0924564  0.05420224
## 296 2009-08-30 - 2009-09-05 3.719694 0.4116866 Train 1.3136413 -1.70337580
## 295 2009-08-23 - 2009-08-29 2.471660 0.3466135 Train 0.9048899 -0.83932764
## 279 2009-05-03 - 2009-05-09 2.437022 0.5551129 Train 0.8907770  0.76701532
## 213 2008-01-27 - 2008-02-02 4.433810 0.4143426 Train 1.4892593  1.83754480
##       Week.bgn   Week.end ILI.2.lag ILI.2.lag.log Week.bgn.POSIX
## 278 2009-04-26 2009-05-02  1.292327     0.2564443     2009-04-26
## 296 2009-08-30 2009-09-05  1.641071     0.4953493     2009-08-30
## 295 2009-08-23 2009-08-29  1.161419     0.1496424     2009-08-23
## 279 2009-05-03 2009-05-09  1.271641     0.2403083     2009-05-03
## 213 2008-01-27 2008-02-02  2.359343     0.8583831     2008-01-27
##     Week.bgn.year.fctr Week.bgn.month.fctr Week.bgn.date.fctr
## 278               2009                  04            (25,31]
## 296               2009                  08            (25,31]
## 295               2009                  08            (19,25]
## 279               2009                  05           (0.97,7]
## 213               2008                  01            (25,31]
##     Week.bgn.wkday.fctr Week.bgn.wkend Week.bgn.hour.fctr
## 278                   0              1                  0
## 296                   0              1                  0
## 295                   0              1                  0
## 279                   0              1                  0
## 213                   0              1                  0
##     Week.bgn.minute.fctr Week.bgn.second.fctr Week.end.POSIX
## 278                    0                    0     2009-04-26
## 296                    0                    0     2009-08-30
## 295                    0                    0     2009-08-23
## 279                    0                    0     2009-05-03
## 213                    0                    0     2008-01-27
##     Week.end.year.fctr Week.end.month.fctr Week.end.date.fctr
## 278               2009                  04            (25,31]
## 296               2009                  08            (25,31]
## 295               2009                  08            (19,25]
## 279               2009                  05           (0.97,7]
## 213               2008                  01            (25,31]
##     Week.end.wkday.fctr Week.end.wkend Week.end.hour.fctr
## 278                   0              1                  0
## 296                   0              1                  0
## 295                   0              1                  0
## 279                   0              1                  0
## 213                   0              1                  0
##     Week.end.minute.fctr Week.end.second.fctr Week.bgn.zoo
## 278                    0                    0   1073192400
## 296                    0                    0   1073797200
## 295                    0                    0   1074402000
## 279                    0                    0   1075006800
## 213                    0                    0   1075611600
##     Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
## 278           13.31265            14.0058            15.61464
## 296           13.31265            14.0058            15.61524
## 295           13.31265            14.0058            15.61524
## 279           13.31265            14.0058            15.61464
## 213           13.31265            14.0058            15.61524
##     Week.bgn.last100.log Week.end.zoo Week.end.last1.log
## 278             17.91782   1073192400           13.31265
## 296             17.91782   1073797200           13.31265
## 295             17.91782   1074402000           13.31265
## 279             17.91782   1075006800           13.31265
## 213             17.91782   1075611600           13.31265
##     Week.end.last2.log Week.end.last10.log Week.end.last100.log
## 278            14.0058            15.61464             17.91782
## 296            14.0058            15.61524             17.91782
## 295            14.0058            15.61524             17.91782
## 279            14.0058            15.61464             17.91782
## 213            14.0058            15.61524             17.91782
##     ILI.2.lag.log.nonNA Week.bgn.year.fctr.nonNA Week.bgn.month.fctr.nonNA
## 278           0.2564443                     2009                        04
## 296           0.4953493                     2009                        08
## 295           0.1496424                     2009                        08
## 279           0.2403083                     2009                        05
## 213           0.8583831                     2008                        01
##     Week.bgn.date.fctr.nonNA ILI.log.predict.Final.lm
## 278                  (25,31]                0.1623063
## 296                  (25,31]                0.5758465
## 295                  (19,25]                0.2119733
## 279                 (0.97,7]                0.2955065
## 213                  (25,31]                0.9013407
##     ILI.log.predict.Final.lm.err                  .label
## 278                    0.9301501 2009-04-26 - 2009-05-02
## 296                    0.7377949 2009-08-30 - 2009-09-05
## 295                    0.6929166 2009-08-23 - 2009-08-29
## 279                    0.5952705 2009-05-03 - 2009-05-09
## 213                    0.5879186 2008-01-27 - 2008-02-02
```

![](Google_Flu_template2_files/figure-html/fit.data.training_1-7.png) 

```r
dsp_feats_vctr <- c(NULL)
for(var in grep(".importance", names(glb_feats_df), fixed=TRUE, value=TRUE))
    dsp_feats_vctr <- union(dsp_feats_vctr, 
                            glb_feats_df[!is.na(glb_feats_df[, var]), "id"])

# print(glb_trnobs_df[glb_trnobs_df$UniqueID %in% FN_OOB_ids, 
#                     grep(glb_rsp_var, names(glb_trnobs_df), value=TRUE)])

print(setdiff(names(glb_trnobs_df), names(glb_allobs_df)))
```

```
## [1] "ILI.log.predict.Final.lm"     "ILI.log.predict.Final.lm.err"
```

```r
for (col in setdiff(names(glb_trnobs_df), names(glb_allobs_df)))
    # Merge or cbind ?
    glb_allobs_df[glb_allobs_df$.src == "Train", col] <- glb_trnobs_df[, col]

print(setdiff(names(glb_fitobs_df), names(glb_allobs_df)))
```

```
## character(0)
```

```r
print(setdiff(names(glb_OOBobs_df), names(glb_allobs_df)))
```

```
## character(0)
```

```r
for (col in setdiff(names(glb_OOBobs_df), names(glb_allobs_df)))
    # Merge or cbind ?
    glb_allobs_df[glb_allobs_df$.lcn == "OOB", col] <- glb_OOBobs_df[, col]
    
print(setdiff(names(glb_newobs_df), names(glb_allobs_df)))
```

```
## character(0)
```

```r
if (glb_save_envir)
    save(glb_feats_df, glb_allobs_df, 
         #glb_trnobs_df, glb_fitobs_df, glb_OOBobs_df, glb_newobs_df,
         glb_models_df, dsp_models_df, glb_models_lst, glb_model_type,
         glb_sel_mdl, glb_sel_mdl_id,
         glb_fin_mdl, glb_fin_mdl_id,
        file=paste0(glb_out_pfx, "dsk.RData"))

replay.petrisim(pn=glb_analytics_pn, 
    replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
        "data.training.all.prediction","model.final")), flip_coord=TRUE)
```

```
## time	trans	 "bgn " "fit.data.training.all " "predict.data.new " "end " 
## 0.0000 	multiple enabled transitions:  data.training.all data.new model.selected 	firing:  data.training.all 
## 1.0000 	 1 	 2 1 0 0 
## 1.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction 	firing:  data.new 
## 2.0000 	 2 	 1 1 1 0 
## 2.0000 	multiple enabled transitions:  data.training.all data.new model.selected model.final data.training.all.prediction data.new.prediction 	firing:  model.selected 
## 3.0000 	 3 	 0 2 1 0 
## 3.0000 	multiple enabled transitions:  model.final data.training.all.prediction data.new.prediction 	firing:  data.training.all.prediction 
## 4.0000 	 5 	 0 1 1 1 
## 4.0000 	multiple enabled transitions:  model.final data.training.all.prediction data.new.prediction 	firing:  model.final 
## 5.0000 	 4 	 0 0 2 1
```

![](Google_Flu_template2_files/figure-html/fit.data.training_1-8.png) 

```r
glb_chunks_df <- myadd_chunk(glb_chunks_df, "predict.data.new", major.inc=TRUE)
```

```
##                label step_major step_minor     bgn     end elapsed
## 15 fit.data.training          8          1 188.577 193.019   4.442
## 16  predict.data.new          9          0 193.020      NA      NA
```

## Step `9.0: predict data new`

```r
# Compute final model predictions
# sav_newobs_df <- glb_newobs_df
glb_newobs_df <- glb_get_predictions(glb_newobs_df, mdl_id=glb_fin_mdl_id, 
                                     rsp_var_out=glb_rsp_var_out,
    prob_threshold_def=ifelse(glb_is_classification && glb_is_binomial, 
        glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, 
                      "opt.prob.threshold.OOB"], NULL))
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

![](Google_Flu_template2_files/figure-html/predict.data.new-1.png) 

```
##                        Week       ILI   Queries .src     ILI.log
## 419 2012-01-08 - 2012-01-14 1.5434005 0.4993360 Test  0.43398812
## 430 2012-03-25 - 2012-03-31 1.7423860 0.3652058 Test  0.55525545
## 464 2012-11-18 - 2012-11-24 2.3046254 0.5112882 Test  0.83491815
## 445 2012-07-08 - 2012-07-14 0.9281519 0.2656042 Test -0.07455986
## 420 2012-01-15 - 2012-01-21 1.6476154 0.5006640 Test  0.49932902
## 441 2012-06-10 - 2012-06-16 1.0861211 0.2509960 Test  0.08261272
##         .rnorm   Week.bgn   Week.end ILI.2.lag ILI.2.lag.log
## 419 -1.3166927 2012-01-08 2012-01-14  2.124130    0.75336227
## 430 -0.2262124 2012-03-25 2012-03-31  2.293422    0.83004483
## 464 -0.4158885 2012-11-18 2012-11-24  1.610915    0.47680233
## 445  0.6476146 2012-07-08 2012-07-14  1.078713    0.07576853
## 420  0.8742224 2012-01-15 2012-01-21  1.766707    0.56911730
## 441 -1.6211072 2012-06-10 2012-06-16  1.299020    0.26160987
##     Week.bgn.POSIX Week.bgn.year.fctr Week.bgn.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
## 441     2012-06-10               2012                  06
##     Week.bgn.date.fctr Week.bgn.wkday.fctr Week.bgn.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
## 441             (7,13]                   0              1
##     Week.bgn.hour.fctr Week.bgn.minute.fctr Week.bgn.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
## 441                  0                    0                    0
##     Week.end.POSIX Week.end.year.fctr Week.end.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
## 441     2012-06-10               2012                  06
##     Week.end.date.fctr Week.end.wkday.fctr Week.end.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
## 441             (7,13]                   0              1
##     Week.end.hour.fctr Week.end.minute.fctr Week.end.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
## 441                  0                    0                    0
##     Week.bgn.zoo Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
## 419   1325394000           13.31265           14.00580            15.61583
## 430   1325998800           13.31265           14.00282            15.61464
## 464   1326603600           13.31265           14.00877            15.61583
## 445   1327208400           13.31265           14.00580            15.61524
## 420   1327813200           13.31265           14.00580            15.61583
## 441   1328418000           13.31265           14.00580            15.61524
##     Week.bgn.last100.log Week.end.zoo Week.end.last1.log
## 419             17.91782   1325394000           13.31265
## 430             17.91782   1325998800           13.31265
## 464             17.91782   1326603600           13.31265
## 445             17.91782   1327208400           13.31265
## 420             17.91782   1327813200           13.31265
## 441             17.91782   1328418000           13.31265
##     Week.end.last2.log Week.end.last10.log Week.end.last100.log
## 419           14.00580            15.61583             17.91782
## 430           14.00282            15.61464             17.91782
## 464           14.00877            15.61583             17.91782
## 445           14.00580            15.61524             17.91782
## 420           14.00580            15.61583             17.91782
## 441           14.00580            15.61524             17.91782
##     ILI.2.lag.log.nonNA Week.bgn.year.fctr.nonNA Week.bgn.month.fctr.nonNA
## 419          0.75336227                     2012                        01
## 430          0.83004483                     2012                        03
## 464          0.47680233                     2012                        11
## 445          0.07576853                     2012                        07
## 420          0.56911730                     2012                        01
## 441          0.26160987                     2012                        06
##     Week.bgn.date.fctr.nonNA ILI.log.predict.Final.lm
## 419                   (7,13]                0.9316064
## 430                  (19,25]                0.8210195
## 464                  (13,19]                0.5772443
## 445                   (7,13]                0.1496508
## 420                  (13,19]                0.7170902
## 441                   (7,13]                0.2986348
##     ILI.log.predict.Final.lm.err
## 419                    0.4976183
## 430                    0.2657641
## 464                    0.2576738
## 445                    0.2242107
## 420                    0.2177612
## 441                    0.2160221
```

```r
if (glb_is_classification && glb_is_binomial)
    glb_analytics_diag_plots(obs_df=glb_newobs_df, mdl_id=glb_fin_mdl_id, 
            prob_threshold=glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, 
                                         "opt.prob.threshold.OOB"]) else
    glb_analytics_diag_plots(obs_df=glb_newobs_df, mdl_id=glb_fin_mdl_id)                  
```

```
## Warning in glb_analytics_diag_plots(obs_df = glb_newobs_df, mdl_id =
## glb_fin_mdl_id): Limiting important feature scatter plots to 5 out of 7
```

![](Google_Flu_template2_files/figure-html/predict.data.new-2.png) ![](Google_Flu_template2_files/figure-html/predict.data.new-3.png) ![](Google_Flu_template2_files/figure-html/predict.data.new-4.png) ![](Google_Flu_template2_files/figure-html/predict.data.new-5.png) ![](Google_Flu_template2_files/figure-html/predict.data.new-6.png) 

```
##                        Week       ILI   Queries .src     ILI.log
## 419 2012-01-08 - 2012-01-14 1.5434005 0.4993360 Test  0.43398812
## 430 2012-03-25 - 2012-03-31 1.7423860 0.3652058 Test  0.55525545
## 464 2012-11-18 - 2012-11-24 2.3046254 0.5112882 Test  0.83491815
## 445 2012-07-08 - 2012-07-14 0.9281519 0.2656042 Test -0.07455986
## 420 2012-01-15 - 2012-01-21 1.6476154 0.5006640 Test  0.49932902
##         .rnorm   Week.bgn   Week.end ILI.2.lag ILI.2.lag.log
## 419 -1.3166927 2012-01-08 2012-01-14  2.124130    0.75336227
## 430 -0.2262124 2012-03-25 2012-03-31  2.293422    0.83004483
## 464 -0.4158885 2012-11-18 2012-11-24  1.610915    0.47680233
## 445  0.6476146 2012-07-08 2012-07-14  1.078713    0.07576853
## 420  0.8742224 2012-01-15 2012-01-21  1.766707    0.56911730
##     Week.bgn.POSIX Week.bgn.year.fctr Week.bgn.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
##     Week.bgn.date.fctr Week.bgn.wkday.fctr Week.bgn.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
##     Week.bgn.hour.fctr Week.bgn.minute.fctr Week.bgn.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
##     Week.end.POSIX Week.end.year.fctr Week.end.month.fctr
## 419     2012-01-08               2012                  01
## 430     2012-03-25               2012                  03
## 464     2012-11-18               2012                  11
## 445     2012-07-08               2012                  07
## 420     2012-01-15               2012                  01
##     Week.end.date.fctr Week.end.wkday.fctr Week.end.wkend
## 419             (7,13]                   0              1
## 430            (19,25]                   0              1
## 464            (13,19]                   0              1
## 445             (7,13]                   0              1
## 420            (13,19]                   0              1
##     Week.end.hour.fctr Week.end.minute.fctr Week.end.second.fctr
## 419                  0                    0                    0
## 430                  0                    0                    0
## 464                  0                    0                    0
## 445                  0                    0                    0
## 420                  0                    0                    0
##     Week.bgn.zoo Week.bgn.last1.log Week.bgn.last2.log Week.bgn.last10.log
## 419   1325394000           13.31265           14.00580            15.61583
## 430   1325998800           13.31265           14.00282            15.61464
## 464   1326603600           13.31265           14.00877            15.61583
## 445   1327208400           13.31265           14.00580            15.61524
## 420   1327813200           13.31265           14.00580            15.61583
##     Week.bgn.last100.log Week.end.zoo Week.end.last1.log
## 419             17.91782   1325394000           13.31265
## 430             17.91782   1325998800           13.31265
## 464             17.91782   1326603600           13.31265
## 445             17.91782   1327208400           13.31265
## 420             17.91782   1327813200           13.31265
##     Week.end.last2.log Week.end.last10.log Week.end.last100.log
## 419           14.00580            15.61583             17.91782
## 430           14.00282            15.61464             17.91782
## 464           14.00877            15.61583             17.91782
## 445           14.00580            15.61524             17.91782
## 420           14.00580            15.61583             17.91782
##     ILI.2.lag.log.nonNA Week.bgn.year.fctr.nonNA Week.bgn.month.fctr.nonNA
## 419          0.75336227                     2012                        01
## 430          0.83004483                     2012                        03
## 464          0.47680233                     2012                        11
## 445          0.07576853                     2012                        07
## 420          0.56911730                     2012                        01
##     Week.bgn.date.fctr.nonNA ILI.log.predict.Final.lm
## 419                   (7,13]                0.9316064
## 430                  (19,25]                0.8210195
## 464                  (13,19]                0.5772443
## 445                   (7,13]                0.1496508
## 420                  (13,19]                0.7170902
##     ILI.log.predict.Final.lm.err                  .label
## 419                    0.4976183 2012-01-08 - 2012-01-14
## 430                    0.2657641 2012-03-25 - 2012-03-31
## 464                    0.2576738 2012-11-18 - 2012-11-24
## 445                    0.2242107 2012-07-08 - 2012-07-14
## 420                    0.2177612 2012-01-15 - 2012-01-21
```

![](Google_Flu_template2_files/figure-html/predict.data.new-7.png) 

```r
if (glb_is_classification && glb_is_binomial) {
    submit_df <- glb_newobs_df[, c(glb_id_var, 
                                   paste0(glb_rsp_var_out, glb_fin_mdl_id, ".prob"))]
    names(submit_df)[2] <- "Probability1"
#     submit_df <- glb_newobs_df[, c(paste0(glb_rsp_var_out, glb_fin_mdl_id)), FALSE]
#     names(submit_df)[1] <- "BDscience"
#     submit_df$BDscience <- as.numeric(submit_df$BDscience) - 1
#     #submit_df <-rbind(submit_df, data.frame(bdanalytics=c(" ")))
#     print("Submission Stats:")
#     print(table(submit_df$BDscience, useNA = "ifany"))
} else submit_df <- glb_newobs_df[, c(glb_id_var, 
                                   paste0(glb_rsp_var_out, glb_fin_mdl_id))]
submit_fname <- paste0(gsub(".", "_", paste0(glb_out_pfx, glb_fin_mdl_id), fixed=TRUE), 
                    "_submit.csv")
write.csv(submit_df, submit_fname, quote=FALSE, row.names=FALSE)
#cat(" ", "\n", file=submit_fn, append=TRUE)

# print(orderBy(~ -max.auc.OOB, glb_models_df[, c("model_id", 
#             "max.auc.OOB", "max.Accuracy.OOB")]))
if (glb_is_classification && glb_is_binomial)
    print(glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, 
                        "opt.prob.threshold.OOB"])
print(sprintf("glb_sel_mdl_id: %s", glb_sel_mdl_id))
```

```
## [1] "glb_sel_mdl_id: Interact.High.cor.Y.lm"
```

```r
print(sprintf("glb_fin_mdl_id: %s", glb_fin_mdl_id))
```

```
## [1] "glb_fin_mdl_id: Final.lm"
```

```r
print(dim(glb_fitobs_df))
```

```
## [1] 417  42
```

```r
print(dsp_models_df)
```

```
##                     model_id min.RMSE.OOB max.R.sq.OOB max.Adj.R.sq.fit
## 6     Interact.High.cor.Y.lm    0.1361616   0.88690870     0.8818277999
## 7               Low.cor.X.lm    0.1402968   0.87993525     0.8979018265
## 13             Flu.Trend2.lm    0.1533303   0.85659116     0.9028365494
## 5               Max.cor.Y.lm    0.1806841   0.80085932     0.8466884911
## 12         All.X.no.rnorm.rf    0.1885728   0.78338031               NA
## 3  Max.cor.Y.cv.0.cp.0.rpart    0.2047780   0.74420816               NA
## 4            Max.cor.Y.rpart    0.2650357   0.57152193               NA
## 11      All.X.no.rnorm.rpart    0.2650357   0.57152193               NA
## 10            All.X.bayesglm    0.3306004   0.33330588               NA
## 8                   All.X.lm    0.3320580   0.32741401     0.9266370734
## 9                  All.X.glm    0.3320580   0.32741401               NA
## 2       Max.cor.Y.cv.0.rpart    0.4048928   0.00000000               NA
## 1                     MFO.lm    0.4080698  -0.01575461    -0.0008029752
```

```r
if (glb_is_regression) {
    print(sprintf("%s OOB RMSE: %0.4f", glb_sel_mdl_id,
                  glb_models_df[glb_models_df$model_id == glb_sel_mdl_id, "min.RMSE.OOB"]))

    if (!is.null(glb_category_vars)) {
        stop("not implemented yet")
        tmp_OOBobs_df <- glb_OOBobs_df[, c(glb_category_vars, predct_accurate_var_name)]
        names(tmp_OOBobs_df)[length(names(tmp_OOBobs_df))] <- "accurate.OOB"
        aOOB_ctgry_df <- mycreate_xtab_df(tmp_OOBobs_df, names(tmp_OOBobs_df)) 
        aOOB_ctgry_df[is.na(aOOB_ctgry_df)] <- 0
        aOOB_ctgry_df <- mutate(aOOB_ctgry_df, 
                                .n.OOB = accurate.OOB.FALSE + accurate.OOB.TRUE,
                                max.accuracy.OOB = accurate.OOB.TRUE / .n.OOB)
        #intersect(names(glb_ctgry_df), names(aOOB_ctgry_df))
        glb_ctgry_df <- merge(glb_ctgry_df, aOOB_ctgry_df, all=TRUE)
        print(orderBy(~-accurate.OOB.FALSE, glb_ctgry_df))
    }
    
    if ((glb_rsp_var %in% names(glb_newobs_df)) &&
        !(any(is.na(glb_newobs_df[, glb_rsp_var])))) {
            pred_stats_df <- 
                mypredict_mdl(mdl=glb_models_lst[[glb_fin_mdl_id]], 
                              df=glb_newobs_df, 
                              rsp_var=glb_rsp_var, 
                              rsp_var_out=glb_rsp_var_out, 
                              model_id_method=glb_fin_mdl_id, 
                              label="new",
						      model_summaryFunction=glb_sel_mdl$control$summaryFunction, 
						      model_metric=glb_sel_mdl$metric,
						      model_metric_maximize=glb_sel_mdl$maximize,
						      ret_type="stats")        
            print(sprintf("%s prediction stats for glb_newobs_df:", glb_fin_mdl_id))
            print(pred_stats_df)
    }    
}    
```

```
## [1] "Interact.High.cor.Y.lm OOB RMSE: 0.1362"
```

```
## Warning: contrasts dropped from factor Week.bgn.year.fctr.nonNA
```

```
## Warning: contrasts dropped from factor Week.bgn.month.fctr.nonNA
```

```
## Warning in predict.lm(modelFit, newdata): prediction from a rank-deficient
## fit may be misleading
```

```
## [1] "Final.lm prediction stats for glb_newobs_df:"
##   model_id max.R.sq.new min.RMSE.new
## 1 Final.lm    0.8869087    0.1361616
```

```r
if (glb_is_classification) {
    print(sprintf("%s OOB confusion matrix & accuracy: ", glb_sel_mdl_id))
    print(t(confusionMatrix(glb_OOBobs_df[, paste0(glb_rsp_var_out, glb_sel_mdl_id)], 
                            glb_OOBobs_df[, glb_rsp_var])$table))

    if (!is.null(glb_category_vars)) {
        tmp_OOBobs_df <- glb_OOBobs_df[, c(glb_category_vars, predct_accurate_var_name)]
        names(tmp_OOBobs_df)[length(names(tmp_OOBobs_df))] <- "accurate.OOB"
        aOOB_ctgry_df <- mycreate_xtab_df(tmp_OOBobs_df, names(tmp_OOBobs_df)) 
        aOOB_ctgry_df[is.na(aOOB_ctgry_df)] <- 0
        aOOB_ctgry_df <- mutate(aOOB_ctgry_df, 
                                .n.OOB = accurate.OOB.FALSE + accurate.OOB.TRUE,
                                max.accuracy.OOB = accurate.OOB.TRUE / .n.OOB)
        #intersect(names(glb_ctgry_df), names(aOOB_ctgry_df))
        glb_ctgry_df <- merge(glb_ctgry_df, aOOB_ctgry_df, all=TRUE)
        print(orderBy(~-accurate.OOB.FALSE, glb_ctgry_df))
    }
    
    if ((glb_rsp_var %in% names(glb_newobs_df)) &&
        !(any(is.na(glb_newobs_df[, glb_rsp_var])))) {
        print(sprintf("%s new confusion matrix & accuracy: ", glb_fin_mdl_id))
        print(t(confusionMatrix(glb_newobs_df[, paste0(glb_rsp_var_out, glb_fin_mdl_id)], 
                                glb_newobs_df[, glb_rsp_var])$table))
    }    

}    

dsp_myCategory_conf_mtrx <- function(myCategory) {
    print(sprintf("%s OOB::myCategory=%s confusion matrix & accuracy: ", 
                  glb_sel_mdl_id, myCategory))
    print(t(confusionMatrix(
        glb_OOBobs_df[glb_OOBobs_df$myCategory == myCategory, 
                      paste0(glb_rsp_var_out, glb_sel_mdl_id)], 
        glb_OOBobs_df[glb_OOBobs_df$myCategory == myCategory, glb_rsp_var])$table))
    print(sum(glb_OOBobs_df[glb_OOBobs_df$myCategory == myCategory, 
                            predct_accurate_var_name]) / 
         nrow(glb_OOBobs_df[glb_OOBobs_df$myCategory == myCategory, ]))
    err_ids <- glb_OOBobs_df[(glb_OOBobs_df$myCategory == myCategory) & 
                             (!glb_OOBobs_df[, predct_accurate_var_name]), glb_id_var]

    OOB_FNerr_df <- glb_OOBobs_df[(glb_OOBobs_df$UniqueID %in% err_ids) & 
                               (glb_OOBobs_df$Popular == 1), 
                        c(
                            ".clusterid", 
                            "Popular", "Headline", "Snippet", "Abstract")]
    print(sprintf("%s OOB::myCategory=%s FN errors: %d", glb_sel_mdl_id, myCategory,
                  nrow(OOB_FNerr_df)))
    print(OOB_FNerr_df)

    OOB_FPerr_df <- glb_OOBobs_df[(glb_OOBobs_df$UniqueID %in% err_ids) & 
                               (glb_OOBobs_df$Popular == 0), 
                        c(
                            ".clusterid", 
                            "Popular", "Headline", "Snippet", "Abstract")]
    print(sprintf("%s OOB::myCategory=%s FP errors: %d", glb_sel_mdl_id, myCategory,
                  nrow(OOB_FPerr_df)))
    print(OOB_FPerr_df)
}
#dsp_myCategory_conf_mtrx(myCategory="OpEd#Opinion#")
#dsp_myCategory_conf_mtrx(myCategory="Business#Business Day#Dealbook")
#dsp_myCategory_conf_mtrx(myCategory="##")

# if (glb_is_classification) {
#     print("FN_OOB_ids:")
#     print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                         grep(glb_rsp_var, names(glb_OOBobs_df), value=TRUE)])
#     print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                         glb_txt_vars])
#     print(dsp_vctr <- colSums(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                         setdiff(grep("[HSA].", names(glb_OOBobs_df), value=TRUE),
#                                 union(myfind_chr_cols_df(glb_OOBobs_df),
#                     grep(".fctr", names(glb_OOBobs_df), fixed=TRUE, value=TRUE)))]))
# }

dsp_hdlpfx_results <- function(hdlpfx) {
    print(hdlpfx)
    print(glb_OOBobs_df[glb_OOBobs_df$Headline.pfx %in% c(hdlpfx), 
                        grep(glb_rsp_var, names(glb_OOBobs_df), value=TRUE)])
    print(glb_newobs_df[glb_newobs_df$Headline.pfx %in% c(hdlpfx), 
                        grep(glb_rsp_var, names(glb_newobs_df), value=TRUE)])
    print(dsp_vctr <- colSums(glb_newobs_df[glb_newobs_df$Headline.pfx %in% c(hdlpfx), 
                        setdiff(grep("[HSA]\\.", names(glb_newobs_df), value=TRUE),
                                union(myfind_chr_cols_df(glb_newobs_df),
                    grep(".fctr", names(glb_newobs_df), fixed=TRUE, value=TRUE)))]))
    print(dsp_vctr <- dsp_vctr[dsp_vctr != 0])
    print(glb_newobs_df[glb_newobs_df$Headline.pfx %in% c(hdlpfx), 
                        union(names(dsp_vctr), myfind_chr_cols_df(glb_newobs_df))])
}
#dsp_hdlpfx_results(hdlpfx="Ask Well::")

# print("myMisc::|OpEd|blank|blank|1:")
# print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% c(6446), 
#                     grep(glb_rsp_var, names(glb_OOBobs_df), value=TRUE)])

# print(glb_OOBobs_df[glb_OOBobs_df$UniqueID %in% FN_OOB_ids, 
#                     c("WordCount", "WordCount.log", "myMultimedia",
#                       "NewsDesk", "SectionName", "SubsectionName")])
# print(mycreate_sqlxtab_df(glb_allobs_df[sel_obs(Headline.contains="[Vv]ideo"), ], 
#                           c(glb_rsp_var, "myMultimedia")))
# dsp_chisq.test(Headline.contains="[Vi]deo")
# print(glb_allobs_df[sel_obs(Headline.contains="[Vv]ideo"), 
#                           c(glb_rsp_var, "Popular", "myMultimedia", "Headline")])
# print(glb_allobs_df[sel_obs(Headline.contains="[Ee]bola", Popular=1), 
#                           c(glb_rsp_var, "Popular", "myMultimedia", "Headline",
#                             "NewsDesk", "SectionName", "SubsectionName")])
# print(subset(glb_feats_df, !is.na(importance))[,
#     c("is.ConditionalX.y", 
#       grep("importance", names(glb_feats_df), fixed=TRUE, value=TRUE))])
# print(subset(glb_feats_df, is.ConditionalX.y & is.na(importance))[,
#     c("is.ConditionalX.y", 
#       grep("importance", names(glb_feats_df), fixed=TRUE, value=TRUE))])
# print(subset(glb_feats_df, !is.na(importance))[,
#     c("zeroVar", "nzv", "myNearZV", 
#       grep("importance", names(glb_feats_df), fixed=TRUE, value=TRUE))])
# print(subset(glb_feats_df, is.na(importance))[,
#     c("zeroVar", "nzv", "myNearZV", 
#       grep("importance", names(glb_feats_df), fixed=TRUE, value=TRUE))])
print(orderBy(as.formula(paste0("~ -", glb_sel_mdl_id, ".importance")), glb_featsimp_df))
```

```
##                                                   Interact.High.cor.Y.lm.importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`                         100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`                          79.081237
## ILI.2.lag.log.nonNA                                                       75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`                         41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`                         39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`                          38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`                                  24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`                                24.897758
## Week.bgn.last100.log                                                      22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`                           22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`                          19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`                           17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`                           15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`                           15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`                           14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`                          14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`                         13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`                           13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`                          13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`                           10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`                          10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`                           5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`                                   0.000000
##                                                   importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`  100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`   79.081237
## ILI.2.lag.log.nonNA                                75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`  41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`  39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`   38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`           24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`         24.897758
## Week.bgn.last100.log                               22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`    22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`   19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`    17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`    15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`    15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`    14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`   14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`  13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`    13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`   13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`    10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`   10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`    5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`            0.000000
##                                                   Final.lm.importance
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA4`           100.000000
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA3`            79.081237
## ILI.2.lag.log.nonNA                                         75.169274
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA10`           41.185653
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA12`           39.486108
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA6`            38.096300
## `ILI.2.lag.log.nonNA:Week.bgn.last2.log`                    24.993131
## `ILI.2.lag.log.nonNA:Week.bgn.last100.log`                  24.897758
## Week.bgn.last100.log                                        22.448374
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA2`             22.305442
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA8`            19.684026
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA6`             17.556280
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA7`             15.831464
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA8`             15.178728
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA4`             14.912604
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA5`            14.258897
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA11`           13.972025
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA3`             13.438086
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA9`            13.372715
## `ILI.2.lag.log.nonNA:Week.bgn.year.fctr.nonNA5`             10.658250
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA7`            10.032455
## `ILI.2.lag.log.nonNA:Week.bgn.month.fctr.nonNA2`             5.653516
## `ILI.2.lag.log.nonNA:Week.bgn.last1.log`                     0.000000
```

```r
# players_df <- data.frame(id=c("Chavez", "Giambi", "Menechino", "Myers", "Pena"),
#                          OBP=c(0.338, 0.391, 0.369, 0.313, 0.361),
#                          SLG=c(0.540, 0.450, 0.374, 0.447, 0.500),
#                         cost=c(1400000, 1065000, 295000, 800000, 300000))
# players_df$RS.predict <- predict(glb_models_lst[[csm_mdl_id]], players_df)
# print(orderBy(~ -RS.predict, players_df))

if (length(diff <- setdiff(names(glb_trnobs_df), names(glb_allobs_df))) > 0)   
    print(diff)
for (col in setdiff(names(glb_trnobs_df), names(glb_allobs_df)))
    # Merge or cbind ?
    glb_allobs_df[glb_allobs_df$.src == "Train", col] <- glb_trnobs_df[, col]

if (length(diff <- setdiff(names(glb_fitobs_df), names(glb_allobs_df))) > 0)   
    print(diff)
if (length(diff <- setdiff(names(glb_OOBobs_df), names(glb_allobs_df))) > 0)   
    print(diff)

for (col in setdiff(names(glb_OOBobs_df), names(glb_allobs_df)))
    # Merge or cbind ?
    glb_allobs_df[glb_allobs_df$.lcn == "OOB", col] <- glb_OOBobs_df[, col]
    
if (length(diff <- setdiff(names(glb_newobs_df), names(glb_allobs_df))) > 0)   
    print(diff)

if (glb_save_envir)
    save(glb_feats_df, glb_allobs_df, 
         #glb_trnobs_df, glb_fitobs_df, glb_OOBobs_df, glb_newobs_df,
         glb_models_df, dsp_models_df, glb_models_lst, glb_model_type,
         glb_sel_mdl, glb_sel_mdl_id,
         glb_fin_mdl, glb_fin_mdl_id,
        file=paste0(glb_out_pfx, "prdnew_dsk.RData"))

rm(submit_df, tmp_OOBobs_df)
```

```
## Warning in rm(submit_df, tmp_OOBobs_df): object 'tmp_OOBobs_df' not found
```

```r
# tmp_replay_lst <- replay.petrisim(pn=glb_analytics_pn, 
#     replay.trans=(glb_analytics_avl_objs <- c(glb_analytics_avl_objs, 
#         "data.new.prediction")), flip_coord=TRUE)
# print(ggplot.petrinet(tmp_replay_lst[["pn"]]) + coord_flip())

glb_chunks_df <- myadd_chunk(glb_chunks_df, "display.session.info", major.inc=TRUE)
```

```
##                   label step_major step_minor     bgn     end elapsed
## 16     predict.data.new          9          0 193.020 197.287   4.267
## 17 display.session.info         10          0 197.288      NA      NA
```

Null Hypothesis ($\sf{H_{0}}$): mpg is not impacted by am_fctr.  
The variance by am_fctr appears to be independent. 
#```{r q1, cache=FALSE}
# print(t.test(subset(cars_df, am_fctr == "automatic")$mpg, 
#              subset(cars_df, am_fctr == "manual")$mpg, 
#              var.equal=FALSE)$conf)
#```
We reject the null hypothesis i.e. we have evidence to conclude that am_fctr impacts mpg (95% confidence). Manual transmission is better for miles per gallon versus automatic transmission.


```
##                      label step_major step_minor     bgn     end elapsed
## 5         extract.features          3          0  10.224 109.764  99.540
## 11              fit.models          7          1 145.145 168.972  23.828
## 6             cluster.data          4          0 109.765 127.266  17.501
## 10              fit.models          7          0 129.347 145.145  15.798
## 12              fit.models          7          2 168.973 180.142  11.169
## 14       fit.data.training          8          0 184.108 188.576   4.468
## 15       fit.data.training          8          1 188.577 193.019   4.442
## 16        predict.data.new          9          0 193.020 197.287   4.267
## 13              fit.models          7          3 180.142 184.107   3.965
## 2             inspect.data          2          0   7.173   9.516   2.343
## 7      manage.missing.data          4          1 127.267 128.257   0.990
## 8          select.features          5          0 128.258 128.999   0.742
## 3               scrub.data          2          1   9.516  10.161   0.645
## 9  partition.data.training          6          0 129.000 129.347   0.347
## 1              import.data          1          0   6.863   7.173   0.310
## 4           transform.data          2          2  10.161  10.224   0.063
##    duration
## 5    99.540
## 11   23.827
## 6    17.501
## 10   15.798
## 12   11.169
## 14    4.468
## 15    4.442
## 16    4.267
## 13    3.965
## 2     2.343
## 7     0.990
## 8     0.741
## 3     0.645
## 9     0.347
## 1     0.310
## 4     0.063
## [1] "Total Elapsed Time: 197.287 secs"
```

![](Google_Flu_template2_files/figure-html/display.session.info-1.png) 

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.10.3 (Yosemite)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] gdata_2.16.1        randomForest_4.6-10 arm_1.8-5          
##  [4] lme4_1.1-7          Matrix_1.2-1        MASS_7.3-40        
##  [7] rpart.plot_1.5.2    rpart_4.1-9         reshape2_1.4.1     
## [10] mice_2.22           Rcpp_0.11.6         XML_3.98-1.2       
## [13] dplyr_0.4.1         plyr_1.8.2          zoo_1.7-12         
## [16] doMC_1.3.3          iterators_1.0.7     foreach_1.4.2      
## [19] doBy_4.5-13         survival_2.38-1     caret_6.0-47       
## [22] ggplot2_1.0.1       lattice_0.20-31    
## 
## loaded via a namespace (and not attached):
##  [1] gtools_3.5.0        splines_3.2.0       colorspace_1.2-6   
##  [4] htmltools_0.2.6     yaml_2.1.13         mgcv_1.8-6         
##  [7] nloptr_1.0.4        DBI_0.3.1           RColorBrewer_1.1-2 
## [10] stringr_1.0.0       munsell_0.4.2       gtable_0.1.2       
## [13] codetools_0.2-11    coda_0.17-1         evaluate_0.7       
## [16] labeling_0.3        knitr_1.10.5        SparseM_1.6        
## [19] quantreg_5.11       pbkrtest_0.4-2      proto_0.3-10       
## [22] scales_0.2.4        formatR_1.2         BradleyTerry2_1.0-6
## [25] abind_1.4-3         digest_0.6.8        stringi_0.4-1      
## [28] brglm_0.5-9         tools_3.2.0         magrittr_1.5       
## [31] lazyeval_0.1.10     car_2.0-25          assertthat_0.1     
## [34] minqa_1.2.4         rmarkdown_0.6.1     nnet_7.3-9         
## [37] nlme_3.1-120        compiler_3.2.0
```
