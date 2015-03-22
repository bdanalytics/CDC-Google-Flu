# CDC+Google: ILI.log regression
bdanalytics  

**  **    
**Date: (Sun) Mar 22, 2015**    

# Introduction:  

Data: 
Source: 
    Training:   https://courses.edx.org/c4x/MITx/15.071x_2/asset/FluTrain.csv  
    New:        https://courses.edx.org/c4x/MITx/15.071x_2/asset/FluTest.csv  
Time period: 



# Synopsis:

Based on analysis utilizing <> techniques, <conclusion heading>:  

### ![](<filename>.png)

## Potential next steps include:

# Analysis: 

```r
rm(list=ls())
set.seed(12345)
options(stringsAsFactors=FALSE)
source("~/Dropbox/datascience/R/mydsutils.R")
source("~/Dropbox/datascience/R/myplot.R")
source("~/Dropbox/datascience/R/mypetrinet.R")
# Gather all package requirements here
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(zoo))

#require(sos); findFn("pinv", maxPages=2, sortby="MaxScore")

# Analysis control global variables
glb_is_separate_predict_dataset <- TRUE    # or TRUE
glb_predct_var <- "ILI.log"           # or NULL
glb_predct_var_name <- paste0(glb_predct_var, ".predict")
glb_id_vars <- c("Week")                # or NULL

glb_exclude_vars_as_features <- glb_id_vars     # or NULL
# List chrs converted into factors; num/int transformed 
glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
                                      c("ILI", "ILI.log",
                                        "ILI.Lag.2")     # or NULL
                                      )
# List feats that shd be excluded due to known causation by prediction variable
# glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
#                                       c("<col_name>")     # or NULL
#                                       )

glb_is_regression <- TRUE; glb_is_classification <- FALSE

glb_mdl <- glb_sel_mdl <- NULL
glb_models_df <- data.frame()

script_df <- data.frame(chunk_label="import_data", chunk_step_major=1, chunk_step_minor=0)
print(script_df)
```

```
##   chunk_label chunk_step_major chunk_step_minor
## 1 import_data                1                0
```

## Step `1`: import data

```r
glb_entity_df <- myimport_data(
    url="https://courses.edx.org/c4x/MITx/15.071x_2/asset/FluTrain.csv", 
    comment="glb_entity_df", force_header=TRUE, print_diagn=TRUE)
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
##                        Week      ILI   Queries
## 69  2005-04-24 - 2005-04-30 1.207025 0.1221780
## 189 2007-08-12 - 2007-08-18 0.618075 0.1221780
## 301 2009-10-04 - 2009-10-10 5.660867 0.7436919
## 316 2010-01-17 - 2010-01-23 1.926056 0.3784861
## 365 2010-12-26 - 2011-01-01 3.431723 0.8061089
## 367 2011-01-09 - 2011-01-15 2.910629 0.5670651
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
##  - attr(*, "comment")= chr "glb_entity_df"
## NULL
```

```r
if (glb_is_separate_predict_dataset) {
    glb_predct_df <- myimport_data(
        url="https://courses.edx.org/c4x/MITx/15.071x_2/asset/FluTest.csv", 
        comment="glb_predct_df", force_header=TRUE, print_diagn=TRUE)
} else {
#     glb_predct_df <- subset(glb_entity_df, <condition>)
    glb_predct_df <- glb_entity_df[sample(1:nrow(glb_entity_df), 
                                          max(2, nrow(glb_entity_df) / 1000)),]
    comment(glb_predct_df) <- "glb_predct_df"
    myprint_df(glb_predct_df)
    str(glb_predct_df)
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
## 2  2012-01-08 - 2012-01-14 1.543401 0.4993360
## 8  2012-02-19 - 2012-02-25 2.103851 0.5006640
## 17 2012-04-22 - 2012-04-28 1.288493 0.3120850
## 26 2012-06-24 - 2012-06-30 1.078713 0.2390438
## 37 2012-09-09 - 2012-09-15 1.186489 0.3585657
## 49 2012-12-02 - 2012-12-08 2.978047 0.6719788
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
##  - attr(*, "comment")= chr "glb_predct_df"
## NULL
```

```r
# glb_entity_df <- subset(glb_entity_df, !<condition>)
# print(dim(glb_entity_df))

script_df <- rbind(script_df,
                   data.frame(chunk_label="cleanse_data", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##    chunk_label chunk_step_major chunk_step_minor
## 1  import_data                1                0
## 2 cleanse_data                2                0
```

## Step `2`: cleanse data

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="inspect_explore_data", 
                              chunk_step_major=max(script_df$chunk_step_major), 
                              chunk_step_minor=1))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
```

### Step `2`.`1`: inspect/explore data

```r
#print(str(glb_entity_df))
#View(glb_entity_df)

# List info gathered for various columns
# <col_name>:   <description>; <notes>

# "Week" - The range of dates represented by this observation, in year/month/day format.
# "ILI" - This column lists the percentage of ILI-related physician visits for the corresponding week.
# "Queries" - This column lists the fraction of queries that are ILI-related for the corresponding week, adjusted to be between 0 and 1 (higher values correspond to more ILI-related search queries).

# Create new features that help diagnostics
#   Convert factors to dummy variables
#   Build splines   require(splines); bsBasis <- bs(training$age, df=3)

add_new_diag_feats <- function(obs_df, obs_twin_df) {
    obs_df <- mutate(obs_df,
#         <col_name>.NA=is.na(<col_name>)

#         <col_name>.fctr=factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))) 
#         <col_name>.fctr=relevel(factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))),
#                                   "<max_n_val>") 

          # This doesn't work - use sapply instead
#         <col_name>.fctr_num=grep(<col_name>, levels(<col_name>.fctr)), 
#         
#         Date.my=as.Date(strptime(Date, "%m/%d/%y %H:%M")),
#         Year=year(Date.my),
#         Month=months(Date.my),
#         Weekday=weekdays(Date.my)
        
        ILI.log=log(ILI)
                        )

    # If levels of a factor are different across obs_df & glb_predct_df; predict.glm fails  
    # Transformations not handled by mutate
#     obs_df$<col_name>.fctr.num <- sapply(1:nrow(obs_df), 
#         function(row_ix) grep(obs_df[row_ix, "<col_name>"],
#                               levels(obs_df[row_ix, "<col_name>.fctr"])))
    
    print(summary(obs_df))
    print(sapply(names(obs_df), function(col) sum(is.na(obs_df[, col]))))
    return(obs_df)
}

glb_entity_df <- add_new_diag_feats(glb_entity_df, glb_predct_df)
```

```
##      Week                ILI            Queries           ILI.log       
##  Length:417         Min.   :0.5341   Min.   :0.04117   Min.   :-0.6272  
##  Class :character   1st Qu.:0.9025   1st Qu.:0.15671   1st Qu.:-0.1026  
##  Mode  :character   Median :1.2526   Median :0.28154   Median : 0.2252  
##                     Mean   :1.6769   Mean   :0.28603   Mean   : 0.3477  
##                     3rd Qu.:2.0587   3rd Qu.:0.37849   3rd Qu.: 0.7221  
##                     Max.   :7.6189   Max.   :1.00000   Max.   : 2.0306  
##    Week     ILI Queries ILI.log 
##       0       0       0       0
```

```r
glb_predct_df <- add_new_diag_feats(glb_predct_df, glb_entity_df)
```

```
##      Week                ILI            Queries          ILI.log       
##  Length:52          Min.   :0.9018   Min.   :0.2390   Min.   :-0.1034  
##  Class :character   1st Qu.:1.1535   1st Qu.:0.2772   1st Qu.: 0.1428  
##  Mode  :character   Median :1.3592   Median :0.3924   Median : 0.3069  
##                     Mean   :1.6638   Mean   :0.4094   Mean   : 0.4139  
##                     3rd Qu.:1.8637   3rd Qu.:0.4874   3rd Qu.: 0.6226  
##                     Max.   :6.0336   Max.   :0.8054   Max.   : 1.7973  
##    Week     ILI Queries ILI.log 
##       0       0       0       0
```

```r
#pairs(subset(glb_entity_df, select=-c(col_symbol)))

#   Histogram of predictor in glb_entity_df & glb_predct_df
# Check for glb_predct_df & glb_entity_df features range mismatches

# Other diagnostics:
# print(subset(glb_entity_df, <col1_name> == max(glb_entity_df$<col1_name>, na.rm=TRUE) & 
#                         <col2_name> <= mean(glb_entity_df$<col1_name>, na.rm=TRUE)))

print(glb_entity_df[which.max(glb_entity_df$ILI),])
```

```
##                        Week      ILI Queries  ILI.log
## 303 2009-10-18 - 2009-10-24 7.618892       1 2.030631
```

```r
print(glb_entity_df[which.max(glb_entity_df$Queries),])
```

```
##                        Week      ILI Queries  ILI.log
## 303 2009-10-18 - 2009-10-24 7.618892       1 2.030631
```

```r
# print(<col_name>_freq_glb_entity_df <- mycreate_tbl_df(glb_entity_df, "<col_name>"))
# print(which.min(table(glb_entity_df$<col_name>)))
# print(which.max(table(glb_entity_df$ILI)))
# print(which.max(table(glb_entity_df$<col1_name>, glb_entity_df$<col2_name>)[, 2]))
# print(table(glb_entity_df$<col1_name>, glb_entity_df$<col2_name>))
# print(table(is.na(glb_entity_df$<col1_name>), glb_entity_df$<col2_name>))
# print(xtabs(~ <col1_name>, glb_entity_df))
# print(xtabs(~ <col1_name> + <col2_name>, glb_entity_df))
# print(<col1_name>_<col2_name>_xtab_glb_entity_df <- 
#   mycreate_xtab(glb_entity_df, c("<col1_name>", "<col2_name>")))
# <col1_name>_<col2_name>_xtab_glb_entity_df[is.na(<col1_name>_<col2_name>_xtab_glb_entity_df)] <- 0
# print(<col1_name>_<col2_name>_xtab_glb_entity_df <- 
#   mutate(<col1_name>_<col2_name>_xtab_glb_entity_df, 
#             <col3_name>=(<col1_name> * 1.0) / (<col1_name> + <col2_name>))) 

# print(<col2_name>_min_entity_arr <- 
#    sort(tapply(glb_entity_df$<col1_name>, glb_entity_df$<col2_name>, min, na.rm=TRUE)))
# print(<col1_name>_na_by_<col2_name>_arr <- 
#    sort(tapply(glb_entity_df$<col1_name>.NA, glb_entity_df$<col2_name>, mean, na.rm=TRUE)))


# Other plots:
print(myplot_histogram(glb_entity_df, "ILI"))
```

```
## stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
```

![](CDC_Google_Flu_files/figure-html/inspect_explore_data_1-1.png) 

```r
# print(myplot_box(df=glb_entity_df, ycol_names="<col1_name>"))
# print(myplot_box(df=glb_entity_df, ycol_names="<col1_name>", xcol_name="<col2_name>"))
# print(myplot_line(subset(glb_entity_df, Symbol %in% c("KO", "PG")), 
#                   "Date.my", "StockPrice", facet_row_colnames="Symbol") + 
#     geom_vline(xintercept=as.numeric(as.Date("2003-03-01"))) +
#     geom_vline(xintercept=as.numeric(as.Date("1983-01-01")))        
#         )
print(myplot_scatter(glb_entity_df, "ILI.log", "Queries", smooth=TRUE))
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![](CDC_Google_Flu_files/figure-html/inspect_explore_data_1-2.png) 

```r
script_df <- rbind(script_df, 
    data.frame(chunk_label="manage_missing_data", 
        chunk_step_major=max(script_df$chunk_step_major), 
        chunk_step_minor=script_df[nrow(script_df), "chunk_step_minor"]+1))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
```

### Step `2`.`2`: manage missing data

```r
# print(sapply(names(glb_entity_df), function(col) sum(is.na(glb_entity_df[, col]))))
# print(sapply(names(glb_predct_df), function(col) sum(is.na(glb_predct_df[, col]))))
# glb_entity_df <- na.omit(glb_entity_df)
# glb_predct_df <- na.omit(glb_predct_df)

script_df <- rbind(script_df, 
    data.frame(chunk_label="encode_retype_data", 
        chunk_step_major=max(script_df$chunk_step_major), 
        chunk_step_minor=script_df[nrow(script_df), "chunk_step_minor"]+1))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
## 5   encode_retype_data                2                3
```

### Step `2`.`3`: encode/retype data

```r
# map_<col_name>_df <- myimport_data(
#     url="<map_url>", 
#     comment="map_<col_name>_df", print_diagn=TRUE)
# map_<col_name>_df <- read.csv(paste0(getwd(), "/data/<file_name>.csv"), strip.white=TRUE)

# glb_entity_df <- mymap_codes(glb_entity_df, "<from_col_name>", "<to_col_name>", 
#     map_<to_col_name>_df, map_join_col_name="<map_join_col_name>", 
#                           map_tgt_col_name="<to_col_name>")
# glb_predct_df <- mymap_codes(glb_predct_df, "<from_col_name>", "<to_col_name>", 
#     map_<to_col_name>_df, map_join_col_name="<map_join_col_name>", 
#                           map_tgt_col_name="<to_col_name>")
    					
# glb_entity_df$<col_name>.fctr <- factor(glb_entity_df$<col_name>, 
#                     as.factor(union(glb_entity_df$<col_name>, glb_predct_df$<col_name>)))
# glb_predct_df$<col_name>.fctr <- factor(glb_predct_df$<col_name>, 
#                     as.factor(union(glb_entity_df$<col_name>, glb_predct_df$<col_name>)))

script_df <- rbind(script_df, 
                   data.frame(chunk_label="extract_features", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
## 5   encode_retype_data                2                3
## 6     extract_features                3                0
```

## Step `3`: extract features

```r
# Create new features that help prediction
ILI.lag.2 <- lag(zoo(glb_entity_df$ILI), -2, na.pad=TRUE)
glb_entity_df[, "ILI.Lag.2"] <- coredata(ILI.lag.2)
ILI.lag.2 <- lag(zoo(glb_predct_df$ILI), -2, na.pad=TRUE)
glb_predct_df[, "ILI.Lag.2"] <- coredata(ILI.lag.2)

glb_predct_df[1, "ILI.Lag.2"] <- glb_entity_df[nrow(glb_entity_df) - 1, 
                                                   "ILI"]
glb_predct_df[2, "ILI.Lag.2"] <- glb_entity_df[nrow(glb_entity_df), 
                                                   "ILI"]

glb_entity_df <- mutate(glb_entity_df,
    ILI.Lag.2.log=log(ILI.Lag.2)
                    )

glb_predct_df <- mutate(glb_predct_df,
    ILI.Lag.2.log=log(ILI.Lag.2)        
                    )

print(sapply(names(glb_entity_df), function(col) sum(is.na(glb_entity_df[, col]))))
```

```
##          Week           ILI       Queries       ILI.log     ILI.Lag.2 
##             0             0             0             0             2 
## ILI.Lag.2.log 
##             2
```

```r
print(sapply(names(glb_predct_df), function(col) sum(is.na(glb_predct_df[, col]))))
```

```
##          Week           ILI       Queries       ILI.log     ILI.Lag.2 
##             0             0             0             0             0 
## ILI.Lag.2.log 
##             0
```

```r
# print(summary(glb_entity_df))
# print(summary(glb_predct_df))

print(myplot_scatter(glb_entity_df, "ILI.Lag.2.log", "ILI.log", smooth=TRUE))
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](CDC_Google_Flu_files/figure-html/extract_features-1.png) 

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="select_features", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
## 5   encode_retype_data                2                3
## 6     extract_features                3                0
## 7      select_features                4                0
```

## Step `4`: select features

```r
print(glb_feats_df <- myselect_features())
```

```
##                          id     cor.y cor.y.abs
## ILI.Lag.2.log ILI.Lag.2.log 0.9214355 0.9214355
## Queries             Queries 0.8420333 0.8420333
```

```r
script_df <- rbind(script_df, 
    data.frame(chunk_label="remove_correlated_features", 
        chunk_step_major=max(script_df$chunk_step_major),
        chunk_step_minor=script_df[nrow(script_df), "chunk_step_minor"]+1))        
print(script_df)
```

```
##                  chunk_label chunk_step_major chunk_step_minor
## 1                import_data                1                0
## 2               cleanse_data                2                0
## 3       inspect_explore_data                2                1
## 4        manage_missing_data                2                2
## 5         encode_retype_data                2                3
## 6           extract_features                3                0
## 7            select_features                4                0
## 8 remove_correlated_features                4                1
```

### Step `4`.`1`: remove correlated features

```r
print(glb_feats_df <- orderBy(~-cor.y, 
                    merge(glb_feats_df, mydelete_cor_features(), all.x=TRUE)))
```

```
## Loading required package: reshape2
```

```
##               ILI.Lag.2.log   Queries
## ILI.Lag.2.log     1.0000000 0.7426629
## Queries           0.7426629 1.0000000
##               ILI.Lag.2.log   Queries
## ILI.Lag.2.log     0.0000000 0.7426629
## Queries           0.7426629 0.0000000
## [1] "cor(ILI.Lag.2.log, Queries)=0.7427"
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](CDC_Google_Flu_files/figure-html/remove_correlated_features-1.png) 

```
## [1] "cor(ILI.log, ILI.Lag.2.log)=0.9214"
## [1] "cor(ILI.log, Queries)=0.8420"
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning in mydelete_cor_features(): Dropping Queries as a feature
```

![](CDC_Google_Flu_files/figure-html/remove_correlated_features-2.png) 

```
##                          id     cor.y cor.y.abs
## ILI.Lag.2.log ILI.Lag.2.log 0.9214355 0.9214355
##              id     cor.y cor.y.abs cor.low
## 1 ILI.Lag.2.log 0.9214355 0.9214355       1
## 2       Queries 0.8420333 0.8420333      NA
```

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="run_models", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##                  chunk_label chunk_step_major chunk_step_minor
## 1                import_data                1                0
## 2               cleanse_data                2                0
## 3       inspect_explore_data                2                1
## 4        manage_missing_data                2                2
## 5         encode_retype_data                2                3
## 6           extract_features                3                0
## 7            select_features                4                0
## 8 remove_correlated_features                4                1
## 9                 run_models                5                0
```

## Step `5`: run models

```r
max_cor_y_x_var <- subset(glb_feats_df, cor.low == 1)[1, "id"]

#   Regression:
if (glb_is_regression) {
    #   Linear:
    myrun_mdl_fn <- myrun_mdl_lm
}    

#   Classification:
if (glb_is_classification) {
    #   Logit Regression:
    myrun_mdl_fn <- myrun_mdl_glm
}    
    
# Highest cor.y
ret_lst <- myrun_mdl_fn(indep_vars_vctr=max_cor_y_x_var,
                        fit_df=glb_entity_df, OOB_df=glb_predct_df)
```

```
## [1] 1.71619
## [1] 0.7989922
## 
## Call:
## lm(formula = reformulate(indep_vars_vctr, response = glb_predct_var), 
##     data = fit_df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.45124 -0.15652 -0.01877  0.11569  0.83034 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    0.02707    0.01249   2.166   0.0308 *  
## ILI.Lag.2.log  0.92104    0.01911  48.196   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.216 on 413 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.849,	Adjusted R-squared:  0.8487 
## F-statistic:  2323 on 1 and 413 DF,  p-value: < 2.2e-16
## 
##           feats n.fit  R.sq.fit  R.sq.OOB Adj.R.sq.fit SSE.fit SSE.OOB
## 1 ILI.Lag.2.log   417 0.8490434 0.7989922    0.8490434 19.2625 1.71619
##   f.score.OOB
## 1          NA
```

```r
# Enhance Highest cor.y model with additions of interaction terms that were 
#   dropped due to high correlations
if (nrow(subset(glb_feats_df, is.na(cor.low))) > 0)
    ret_lst <- myrun_mdl_fn(indep_vars_vctr=c(max_cor_y_x_var, 
        paste(max_cor_y_x_var, subset(glb_feats_df, is.na(cor.low))[, "id"], sep=":")),
                            fit_df=glb_entity_df, OOB_df=glb_predct_df)    
```

```
## [1] 1.132991
## [1] 0.867299
## 
## Call:
## lm(formula = reformulate(indep_vars_vctr, response = glb_predct_var), 
##     data = fit_df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.50449 -0.14370 -0.02792  0.11351  0.84210 
## 
## Coefficients:
##                       Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            0.02013    0.01209   1.665   0.0966 .  
## ILI.Lag.2.log          0.71522    0.04002  17.872  < 2e-16 ***
## ILI.Lag.2.log:Queries  0.47629    0.08225   5.791 1.39e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.2079 on 412 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.8604,	Adjusted R-squared:  0.8597 
## F-statistic:  1270 on 2 and 412 DF,  p-value: < 2.2e-16
## 
##                                  feats n.fit  R.sq.fit  R.sq.OOB
## 2 ILI.Lag.2.log, ILI.Lag.2.log:Queries   417 0.8604061 0.8672990
## 1                        ILI.Lag.2.log   417 0.8490434 0.7989922
##   Adj.R.sq.fit  SSE.fit  SSE.OOB f.score.OOB
## 2    0.8604061 17.81258 1.132991          NA
## 1    0.8490434 19.26250 1.716190          NA
```

```r
# Low correlated X
ret_lst <- myrun_mdl_fn(indep_vars_vctr=subset(glb_feats_df, cor.low == 1)[, "id"],
                        fit_df=glb_entity_df, OOB_df=glb_predct_df)
```

```
## [1] 1.71619
## [1] 0.7989922
## 
## Call:
## lm(formula = reformulate(indep_vars_vctr, response = glb_predct_var), 
##     data = fit_df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.45124 -0.15652 -0.01877  0.11569  0.83034 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    0.02707    0.01249   2.166   0.0308 *  
## ILI.Lag.2.log  0.92104    0.01911  48.196   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.216 on 413 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.849,	Adjusted R-squared:  0.8487 
## F-statistic:  2323 on 1 and 413 DF,  p-value: < 2.2e-16
## 
##                                  feats n.fit  R.sq.fit  R.sq.OOB
## 2 ILI.Lag.2.log, ILI.Lag.2.log:Queries   417 0.8604061 0.8672990
## 1                        ILI.Lag.2.log   417 0.8490434 0.7989922
## 3                        ILI.Lag.2.log   417 0.8490434 0.7989922
##   Adj.R.sq.fit  SSE.fit  SSE.OOB f.score.OOB
## 2    0.8604061 17.81258 1.132991          NA
## 1    0.8490434 19.26250 1.716190          NA
## 3    0.8490434 19.26250 1.716190          NA
```

```r
glb_sel_mdl <- glb_mdl                        

# All X that is not user excluded
ret_lst <- myrun_mdl_fn(indep_vars_vctr=setdiff(setdiff(names(glb_entity_df),
                                                        glb_predct_var),
                                                glb_exclude_vars_as_features),
                        fit_df=glb_entity_df, OOB_df=glb_predct_df)
```

```
## [1] 1.212356
## [1] 0.8580035
## 
## Call:
## lm(formula = reformulate(indep_vars_vctr, response = glb_predct_var), 
##     data = fit_df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.52209 -0.11082 -0.01819  0.08143  0.76785 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   -0.24064    0.01953  -12.32   <2e-16 ***
## Queries        1.25578    0.07910   15.88   <2e-16 ***
## ILI.Lag.2.log  0.65569    0.02251   29.14   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1703 on 412 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.9063,	Adjusted R-squared:  0.9059 
## F-statistic:  1993 on 2 and 412 DF,  p-value: < 2.2e-16
## 
##                                  feats n.fit  R.sq.fit  R.sq.OOB
## 2 ILI.Lag.2.log, ILI.Lag.2.log:Queries   417 0.8604061 0.8672990
## 4               Queries, ILI.Lag.2.log   417 0.9063414 0.8580035
## 1                        ILI.Lag.2.log   417 0.8490434 0.7989922
## 3                        ILI.Lag.2.log   417 0.8490434 0.7989922
##   Adj.R.sq.fit  SSE.fit  SSE.OOB f.score.OOB
## 2    0.8604061 17.81258 1.132991          NA
## 4    0.9063414 11.95110 1.212356          NA
## 1    0.8490434 19.26250 1.716190          NA
## 3    0.8490434 19.26250 1.716190          NA
```

```r
glb_sel_mdl <- glb_mdl

# User specified
ret_lst <- myrun_mdl_fn(indep_vars_vctr=c("Queries"),
                        fit_df=glb_entity_df, OOB_df=glb_predct_df)
```

```
## [1] 6.536337
## [1] 0.2332548
## 
## Call:
## lm(formula = reformulate(indep_vars_vctr, response = glb_predct_var), 
##     data = fit_df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.76003 -0.19696 -0.01657  0.18685  1.06450 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -0.49934    0.03041  -16.42   <2e-16 ***
## Queries      2.96129    0.09312   31.80   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.2995 on 415 degrees of freedom
## Multiple R-squared:  0.709,	Adjusted R-squared:  0.7083 
## F-statistic:  1011 on 1 and 415 DF,  p-value: < 2.2e-16
## 
##                                  feats n.fit  R.sq.fit  R.sq.OOB
## 2 ILI.Lag.2.log, ILI.Lag.2.log:Queries   417 0.8604061 0.8672990
## 4               Queries, ILI.Lag.2.log   417 0.9063414 0.8580035
## 1                        ILI.Lag.2.log   417 0.8490434 0.7989922
## 3                        ILI.Lag.2.log   417 0.8490434 0.7989922
## 5                              Queries   417 0.7090201 0.2332548
##   Adj.R.sq.fit  SSE.fit  SSE.OOB f.score.OOB
## 2    0.8604061 17.81258 1.132991          NA
## 4    0.9063414 11.95110 1.212356          NA
## 1    0.8490434 19.26250 1.716190          NA
## 3    0.8490434 19.26250 1.716190          NA
## 5    0.7090201 37.23121 6.536337          NA
```

```r
# Simplify a model
# fit_df <- glb_entity_df; glb_mdl <- step(<complex>_mdl)


if (glb_is_regression)
    print(myplot_scatter(glb_models_df, "Adj.R.sq.fit", "R.sq.OOB") + 
          geom_text(aes(label=feats), data=glb_models_df, color="NavyBlue", 
                    size=3.5))
```

![](CDC_Google_Flu_files/figure-html/run_models-1.png) 

```r
if (glb_is_classification) {
    plot_models_df <- mutate(glb_models_df, feats.label=substr(feats, 1, 20))
    print(myplot_hbar(df=plot_models_df, xcol_name="feats.label", 
                      ycol_names="f.score.OOB"))
}

script_df <- rbind(script_df, 
                   data.frame(chunk_label="fit_training.all", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##                   chunk_label chunk_step_major chunk_step_minor
## 1                 import_data                1                0
## 2                cleanse_data                2                0
## 3        inspect_explore_data                2                1
## 4         manage_missing_data                2                2
## 5          encode_retype_data                2                3
## 6            extract_features                3                0
## 7             select_features                4                0
## 8  remove_correlated_features                4                1
## 9                  run_models                5                0
## 10           fit_training.all                6                0
```

## Step `6`: fit training.all

```r
print(mdl_feats_df <- myextract_mdl_feats())
```

```
##                Estimate Std. Error  t value          Pr.z            id
## ILI.Lag.2.log 0.6556896 0.02250504 29.13523 4.082764e-102 ILI.Lag.2.log
## Queries       1.2557762 0.07909836 15.87613  1.253513e-44       Queries
##               fit.feat
## ILI.Lag.2.log     TRUE
## Queries           TRUE
```

```r
if (glb_is_regression) {
    ret_lst <- myrun_mdl_lm(indep_vars_vctr=mdl_feats_df$id, fit_df=glb_entity_df)
    glb_sel_mdl <- glb_mdl    
    glb_entity_df[, glb_predct_var_name] <- predict(glb_sel_mdl, newdata=glb_entity_df)
    print(myplot_scatter(glb_entity_df, glb_predct_var, glb_predct_var_name, 
                         smooth=TRUE))    
}    
```

```
## 
## Call:
## lm(formula = reformulate(indep_vars_vctr, response = glb_predct_var), 
##     data = fit_df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.52209 -0.11082 -0.01819  0.08143  0.76785 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   -0.24064    0.01953  -12.32   <2e-16 ***
## ILI.Lag.2.log  0.65569    0.02251   29.14   <2e-16 ***
## Queries        1.25578    0.07910   15.88   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1703 on 412 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.9063,	Adjusted R-squared:  0.9059 
## F-statistic:  1993 on 2 and 412 DF,  p-value: < 2.2e-16
## 
##                                  feats n.fit  R.sq.fit  R.sq.OOB
## 2 ILI.Lag.2.log, ILI.Lag.2.log:Queries   417 0.8604061 0.8672990
## 4               Queries, ILI.Lag.2.log   417 0.9063414 0.8580035
## 1                        ILI.Lag.2.log   417 0.8490434 0.7989922
## 3                        ILI.Lag.2.log   417 0.8490434 0.7989922
## 5                              Queries   417 0.7090201 0.2332548
## 6               ILI.Lag.2.log, Queries   417 0.9063414        NA
##   Adj.R.sq.fit  SSE.fit  SSE.OOB f.score.OOB
## 2    0.8604061 17.81258 1.132991          NA
## 4    0.9063414 11.95110 1.212356          NA
## 1    0.8490434 19.26250 1.716190          NA
## 3    0.8490434 19.26250 1.716190          NA
## 5    0.7090201 37.23121 6.536337          NA
## 6    0.9063414 11.95110       NA          NA
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
```

```
## Warning: Removed 2 rows containing missing values (stat_smooth).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](CDC_Google_Flu_files/figure-html/fit_training.all-1.png) 

```r
if (glb_is_classification) {
    ret_lst <- myrun_mdl_glm(indep_vars_vctr=mdl_feats_df$id, fit_df=glb_entity_df)
    glb_sel_mdl <- glb_mdl        
    glb_entity_df[, glb_predct_var_name] <- (predict(glb_sel_mdl, 
                        newdata=glb_entity_df, type="response") >= 0.5) * 1.0
    print(xtabs(reformulate(paste(glb_predct_var, glb_predct_var_name, sep=" + ")),
                glb_entity_df))                        
}    

print(glb_feats_df <- mymerge_feats_Pr.z())
```

```
##              id     cor.y cor.y.abs cor.low          Pr.z
## 1 ILI.Lag.2.log 0.9214355 0.9214355       1 4.082764e-102
## 2       Queries 0.8420333 0.8420333      NA  1.253513e-44
```

```r
# Most of this code is used again in predict_newdata chunk
glb_analytics_diag_plots <- function(obs_df) {
    for (var in subset(glb_feats_df, Pr.z < 0.1)$id) {
        plot_df <- melt(obs_df, id.vars=var, 
                        measure.vars=c(glb_predct_var, glb_predct_var_name))
#         if (var == "<feat_name>") print(myplot_scatter(plot_df, var, "value", 
#                                              facet_colcol_name="variable") + 
#                       geom_vline(xintercept=<divider_val>, linetype="dotted")) else     
            print(myplot_scatter(plot_df, var, "value", facet_colcol_name="variable"))
    }
    
    if (glb_is_regression) {
        plot_vars_df <- subset(glb_feats_df, Pr.z < 0.1)
        print(myplot_prediction_regression(obs_df, 
                    ifelse(nrow(plot_vars_df) > 1, plot_vars_df$id[2], ".rownames"), 
                                           plot_vars_df$id[1])
#               + facet_wrap(reformulate(plot_vars_df$id[2])) # if [1,2] is a factor                                                         
#               + geom_point(aes_string(color="<col_name>.fctr")) #  to color the plot
              )
    }    
    
    if (glb_is_classification) {
        plot_vars_df <- subset(glb_feats_df, Pr.z < 0.1)
        print(myplot_prediction_classification(obs_df, 
                    ifelse(nrow(plot_vars_df) > 1, plot_vars_df$id[2], ".rownames"),
                                               plot_vars_df$id[1])
#               + geom_hline(yintercept=<divider_val>, linetype = "dotted")
                )
    }    
}
glb_analytics_diag_plots(glb_entity_df)
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](CDC_Google_Flu_files/figure-html/fit_training.all-2.png) 

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](CDC_Google_Flu_files/figure-html/fit_training.all-3.png) 

```
##                        Week      ILI   Queries   ILI.log ILI.Lag.2
## 282 2009-05-24 - 2009-05-30 4.213152 0.2948207 1.4382111  2.281301
## 296 2009-08-30 - 2009-09-05 3.719694 0.4116866 1.3136413  1.641071
## 213 2008-01-27 - 2008-02-02 4.433810 0.4143426 1.4892593  2.359343
## 295 2009-08-23 - 2009-08-29 2.471660 0.3466135 0.9048899  1.161419
## 281 2009-05-17 - 2009-05-23 3.815720 0.3253652 1.3391293  2.437022
##     ILI.Lag.2.log ILI.log.predict ILI.log.predict.err
## 282     0.8247459       0.6703656           0.7678456
## 296     0.4953493       0.6011410           0.7125003
## 213     0.8583831       0.8425139           0.6467454
## 295     0.1496424       0.2927474           0.6121426
## 281     0.8907770       0.7520185           0.5871109
##                      .label
## 282 2009-05-24 - 2009-05-30
## 296 2009-08-30 - 2009-09-05
## 213 2008-01-27 - 2008-02-02
## 295 2009-08-23 - 2009-08-29
## 281 2009-05-17 - 2009-05-23
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Warning: Removed 2 rows containing missing values (geom_text).
```

![](CDC_Google_Flu_files/figure-html/fit_training.all-4.png) 

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="predict_newdata", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##                   chunk_label chunk_step_major chunk_step_minor
## 1                 import_data                1                0
## 2                cleanse_data                2                0
## 3        inspect_explore_data                2                1
## 4         manage_missing_data                2                2
## 5          encode_retype_data                2                3
## 6            extract_features                3                0
## 7             select_features                4                0
## 8  remove_correlated_features                4                1
## 9                  run_models                5                0
## 10           fit_training.all                6                0
## 11            predict_newdata                7                0
```

## Step `7`: predict newdata

```r
if (glb_is_regression)
    glb_predct_df[, glb_predct_var_name] <- predict(glb_sel_mdl, 
                                        newdata=glb_predct_df, type="response")

if (glb_is_classification)
    glb_predct_df[, glb_predct_var_name] <- (predict(glb_sel_mdl, 
                        newdata=glb_predct_df, type="response") >= 0.5) * 1.0
    
myprint_df(glb_predct_df[, c(glb_id_vars, glb_predct_var, glb_predct_var_name)])
```

```
##                      Week   ILI.log ILI.log.predict
## 1 2012-01-01 - 2012-01-07 0.5691173       0.9091598
## 2 2012-01-08 - 2012-01-14 0.4339881       0.8803854
## 3 2012-01-15 - 2012-01-21 0.4993290       0.7612456
## 4 2012-01-22 - 2012-01-28 0.5213484       0.6459596
## 5 2012-01-29 - 2012-02-04 0.6224787       0.6787968
## 6 2012-02-05 - 2012-02-11 0.6227672       0.7332594
##                       Week   ILI.log ILI.log.predict
## 1  2012-01-01 - 2012-01-07 0.5691173       0.9091598
## 19 2012-05-06 - 2012-05-12 0.2682148       0.2741079
## 20 2012-05-13 - 2012-05-19 0.2365878       0.3401674
## 23 2012-06-03 - 2012-06-09 0.1504465       0.2617214
## 39 2012-09-23 - 2012-09-29 0.2247267       0.4168182
## 48 2012-11-25 - 2012-12-01 0.8002047       0.8854771
##                       Week   ILI.log ILI.log.predict
## 47 2012-11-18 - 2012-11-24 0.8349182       0.7140572
## 48 2012-11-25 - 2012-12-01 0.8002047       0.8854771
## 49 2012-12-02 - 2012-12-08 1.0912677       1.1506614
## 50 2012-12-09 - 2012-12-15 1.2809977       1.1695926
## 51 2012-12-16 - 2012-12-22 1.5145266       1.4638369
## 52 2012-12-23 - 2012-12-29 1.7973462       1.6107246
```

```r
if (glb_is_regression) {
    print(sprintf("Total SSE: %0.4f", 
                  sum((glb_predct_df[, glb_predct_var_name] - 
                        glb_predct_df[, glb_predct_var]) ^ 2)))
    print(myplot_scatter(glb_predct_df, glb_predct_var, glb_predct_var_name, 
                         smooth=TRUE))
    
    glb_predct_df[, "ILI.predict"] <- exp(glb_predct_df[, "ILI.log.predict"])
    head(glb_predct_df)
    print(sprintf("Total Exp SSE: %0.4f", 
                  sum((glb_predct_df[, "ILI.predict"] - 
                        glb_predct_df[, "ILI"]) ^ 2)))
    print(sprintf("Exp RMSE: %0.4f", 
                  (sum((glb_predct_df[, "ILI.predict"] - 
                        glb_predct_df[, "ILI"]) ^ 2) / nrow(glb_predct_df)) ^ 0.5))
}                         
```

```
## [1] "Total SSE: 1.2124"
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![](CDC_Google_Flu_files/figure-html/predict_newdata-1.png) 

```
## [1] "Total Exp SSE: 4.5009"
## [1] "Exp RMSE: 0.2942"
```

```r
if (glb_is_classification)
    print(xtabs(reformulate(paste(glb_predct_var, glb_predct_var_name, sep=" + ")),
                glb_predct_df))
    
glb_analytics_diag_plots(glb_predct_df)
```

![](CDC_Google_Flu_files/figure-html/predict_newdata-2.png) ![](CDC_Google_Flu_files/figure-html/predict_newdata-3.png) 

```
##                       Week       ILI   Queries     ILI.log ILI.Lag.2
## 2  2012-01-08 - 2012-01-14 1.5434005 0.4993360  0.43398812  2.124130
## 1  2012-01-01 - 2012-01-07 1.7667069 0.5936255  0.56911730  1.852736
## 40 2012-09-30 - 2012-10-06 1.1918833 0.4581673  0.17553463  1.214153
## 3  2012-01-15 - 2012-01-21 1.6476154 0.5006640  0.49932902  1.766707
## 28 2012-07-08 - 2012-07-14 0.9281519 0.2656042 -0.07455986  1.078713
##    ILI.Lag.2.log ILI.log.predict ILI.predict ILI.log.predict.err
## 2     0.75336227       0.8803854    2.411829           0.4463973
## 1     0.61666323       0.9091598    2.482236           0.3400425
## 40    0.19404677       0.4619494    1.587165           0.2864148
## 3     0.56911730       0.7612456    2.140941           0.2619166
## 28    0.07576853       0.1425795    1.153245           0.2171393
##                     .label
## 2  2012-01-08 - 2012-01-14
## 1  2012-01-01 - 2012-01-07
## 40 2012-09-30 - 2012-10-06
## 3  2012-01-15 - 2012-01-21
## 28 2012-07-08 - 2012-07-14
```

![](CDC_Google_Flu_files/figure-html/predict_newdata-4.png) 

Null Hypothesis ($\sf{H_{0}}$): mpg is not impacted by am_fctr.  
The variance by am_fctr appears to be independent. 

```r
# print(t.test(subset(cars_df, am_fctr == "automatic")$mpg, 
#              subset(cars_df, am_fctr == "manual")$mpg, 
#              var.equal=FALSE)$conf)
```
We reject the null hypothesis i.e. we have evidence to conclude that am_fctr impacts mpg (95% confidence). Manual transmission is better for miles per gallon versus automatic transmission.


```
## R version 3.1.2 (2014-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] reshape2_1.4.1  zoo_1.7-11      plyr_1.8.1      doBy_4.5-13    
## [5] survival_2.38-1 ggplot2_1.0.0  
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-5 digest_0.6.8     evaluate_0.5.5   formatR_1.0     
##  [5] grid_3.1.2       gtable_0.1.2     htmltools_0.2.6  knitr_1.9       
##  [9] labeling_0.3     lattice_0.20-30  MASS_7.3-39      Matrix_1.1-5    
## [13] munsell_0.4.2    proto_0.3-10     Rcpp_0.11.4      rmarkdown_0.5.1 
## [17] scales_0.2.4     splines_3.1.2    stringr_0.6.2    tcltk_3.1.2     
## [21] tools_3.1.2      yaml_2.1.13
```
