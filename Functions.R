################################################################################
# @Project - Microbial-colonization-programs
# @Description - functions and project configuration
################################################################################

#------------------------------------------------------------------------------
# path declared
#------------------------------------------------------------------------------

## for data stored in separate location
path <- "~/Library/CloudStorage/OneDrive-UniversityofManitoba/CHILD_Nasal_Milk_Project/Analyses/0_Project_Final/" 


#------------------------------------------------------------------------------
# Initial Permutation test for PreTCo (16S data)
#------------------------------------------------------------------------------
### Group_Y shown as healthy group and X as unhealthy or no bm in paper
### N_taxa = number of taxa - first columns in 'change_3m1y_pa_sdf' dataframe
perm_func_dfs <- function(X_df, change_3m1y_pa_sdf, variable, Group_X, Group_Y, N_taxa) {
  
  Group_A <- X_df$X_P_Change
  names(Group_A) <- rep(Group_X, length(Group_A))
  
  Group_B <- X_df$Y_P_Change
  names(Group_B) <- rep(Group_Y, length(Group_B))
  
  original_wilcox <- wilcox.test(Group_A, Group_B, paired = T) # p<0.001
  boxplot(Group_A, Group_B)
  
  shuff_wilcox_p_val <- c()
  diff_list <- list()
  
  set.seed(9999)
  
  it_num <- 1000
  
  for(it in 1:it_num){
    
    change_3m1y_pa_sdf_shuff <- change_3m1y_pa_sdf
    
    ##Permuting the labels
    idx_shuff <- sample(c(1:length(change_3m1y_pa_sdf_shuff[[variable]])))
    shuff_labels <- change_3m1y_pa_sdf_shuff[[variable]][idx_shuff]
    
    change_3m1y_pa_sdf_shuff[[variable]] <- shuff_labels
    
    
    # Make data.frame of prevalence 
    
    X <- change_3m1y_pa_sdf_shuff[change_3m1y_pa_sdf_shuff[[variable]] == Group_X,] %>% .[!is.na(.[[variable]]),]
    Y <- change_3m1y_pa_sdf_shuff[change_3m1y_pa_sdf_shuff[[variable]] == Group_Y,] %>% .[!is.na(.[[variable]]),]
    
    df.pa.nasal.BF_shuff <- data.frame(X_P_Change=colSums(X[1:N_taxa])*100/nrow(X),
                                       Y_P_Change=colSums(Y[1:N_taxa])*100/nrow(Y)) 
    
    
    Group_A_shuff <- df.pa.nasal.BF_shuff$X_P_Change
    names(Group_A_shuff) <- rep("group1", length(Group_A))
    
    Group_B_shuff <- df.pa.nasal.BF_shuff$Y_P_Change
    names(Group_B_shuff) <- rep("group2", length(Group_B_shuff))
    
    shuff_wilcox <- wilcox.test(Group_A_shuff, Group_B_shuff, paired = T) 
    shuff_wilcox_p_val[it] <- shuff_wilcox$p.value
    
    diff_list[[it]] <- df.pa.nasal.BF_shuff$X_P_Change - df.pa.nasal.BF_shuff$Y_P_Change
  }
  
  #length(which(shuff_wilcox_p_val < original_wilcox$p.value))/it_num # 0.012
  
  
  p_val_permut <- c()
  for(k in 1:dim(X_df)[1]){
    
    null_dist <- c()
    for(j in 1:length(diff_list)){
      
      null_dist[j] = diff_list[[j]][k]
      
    }
    real_test_stat <- c(X_df$X_P_Change - X_df$Y_P_Change)[k]
    
    p_val_permut[k] <- length(which(null_dist > real_test_stat))/it 
    
  }
  
  
  names(p_val_permut) <- rownames(X_df)
  sort(1- p_val_permut)  
  
  Res_per_taxa <- data.frame(X_df[,c(1:2)], 
                             X_df$X_P_Change - X_df$Y_P_Change, 
                             1- p_val_permut, p_val_permut)
  names(Res_per_taxa) <- c(Group_X, Group_Y, "original_effect_diff", "p_value", "p_value_2")
  
  return(list(global_res = original_wilcox, per_taxa_res = Res_per_taxa))
  
}

# output: original_wilcox, Res_per_taxa




#------------------------------------------------------------------------------
# functions commonly used
#------------------------------------------------------------------------------

##Summary Function frequently used
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(Mean = mean(x[[col]], na.rm=TRUE),
      SD = sd(x[[col]], na.rm=TRUE), 
      Median = median(x[[col]], na.rm = TRUE), 
      Min = min(x[[col]], na.rm=TRUE), 
      Max = max(x[[col]], na.rm=TRUE), 
      Q1 = quantile(x[[col]], 0.25, type = 7, na.rm=TRUE), 
      Q3 = quantile(x[[col]], 0.75, type = 7, na.rm=TRUE), 
      "N not 0" = sum(x[[col]] != 0, na.rm=TRUE), 
      Prevalence = (sum(x[[col]] != 0, na.rm=TRUE)/NROW(x))*100)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}
#In function use na.rm=TRUE if omiting missing values (NA), otherwise function returns NA if an input is NA
#e.g. usage: data_summary(melt_milk_subset_P1, varname="Relative_Abundance", groupnames=c("Farm", "OTU"))

## Abbreviated Summary Function used for non-normal data

data_summary_V2 <- function(data, varname, groupnames, n_round = 2){
  require(plyr) 
  summary_func <- function(x, col){
    c(mean = round(mean(x[[col]], na.rm=TRUE), n_round),
      median = round(median(x[[col]], na.rm = TRUE), n_round), 
      min = round(min(x[[col]], na.rm=TRUE), n_round), 
      max = round(max(x[[col]], na.rm=TRUE), n_round), 
      Q1 = round(quantile(x[[col]], 0.25, type = 7, na.rm=TRUE), n_round), 
      Q3 = round(quantile(x[[col]], 0.75, type = 7, na.rm=TRUE), n_round), 
      N = NROW(x), 
      "% N not 0" = round((sum(x[[col]] != 0, na.rm=TRUE)/NROW(x))*100,n_round),
      NAs = sum(is.na(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}


prev_summary_2 <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(Prevalence = (sum(x[[col]] != 0, na.rm = TRUE)/sum(!is.na(x[[col]])))*100,
      Count = sum(x[[col]] != 0, na.rm = TRUE),
      N_notNA = sum(!is.na(x[[col]])) 
    )
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}


# function to decrease legend size, - setting defaults to reduce size -can adjust as needed when using addSmallLegend
#https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text 
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 8, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}


# Opposite of %in%
'%ni%' <- Negate('%in%')


## make sure data is binary
make_binary <- function(x) {x2 <- ifelse(x > 0, 1, 0);
return(x2)
}



# Function from CoDaSeq package

codaSeq.clr <- function(x, IQLR=FALSE, aitch=FALSE, samples.by.row=TRUE){
  if(min(x) < 0) stop("only positive real values permitted")
  if (!is.vector(x[,1], mode="numeric") ) stop("counts must be supplied as numbers")
  if ( any( x < 0 ) ) stop("counts cannot be negative")
  if(samples.by.row == TRUE) margin=1
  if(samples.by.row == FALSE) margin=2
  
  if (aitch == FALSE){
    if (IQLR == FALSE){
      if(samples.by.row == T) return( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
      if(samples.by.row == F) return( apply(x, margin, function(x){log(x) - mean(log(x))}) )
    } else if (IQLR == TRUE){
      reads.clr <- t(apply(x, margin, function(x){log(x) - mean(log(x))}))
      reads.var <- apply(reads.clr, 2, var)
      reads.qtl <- quantile(unlist(reads.var))
      mid.set <- which(reads.var < (reads.qtl[4]) & reads.var > (reads.qtl[2]))
      if(samples.by.row == F) return(apply(x, margin, function(x) log(x) - mean(log(x[mid.set]))))
      
      if(samples.by.row == T) return(t(apply(x, margin, function(x) log(x) - mean(log(x[mid.set])))))
    }
  }
  if (aitch == TRUE){
    aitchison.mean <- function( n, log=TRUE ) {
      
      # Input is a vector of non-negative integer counts.
      # Output is a probability vector of expected frequencies.
      # If log-frequencies are requested, the uninformative subspace is removed.
      
      a <- n + 0.5
      sa <- sum(a)
      
      log.p <- digamma(a) - digamma(sa)
      log.p <- log.p - mean(log.p)
      
      if ( log ) return(log.p)
      
      p <- exp( log.p - max(log.p) )
      p <- p / sum(p)
      return(p)
    }
    
    if(samples.by.row == FALSE){
      return(apply(x, margin, aitchison.mean))
    }
    if(samples.by.row == TRUE) x <- t(x)
    return(apply(x, margin, aitchison.mean))
  }
}


# Sum columns, but on a groupwise basis instead of the total column sum
data_sum <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(sum = sum(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}


# # Modified Version of the fromList function in the UpSet package. 
#### need to load UpSet package. This function will show which values/characters/etc. are shared between each set shown in an UpSet plot
fromList2 <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}


## Quicker verions of data summary for binary variables, just showing the "% N not 0" or the "prevalence"
# Only works for binary 1, 0 datasets 

prev_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c("Prevalence" = (sum(x[[col]] != 0, na.rm = TRUE)/sum(!is.na(x[[col]])))*100)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}




#------------------------------------------------------------------------------
# Generate Wilcoxon signed-rank test results (for paired data) and summary table - looped across datasets in a list
## Output dataframe formatted for ggpubr, default p adjustment is BH
## Simplified version, one dataframe per test, one x and one y variable
#------------------------------------------------------------------------------

### column names to test against each other - x_var and y_var (e.g. 3m and 1y, anything paired)

# Extra entries for plotting w/t pubr - defaults just show missing, not necessary
## x_name ## extra, for plotting w/t pubr - must match plot group name
## y_name  ## extra, for plotting w/t pubr - must match plot group name
## variable2 ## Comparison Label (consistent with plot label)

wilcox.pair.df.loops <- function(paired_df_ls, x_var, y_var, 
                                 x_name = NA, y_name = NA, variable2 = NA, # extra if plotting w/t pubr
                                 effect_var = "pa.change.diff",  p.adj.meth = "BH") {
  df_res_ls <- list()
  for(i in seq(paired_df_ls)) {
    test_list <- wilcox.test(paired_df_ls[[i]][[x_var]], paired_df_ls[[i]][[y_var]], paired = TRUE) %>%
      .[["p.value"]] 
    df_res_ls[[i]] <- data.frame(p.value = test_list, 
                                 median_diff =  round(median(paired_df_ls[[i]][[effect_var]]), 1),
                                 ## Extra additions to automatically make it easily plotted 
                                 max = max(c(paired_df_ls[[i]][[x_var]], paired_df_ls[[i]][[y_var]])),
                                 min = min(c(paired_df_ls[[i]][[x_var]], paired_df_ls[[i]][[y_var]])),
                                 n = length(paired_df_ls[[i]][[x_var]]) ## x_var and y_var same length
    )
  }
  names(df_res_ls) <- names(paired_df_ls)
  
  df_res <- bind_rows(df_res_ls, .id = "Taxa") %>% 
    mutate(p.adj = p.adjust(p.value, method=p.adj.meth, n = length(p.value)),
           p.adj.label = signif(p.adj, 3),
           p.cat = ifelse(p.adj < 0.001, "**",
                          ifelse(p.adj < 0.05, "*", "")),
           variable2 = variable2,
           group1 = x_name,
           group2 = y_name
    )
  
  # effect label
  df_res$effect_lab <- paste0(df_res$median_diff, "\nn=", df_res$n)
  
  return(df_res)
}
  

#------------------------------------------------------------------------------
# Wilcoxon signed-rank test results (for paired data) from list of wide dfs
##  default p adjustment is BH
## More complex vs wilcox.pair.df.loops - can enter in multiple different x and y variables, different comparisons 
### doesn't come with effect size output vs. wilcox.pair.df.loops does
#------------------------------------------------------------------------------

## Wilcoxon Signed-Rank test - each variable in var_names_x tested against var_names_y in same position 
multi.wilcox.pair.loops <- function(paired_df_ls, var_names_x, var_names_y, p.adj.meth = "BH", var_names = NA) {
  test_list <- list()
  for(i in seq_along(paired_df_ls)){ ## seq_along(paired_df_ls) - index of each dataframe in paired_df_ls
    test_list[[i]] <- mapply(function(x, y) {wilcox.test(x, y, paired = TRUE)}, 
                             x = as.list(subset(paired_df_ls[[i]], select = var_names_x)), ## won't change row order, matched by subjectid in wide format merge
                             y = as.list(subset(paired_df_ls[[i]], select = var_names_y)), 
                             SIMPLIFY = FALSE) %>% lapply(., function(x) {x[["p.value"]]}) %>% bind_rows() %>% t %>% as.data.frame 
    test_list[[i]] <- test_list[[i]] %>% mutate(p.adj = p.adjust(V1, method=p.adj.meth, n = length(var_names_x)),
                                                variable = var_names,
                                                p.cat = ifelse(p.adj < 0.001, "**", 
                                                               ifelse(p.adj < 0.05, "*", ""))
    )
  }
  ## Order is maintained, adding Lost names 
  names(test_list) <- names(paired_df_ls)
  test_df <- bind_rows(test_list, .id = "Subset")

  return(test_df)
}


#------------------------------------------------------------------------------
# Simple Wilcoxon rank sum test results for unpaired data, from list of wide dfs
## Output dataframe formatted for ggpubr, default p adjustment is BH
## unlike above wilcoxon tests, this requires a vector of continuous and factor vars to test against each other (vs 2 vectors of continuous vars)
#------------------------------------------------------------------------------

## vars_names = 'x' (continuous variable) in wilcox.text (x ~ factor_vars[[i]])
wilcox.unpaired.loops <- function(df_ls, vars_names, factor_vars, p.adj.meth = "BH") {
  test_list <- list()
  for(i in seq_along(df_ls)){ ## seq_along(df_ls) - index of each dataframe in df_ls
    test_list[[i]] <- lapply(as.list(subset(df_ls[[i]], select = vars_names)),
                             function(x) {wilcox.test(x ~ df_ls[[i]][,factor_vars[i]])}) %>% 
      lapply(., function(x) {x[["p.value"]]}) %>% bind_rows() %>% t %>% as.data.frame 
    
    test_list[[i]] <- test_list[[i]] %>% mutate(p.adj = p.adjust(V1, method=p.adj.meth, n = length(vars_names)),
                                                variable = vars_names,
                                                p.cat = ifelse(p.adj < 0.001, "**", 
                                                               ifelse(p.adj < 0.05, "*", "")),
                                                group1 = levels(as.factor(df_ls[[i]][,factor_vars[i]]))[1],
                                                group2 = levels(as.factor(df_ls[[i]][,factor_vars[i]]))[2] ## group1 & group2, extra cols in case of p-value plot annotations 
    )
  }
  ## Order is maintained, adding Lost names 
  names(test_list) <- names(df_ls)
  test_df <- bind_rows(test_list, .id = "Subset")
  
  return(test_df)
}




#------------------------------------------------------------------------------
# Version 2 of Permutation test for PreTCo (metagenomic data)
# not really different from perm_func_dfs, can just use perm_func_dfs in future
#------------------------------------------------------------------------------
### Group_Y shown as healthy group and X as unhealthy or no bm in paper
### N = number of taxa - first columns in 'data' dataframe
perm_func_dfs_v2 <- function(X_df, data, variable, Group_X, Group_Y, N) {
  
  Group_A <- X_df$X_P_Change
  names(Group_A) <- rep(paste(variable, Group_X, sep = "_"), length(Group_A))
  
  Group_B <- X_df$Y_P_Change
  names(Group_B) <- rep(paste(variable, Group_Y, sep = "_"), length(Group_B))
  
  original_wilcox <- wilcox.test(Group_A, Group_B, paired = T)
  boxplot(Group_A, Group_B)
  
  shuff_wilcox_p_val <- c()
  diff_list <- list()
  
  
  set.seed(9999)
  
  it_num <- 1000
  
  for(it in 1:it_num){
    
    change_3m1y_pa_sdf_shuff <- data
    
    ##Permuting the labels
    idx_shuff <- sample(c(1:length(change_3m1y_pa_sdf_shuff[[variable]])))
    shuff_labels <- change_3m1y_pa_sdf_shuff[[variable]][idx_shuff]
    
    change_3m1y_pa_sdf_shuff[[variable]] <- shuff_labels
    
    
    # Make data.frame of prevalence
    X <- change_3m1y_pa_sdf_shuff[change_3m1y_pa_sdf_shuff[[variable]] == Group_X,] %>% .[!is.na(.[[variable]]),]
    Y <- change_3m1y_pa_sdf_shuff[change_3m1y_pa_sdf_shuff[[variable]] == Group_Y,] %>% .[!is.na(.[[variable]]),]
    
    X_df_shuff <- data.frame(X_P_Change=colSums(X[1:N])*100/dim(X)[1],
                             Y_P_Change=colSums(Y[1:N])*100/dim(Y)[1]) 
    
    
    Group_A_shuff <- X_df_shuff$X_P_Change
    names(Group_A_shuff) <- rep(paste(variable, Group_X, sep = "_"), length(Group_A))
    
    Group_B_shuff <- X_df_shuff$Y_P_Change
    names(Group_B_shuff) <- rep(paste(variable, Group_Y, sep = "_"), length(Group_B_shuff))
    
    shuff_wilcox <- wilcox.test(Group_A_shuff, Group_B_shuff, paired = T) 
    shuff_wilcox_p_val[it] <- shuff_wilcox$p.value
    
    diff_list[[it]] <- X_df_shuff$X_P_Change - X_df_shuff$Y_P_Change
  }
  
  length(which(shuff_wilcox_p_val < original_wilcox$p.value))/it_num 
  
  
  p_val_permut <- c()
  for(k in 1:dim(X_df)[1]){
    
    null_dist <- c()
    for(j in 1:length(diff_list)){
      
    null_dist[j] = diff_list[[j]][k]
      
    }
    # rownames(X_df)[k]
    real_test_stat <- c(X_df$X_P_Change - X_df$Y_P_Change)[k]
    
    p_val_permut[k] <- length(which(null_dist > real_test_stat))/it 
    
  }
  
  names(p_val_permut) <- rownames(X_df)
  sort(1- p_val_permut)  
  
  Res_per_taxa <- data.frame(X_df[,c(1:2)], X_df$X_P_Change - X_df$Y_P_Change, 1- p_val_permut, p_val_permut)
  
  
  names(Res_per_taxa) <- c(paste(variable, Group_X, sep = "_"), paste(variable, Group_Y, sep = "_"), 
                           "original_effect_diff", "p_value", "p_value_2")
  
  return(list(global_res = original_wilcox, per_taxa_res = Res_per_taxa))
  
  
}
