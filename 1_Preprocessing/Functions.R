# FUNCTIONS 

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
n_round <- 2 # to modify n_round depending on how many decimals are appropriate for data

data_summary_V2 <- function(data, varname, groupnames){
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




# Opposite of %in%
'%ni%' <- Negate('%in%')


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



