## Script for plotting KStar Results

dat <- read.csv("~/Downloads/Post2xMFSNesikExperiments.tsv", 
                header = T, stringsAsFactors = F, sep = '\t')

library(ggplot2)


plot_dat <- function(dat, residue, indices){
  ## Clean Names
  dat <- dat[,c("Sequence", "K..Lower.Bound", "K..Upper.Bound")]
  names(dat) <- c("Sequence", "K_lower_log10", "K_upper_log10")
  
  ## Trim white space
  dat$Sequence <- trimws(dat$Sequence)
  
  ## Ensure numeric
  dat$K_lower_log10 <- as.numeric(dat$K_lower_log10)
  dat$K_upper_log10 <- as.numeric(dat$K_upper_log10)
  dat$K_lower_log10 <- ifelse(is.infinite(dat$K_lower_log10), 0, dat$K_lower_log10)
  dat$K_upper_log10 <- ifelse(is.infinite(dat$K_upper_log10), 0, dat$K_upper_log10)
  
  ## Clean Sequence into Distinct Columns
  split_seq <- read.table(text = gsub("=", " ", dat$Sequence),
                          col.names = c("tag1","B3",
                                        "tag2","B5",
                                        "tag3","B9",
                                        "tag4","B14"),
                          stringsAsFactors = FALSE)[ , c("B3","B5","B9","B14")]
  
  dat <- cbind(dat, split_seq)
  
  ## All-lower == Wild type
  dat$residue_designation <- ifelse(dat[[residue]] == tolower(dat[[residue]]), 'WT', dat[[residue]])
  
  plot_data <- dat[indices,]
  
  p1 <- ggplot(plot_data,
         aes(y = factor(residue_designation, levels = residue_designation))) +
    geom_linerange(aes(xmin = 0, xmax = K_upper_log10+1), size = 1.2) +
    geom_point(aes(x = K_lower_log10), shape = 21, size = 3) +
    geom_point(aes(x = K_upper_log10), shape = 21, size = 3) +
    labs(x = expression(log[10]~K~score),
         y = NULL,
         title = paste0("Lower vs. Upper log10 K* Scores by Sequence\n (MKI Position: ", residue, ')', sep = '')) +
    theme_classic()
    
  return(p1)
}

p1 <- plot_dat(dat, "B14", 1:4)

p2 <- plot_dat(dat, "B9", c(1,5:7))

p3 <- plot_dat(dat, "B5", c(1,8:10))

p4 <- plot_dat(dat, "B3", c(1, 11:nrow(dat)))

library(gridExtra)

grid.arrange(p4, p3, p2, p1, nrow = 2, ncol=2)
