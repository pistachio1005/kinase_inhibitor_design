## Script for plotting KStar Results
dat <- read.csv("~/Downloads/Post2xMFSNesikExperiments.tsv", 
                header = T, stringsAsFactors = F, sep = '\t')

library(ggplot2)
library(gridExtra)

plot_dat <- function(dat, residue, indices){
  ## Clean Names
  dat <- dat[,c("Sequence", "K..Lower.Bound", "K..Upper.Bound")]
  names(dat) <- c("Sequence", "K_lower", "K_upper")
  
  ## Trim white space
  dat$Sequence <- trimws(dat$Sequence)
  
  ## Ensure numeric
  dat$K_lower <- as.numeric(dat$K_lower)
  dat$K_upper <- as.numeric(dat$K_upper)
  dat$K_lower <- ifelse(is.infinite(dat$K_lower), 0, dat$K_lower)
  dat$K_upper <- ifelse(is.infinite(dat$K_upper), 0, dat$K_upper)
  
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
  
  x_max <- max(dat$K_upper)
  
  plot_data <- dat[indices,]
  
  p1 <- ggplot(plot_data,
         aes(y = factor(residue_designation, levels = residue_designation))) +
    geom_linerange(aes(xmin = K_lower, xmax = K_upper), size = 1.2) +
    geom_point(aes(x = K_lower), shape = 21, size = 3) +
    geom_point(aes(x = K_upper), shape = 21, size = 3) +
    scale_x_continuous(limits = c(-0.05, x_max),
                       expand = expansion(mult = c(0, .02))) +
    labs(x = expression(KStar~score),
         y = NULL,
         title = paste0("Lower vs. Upper K* Scores by Sequence\n (MKI Position: ", residue, ')', sep = '')) +
    theme(
      axis.text  = element_text(face = "bold", colour = "black"),  # tick labels
      axis.title = element_text(face = "bold", colour = "black")   # axis titles
    ) +
    theme_classic()
    
  return(p1)
}

p1 <- plot_dat(dat, "B14", 1:4)

p2 <- plot_dat(dat, "B9", c(1,5:7))

p3 <- plot_dat(dat, "B5", c(1,8:10))

p4 <- plot_dat(dat, "B3", c(1, 11:nrow(dat)))

par(mar = c(5, 10, 2, 1),
    font.lab  = 2,   # axis titles
    font.axis = 2)   # tick labels
## font.lab / font.axis = 2 ->  bold.  (1 = plain, 3 = italic, 4 = bold-italic)


grid.arrange(p4, p3, p2, p1, nrow = 2, ncol = 2)
