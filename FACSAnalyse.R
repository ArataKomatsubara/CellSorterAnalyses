# library
library(ggplot2)
library(dplyr)
library(tidyr)
library(hrbrthemes)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(glue)
library(flowCore)
require(scales) 

setwd("C:/Users/user/Desktop/FACS/230211")

# Get current working directory
wd <- getwd()

# Get list of all .fcs files in the working directory
fcs_files <- list.files(path = wd, pattern = "*.fcs", full.names = TRUE)


# Remove "_Data Source - 1.fcs" from each string in fcs_files list
fcs_files <- gsub("_Data Source - 1.fcs", "", fcs_files)
fcs_filesNew <- gsub("C:/Users/user/Desktop/FACS/230211/", "", fcs_files)


facsanalyze_fun <- function(samplename) 
{ sampleFile = paste(samplename , "_Data Source - 1.fcs", sep = "")
  fcs1 = read.FCS(sampleFile)
  
  maxx = 10^4.2
  mini = 10^4.15
  
  data <- fcs1@exprs
  df <- as.data.frame(data) 
  
  df = df %>%
    mutate(SampleNo = "s07")
  
  
  
  selectHA <- df[df$`FL2-A` < maxx & df$`FL2-A` > mini ,] #4.45~4.5èÊ
  
  
  vline <- selectHA %>%
    group_by(SampleNo) %>%
    summarise(med_intensity = mean(`FL1-A`), num = sum(`FL1-A` >= 0))
  
  vline <- vline  %>%
    mutate(logmid_int = log10(med_intensity))
  
  p2 <- df %>%
    ggplot(aes(x = log10(`FL1-A`),y= log10(`FL2-A`)))+
    scale_fill_continuous(limits = c(0,450) ,type = "viridis") +
    geom_bin2d(bins =300) + 
    facet_grid(. ~ SampleNo) +
    xlim(0,6) + ylim(0,6) +
    ggtitle(samplename) +
    geom_hline(yintercept = log10(maxx),linetype = 2,color = "brown") +
    geom_hline(yintercept = log10(mini),linetype = 2, color = "brown") +
    theme_bw() 
  plot(p2)
  
  #Correlation coefficient
  
  
  p <- selectHA %>%
    ggplot(aes(x = `FL1-A`)) +
    geom_density(aes(alpha = 0.15),fill = "gray") +
    geom_vline(data = vline, aes(xintercept = med_intensity), color = "red") +
    scale_fill_manual(values=c("gray")) + 
    geom_text(data = vline,aes( x = 10^(logmid_int+0.15) , y= 3 ,label = round (med_intensity,0), color="red"), size= 5) +
    geom_text(data = vline,aes( x = 10 , y= 0.15 ,label = glue("N = {num}"), )) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),limits= c(10000, 320000),labels = trans_format("log10", math_format(10^.x))) +
    ggtitle(samplename) + xlab("AF488") + ylab("Density") +
    theme_bw() +
    theme(legend.position = "none")
  
  
  plot(p)
  
  ggsave(file = paste("2D_",paste(samplename , ".png", sep = ""),sep = ""), plot = p2,width = 8, height = 7)
  ggsave(file = paste("3Ddensity_",paste(samplename , ".png", sep = ""),sep = ""), plot = p,width = 7, height = 3)
}


lapply(fcs_filesNew, facsanalyze_fun)


