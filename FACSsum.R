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

#解析するディレクトリ名
directoryname = "C:/Users/user/Desktop/FACS/230113_15F11mut"

#HAの提示量をどの範囲から指定するのかを設定
maxx = 10^4.2
mini = 10^4.15

#ボックスプロットの最小値
minbox = 4.4

setwd(directoryname)

# ディレクトリ内のCSVファイルのリストを取得する
path_to_directory <- directoryname
file_list <- list.files(path_to_directory, pattern = "*.fcs", full.names = TRUE)

# ファイルを読み込んでsampleNoをヘッダーとしてファイル名を追加する
data_list <- lapply(file_list, function(file) {
  fcs <- read.FCS(file)
  fcsdata <- fcs@exprs
  data <- as.data.frame(fcsdata) 
  samplename <- gsub("_Data Source - 1.fcs", "", file)
  data$sampleNo <- basename(samplename)
  colnames(data)[1] <- "old_sampleNo"
  colnames(data)[ncol(data)] <- "sampleNo"
  return(data)
})

# データフレームを縦に結合する
combined_data <- do.call(rbind, data_list)

#HAの提示量で特定の範囲のデータを選択、本コード一番上を参照
selectHA <- combined_data[combined_data$`FL2-A` < maxx & combined_data$`FL2-A` > mini ,] 


#平均値を取っておく(一応作ったがボックスとかバイオリンの場合は使わないかも？)
meanFL1A <- selectHA %>%
  group_by(sampleNo) %>%
  summarise(med_intensity = mean(`FL1-A`), num = sum(`FL1-A` >= 0))

LogVline <- meanFL1A  %>%
  mutate(logmid_int = log10(med_intensity))


#まずは適当にボックスプロットしてみる
p <- selectHA %>% 
  ggplot(aes(x = log10(`FL1-A`), y = sampleNo, fill = sampleNo)) +
  stat_boxplot(geom='errorbar', width = 0.5) +
  geom_boxplot(width = 0.8, color = "Black", outlier.colour = NA) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 1.5, color = "red", fill = "red") +
  geom_text(data = meanFL1A ,aes(x = minbox +0.05, y = sampleNo, label = round(med_intensity, 0)), color = "red", size = 3)+
  theme_bw() + 
  theme(legend.position = "none") +
  xlim(minbox, 5) +
  ylab("SampleNo") +
  xlab("AF488_Intensity(log10X)") +
  scale_fill_manual(values = rep("lightgray", length(unique(selectHA$sampleNo)))) # 全てのカテゴリにlightgrayを使用

plot(p)

ggsave("myplot.png", plot = p, width = 5, height = 5, units = "in")
