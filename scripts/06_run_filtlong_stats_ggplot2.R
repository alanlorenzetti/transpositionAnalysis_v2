#!/usr/bin/Rscript

library(ggplot2) ; theme_set(theme_bw())
library(tidyverse)
library(ggridges)
library(gridExtra)
library(dplyr)
library(scales)
library(ggthemes)
library(ggpubr)
library(ggtext)

args=commandArgs(trailingOnly=T)
outputdir=args[1]
inputdir=args[1]

df = read.delim(paste0(inputdir, "/stats.txt"), header=T, stringsAsFactors=F)
#df = read.delim("../20190729_filtlong_stats/stats.txt", header=T, stringsAsFactors=F)
k = 2500000
palette=ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[c(3,2,6,1,4,5,10)]
palette=rev(palette)

df$barcode[df$barcode == "barcode01"] = "&Delta;*ura3* A"
df$barcode[df$barcode == "barcode02"] = "&Delta;*ura3* B"
df$barcode[df$barcode == "barcode03"] = "&Delta;*ura3* C"
df$barcode[df$barcode == "barcode04"] = "&Delta;*ura3* &Delta;*smap1* A"
df$barcode[df$barcode == "barcode05"] = "&Delta;*ura3* &Delta;*smap1* B"
df$barcode[df$barcode == "barcode06"] = "&Delta;*ura3* &Delta;*smap1* C"
df$barcode[df$barcode == "unclassified"] = "Não classificado"

df$barcode = factor(df$barcode,
                    levels = rev(c("&Delta;*ura3* A", "&Delta;*ura3* B", "&Delta;*ura3* C",
                                   "&Delta;*ura3* &Delta;*smap1* A", "&Delta;*ura3* &Delta;*smap1* B", "&Delta;*ura3* &Delta;*smap1* C",
                                   "Não classificado")))

dfsum = group_by(df, barcode)
dfsum = summarise(dfsum, throughput=sum(length))
dfsum$barcode = factor(dfsum$barcode,
                       levels = rev(c("&Delta;*ura3* A", "&Delta;*ura3* B", "&Delta;*ura3* C",
                                      "&Delta;*ura3* &Delta;*smap1* A", "&Delta;*ura3* &Delta;*smap1* B", "&Delta;*ura3* &Delta;*smap1* C",
                                      "Não classificado")))

p1 = ggplot(df,aes(x=as.factor(barcode), fill=barcode)) +
  geom_bar(color="black") + ylim(c(0,25000)) +
  ylab("Número de reads") +
  coord_flip() +
  theme(text=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.y = element_markdown()) +
  scale_fill_manual(values = palette,
                    guide=F)

p2 = ggplot(df, aes(y=as.factor(barcode), x=log10(length),stat="identity", fill=barcode)) +
  geom_density_ridges2() +
  xlab("Log<sub>10</sub>(Tamanho)") +
  xlim(c(1,5)) +
  theme(text=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_markdown(),
        plot.ma = unit(c(5.5,5.5,5.5,20), "points")) +
  scale_fill_manual(values = palette, guide=F)

p3 = ggplot(df, aes(y=as.factor(barcode), x=meanQuality, stat="identity", fill=barcode)) +
  geom_density_ridges2() +
  xlab("Qualidade média") +
  xlim(c(40,100)) +
  theme(text=element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown()) +
  scale_fill_manual(values = palette, guide=F)

p4 = ggplot(dfsum, aes(x=barcode, y=throughput/k, fill=barcode)) +
  geom_col(color="black") + ylim(c(0,100)) +
  ylab(label="Cobertura estimada") +
  coord_flip() +
  theme(text=element_text(size=12),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        plot.ma = unit(c(5.5,5.5,5.5,20), "points")) +
  scale_fill_manual(values = palette, guide=F)

#png(paste0(outputdir, "/stats_barcoded.png"), width=1400, height=900, units="px")
if(!dir.exists("plots")){dir.create("plots")}

ggsave(filename = "plots/sequencingSummary.png",
       plot = ggarrange(plotlist = list(p1, p2, p3, p4),
                        widths = c(1.25,1),
                        labels = "AUTO"),
       width = 8,
       height = 7,
       dpi = 300)

