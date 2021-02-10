#!/usr/bin/Rscript

library(ggplot2) ; theme_set(theme_bw())
library(tidyverse)
library(ggridges)
library(gridExtra)
library(dplyr)
library(scales)
library(ggthemes)
library(ggpubr)

args=commandArgs(trailingOnly=T)
outputdir=args[1]
inputdir=args[1]

df = read.delim(paste0(inputdir, "/stats.txt"), header=T, stringsAsFactors=F)
#df = read.delim("../20190729_filtlong_stats/stats.txt", header=T, stringsAsFactors=F)
k = 2500000
palette=ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[c(3,2,6,1,4,5,10)]
palette=rev(palette)

df$barcode[df$barcode == "barcode01"] = "Δura3 A"
df$barcode[df$barcode == "barcode02"] = "Δura3 B"
df$barcode[df$barcode == "barcode03"] = "Δura3 C"
df$barcode[df$barcode == "barcode04"] = "Δura3Δlsm A"
df$barcode[df$barcode == "barcode05"] = "Δura3Δlsm B"
df$barcode[df$barcode == "barcode06"] = "Δura3Δlsm C"
df$barcode[df$barcode == "unclassified"] = "Não classificado"

df$barcode = factor(df$barcode,
                    levels = rev(c("Δura3 A", "Δura3 B", "Δura3 C",
                                   "Δura3Δlsm A", "Δura3Δlsm B", "Δura3Δlsm C",
                                   "Não classificado")))

dfsum = group_by(df, barcode)
dfsum = summarise(dfsum, throughput=sum(length))
dfsum$barcode = factor(dfsum$barcode,
                       levels = rev(c("Δura3 A", "Δura3 B", "Δura3 C",
                                      "Δura3Δlsm A", "Δura3Δlsm B", "Δura3Δlsm C",
                                      "Não classificado")))

p1 = ggplot(df,aes(x=as.factor(barcode), fill=barcode)) +
  geom_bar(color="black") + ylim(c(0,25000)) +
  ylab("Número de reads") +
  coord_flip() +
  theme(text=element_text(size=12),
        axis.title.y=element_blank()) +
  scale_fill_manual(values = palette,
                    guide=F) +
  ggtitle("A")

p2 = ggplot(df, aes(y=as.factor(barcode), x=log10(length),stat="identity", fill=barcode)) +
  geom_density_ridges2() +
  xlab(bquote(Log[10](Tamanho))) +
  xlim(c(1,5)) +
  theme(text=element_text(size=12),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_manual(values = palette, guide=F) +
  ggtitle("B")

p3 = ggplot(df, aes(y=as.factor(barcode), x=meanQuality, stat="identity", fill=barcode)) +
  geom_density_ridges2() +
  xlab("Qualidade média") +
  xlim(c(40,100)) +
  theme(text=element_text(size=12),
        axis.title.y=element_blank()) +
  scale_fill_manual(values = palette, guide=F) +
  ggtitle("C")

p4 = ggplot(dfsum, aes(x=barcode, y=throughput/k, fill=barcode)) +
  geom_col(color="black") + ylim(c(0,100)) +
  ylab(label="Cobertura estimada") +
  coord_flip() +
  theme(text=element_text(size=12),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_manual(values = palette, guide=F) +
  ggtitle("D")

#png(paste0(outputdir, "/stats_barcoded.png"), width=1400, height=900, units="px")
ggarrange(plotlist = list(p1, p2, p3, p4),
          widths = c(1.25,1))
#dev.off()

