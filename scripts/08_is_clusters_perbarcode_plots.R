#!/usr/bin/R

# alorenzetti 20190826
# requires samtools installed

# set up ####
library(pacman)

packs = c(
  "ggplot2",
  "gtools",
  "tidyverse",
  "scales",
  "viridis",
  "ggtext",
  "Rmisc",
  "ggpubr")

# loading packages
p_load(char = packs)

theme_set(theme_bw())

# defining color palette
# defaultPal=hue_pal()(6)
# defaultPal[1]="#00BA38"
# defaultPal[3]="#F8766D"

# number of nts used to compute flanking coverage of a region (25 for each side)
flank=25
low=0.1
mid=0.5

# setting sniffles dir
snifflesdir="../20190826_sniffles_ref"
#snifflesdir="20190826_sniffles_ref"

# loading files ####
dfins = read.delim(
  paste0(snifflesdir, "/insertionAnnot/insertion_clusters.txt"),
  header = F,
  stringsAsFactors = F
)
dfdel = read.delim(
  paste0(snifflesdir, "/deletionAnnot/deletion_clusters.txt"),
  header = F,
  stringsAsFactors = F
)

# setting colnames
cols = c(
  "strain",
  "cluster",
  "replicon",
  "ISName",
  "ISFamily",
  "meanStart",
  "sdStart",
  "meanLength",
  "sdLength",
  "count",
  "rnames"
)
colnames(dfins) = cols
colnames(dfdel) = cols

# classifying insertion cluster status according to coverage ####
for(i in unique(dfins$strain)){
  df = read.table(paste0(snifflesdir, "/bam/", i, ".txt"), header=F)
  colnames(df) = c("replicon","position","value")
  firstOfBc=head(which(dfins$strain==i), n = 1)
  lastOfBc=tail(which(dfins$strain==i), n = 1)
    
  for(j in firstOfBc:lastOfBc){
    acc=dfins[j,"replicon"]
    pos=seq(dfins[j,"meanStart"]-flank,dfins[j,"meanStart"]+flank)
    meancov=mean(filter(df, replicon == acc & position %in% pos)$value)
    if((dfins[j,"count"] / meancov) <= low){
      dfins[j,"status"] = "rare"
    }else if((dfins[j,"count"] / meancov) > low & (dfins[j,"count"] / meancov) <= mid){
      dfins[j,"status"] = "common"
    }else if((dfins[j,"count"] / meancov) > mid){
      dfins[j,"status"] = "predominant"
    }
  }
}

# classifying deletion cluster status according to coverage ####
for(i in unique(dfdel$strain)){
  df = read.table(paste0(snifflesdir, "/bam/", i, ".txt"), header=F)
  colnames(df) = c("replicon","position","value")
  firstOfBc=head(which(dfdel$strain==i), n = 1)
  lastOfBc=tail(which(dfdel$strain==i), n = 1)
  
  for(j in firstOfBc:lastOfBc){
    acc=dfdel[j,"replicon"]
    pos=seq(dfdel[j,"meanStart"]-flank,dfdel[j,"meanStart"]+flank)
    meancov=mean(filter(df, replicon == acc & position %in% pos)$value)
    if((dfdel[j,"count"] / meancov) <= low){
      dfdel[j,"status"] = "rare"
    }else if((dfdel[j,"count"] / meancov) > low & (dfdel[j,"count"] / meancov) <= mid){
      dfdel[j,"status"] = "common"
    }else if((dfdel[j,"count"] / meancov) > mid){
      dfdel[j,"status"] = "predominant"
    }
  }
}

# setting strain names ####
# two options of names 1 and 2
i=2

dfins$strain[dfins$strain == "barcode01"] = c("NRC1_0p", "&Delta;*ura3* A")[i]
dfins$strain[dfins$strain == "barcode02"] = c("Mutant1_0p", "&Delta;*ura3* B")[i]
dfins$strain[dfins$strain == "barcode03"] = c("Mutant2_0p", "&Delta;*ura3* C")[i]
dfins$strain[dfins$strain == "barcode04"] = c("Mutant3_0p", "&Delta;*ura3* &Delta;*smap1* A")[i]
dfins$strain[dfins$strain == "barcode05"] = c("Mutant1_20p", "&Delta;*ura3* &Delta;*smap1* B")[i]
dfins$strain[dfins$strain == "barcode06"] = c("Mutant2_20p", "&Delta;*ura3* &Delta;*smap1* C")[i]
lvs = levels(as.factor(dfins$strain))
if(i==1){
  lvs = lvs[c(6,1,3,5,2,4)]
}else{
  lvs = lvs[c(1,2,3,4,5,6)]
}
dfins$strain = factor(dfins$strain, levels=lvs)

dfdel$strain[dfdel$strain == "barcode01"] = c("NRC1_0p", "&Delta;*ura3* A")[i]
dfdel$strain[dfdel$strain == "barcode02"] = c("Mutant1_0p", "&Delta;*ura3* B")[i]
dfdel$strain[dfdel$strain == "barcode03"] = c("Mutant2_0p", "&Delta;*ura3* C")[i]
dfdel$strain[dfdel$strain == "barcode04"] = c("Mutant3_0p", "&Delta;*ura3* &Delta;*smap1* A")[i]
dfdel$strain[dfdel$strain == "barcode05"] = c("Mutant1_20p", "&Delta;*ura3* &Delta;*smap1* B")[i]
dfdel$strain[dfdel$strain == "barcode06"] = c("Mutant2_20p", "&Delta;*ura3* &Delta;*smap1* C")[i]
lvs = levels(as.factor(dfdel$strain))
if(i==1){
  lvs = lvs[c(6,1,3,5,2,4)]
}else{
  lvs = lvs[c(1,2,3,4,5,6)]
}
dfdel$strain = factor(dfdel$strain, levels=lvs)

# adjusting isfamily names ####
dfins$ISFamily = sub("_ssgr.*", "", dfins$ISFamily)
dfins$ISFamily = str_replace(dfins$ISFamily, "IS(.*)$", "IS*\\1*")

dfins[dfins$ISFamily != "IS*H3*" & dfins$ISFamily != "IS*4*","ISFamily"] = "Outras famílias"
lvs = levels(as.factor(dfins$ISFamily))
lvs = lvs[c(1,2,3)]
dfins$ISFamily = factor(dfins$ISFamily, levels=lvs)

dfdel$ISFamily = sub("_ssgr.*", "", dfdel$ISFamily)
dfdel$ISFamily = str_replace(dfdel$ISFamily, "IS(.*)$", "IS*\\1*")

dfdel[dfdel$ISFamily != "IS*H3*" & dfdel$ISFamily != "IS*4*","ISFamily"] = "Outras famílias"
lvs = levels(as.factor(dfdel$ISFamily))
lvs = lvs[c(1,2,3)]
dfdel$ISFamily = factor(dfdel$ISFamily, levels=lvs)

# adjusting ISnames and their levels levels of factor ISname
dfins$ISName = str_replace(dfins$ISName, "IS(.*)$", "IS*\\1*")

lvs = mixedsort(levels(as.factor(dfins$ISName)))
dfins$ISName = factor(dfins$ISName, levels=rev(lvs))

dfdel$ISName = str_replace(dfdel$ISName, "IS(.*)$", "IS*\\1*")

lvs = mixedsort(levels(as.factor(dfdel$ISName)))
dfdel$ISName = factor(dfdel$ISName, levels=rev(lvs))

# adjusting levels of factor replicon
lvs=c("NC_002607.1", "NC_001869.1", "NC_002608.1")
dfins$replicon = factor(dfins$replicon, levels=lvs)

# writing dfins and dfdels to file ####
write.table(file = "../dfins.txt", x = dfins, sep="\t", quote = F, row.names = F, col.names = T)
write.table(file = "../dfdel.txt", x = dfdel, sep="\t", quote = F, row.names = F, col.names = T)

# plots ####
# how many IS per barcode ####
# insertions
isCountPerLib = ggplot(dfins, (aes(ISName, fill = ISFamily))) +
  geom_bar() +
  ylim(c(0, 30)) +
  facet_wrap(strain ~ .) +
  coord_flip() +
  ylab("Número de observações") +
  xlab(label = "") +
  scale_fill_viridis(discrete = T,
                     name = "Família:") +
  theme(legend.position = "bottom",
        strip.text = element_markdown(),
        axis.text.y = element_markdown(),
        legend.text = element_markdown())

# saving
ggsave(filename = "plots/isCountPerLib.png",
       plot = isCountPerLib,
       dpi = 300,
       width = 7,
       height = 5)

# deletions
ggplot(dfdel, (aes(ISName, fill = ISFamily))) +
  geom_bar() +
  ylim(c(0, 30)) +
  facet_wrap(strain ~ .) +
  coord_flip() +
  ylab("Número de observações") +
  xlab(label = "") +
  scale_fill_viridis(discrete = T,
                     name = "Família:") +
  theme(legend.position = "bottom")

# size of insertions observed more than 10 times ####
# considering only mean of size of clusters
# df = dfins %>% group_by(ISName) %>% summarise(length(ISName))
# filter = as.character(df[df$`length(ISName)` >= 10, ]$ISName)
# ggplot(filter(dfins, ISName %in% filter), (aes(x=ISName,y=meanLength,fill=ISFamily))) + 
#   geom_violin() +
#   coord_flip() + 
#   xlab(label="") + 
#   scale_fill_viridis(discrete=T)

# classifying IS clusters according to frequency of observations ####
lvs = unique(as.character(dfins$status))
lvs = lvs[c(2,1,3)]
dfins$status = factor(dfins$status, levels=lvs)

ggplot(dfins, (aes(ISName, fill=ISFamily))) + 
  geom_bar() + 
  facet_grid(status ~ .) + 
  coord_flip() + 
  xlab(label="") + 
  scale_fill_viridis(discrete=T)

lvs = unique(as.character(dfdel$status))
dfdel$status = factor(dfdel$status, levels=lvs)

ggplot(dfdel, (aes(ISName, fill=ISFamily))) +
  geom_bar() + facet_grid(status ~ .) +
  coord_flip() + xlab(label="") +
  scale_fill_viridis(discrete=T)

# merging dfins and dfdel to observe frequency status
dfins$svType = "insertion"
dfdel$svType = "excision"
df = rbind.data.frame(dfins, dfdel)

df = df %>%
  mutate(status = case_when(status == "predominant" ~ "Predominante",
                            status == "common" ~ "Comum",
                            status == "rare" ~ "Raro",
                            TRUE ~ as.character(status)),
         svType = case_when(svType == "insertion" ~ "Inserção",
                            svType == "excision" ~ "Excisão",
                            TRUE ~ as.character(status)))

lvs = mixedsort(levels(as.factor(df$ISName)), decreasing = T)
df$ISName = factor(df$ISName, levels=lvs)
df$status = factor(df$status, levels=c("Predominante", "Comum", "Raro"))
df$svType = factor(df$svType, levels=c("Inserção", "Excisão"))

isCountPerStatus = ggplot(df, (aes(ISName, fill=ISFamily))) +
  geom_bar() +
  facet_grid(status ~ svType, scales = "free_x") +
  ylab("Número de observações") +
  coord_flip() +
  xlab(label="") +
  scale_fill_viridis(discrete=T,
                     name = "Família: ") +
  theme(legend.position = "bottom",
        strip.text = element_markdown(),
        axis.text.y = element_markdown(),
        legend.text = element_markdown())

# saving
ggsave(filename = "plots/isCountPerStatus.png",
       plot = isCountPerStatus,
       dpi = 300,
       width = 7,
       height = 8)

# hotspots of insertion ####
# ggplot(dfins, (aes(x=meanStart, y=meanLength, shape=status, col=ISName))) +
#   geom_point(alpha=0.6) + facet_grid(ISFamily ~ replicon, scales = "free_x") +
#   xlab(label="") + scale_color_viridis(discrete=T)
# 
# ggplot(dfins, (aes(x=meanStart, y=meanLength, shape=status, col=ISName))) +
#   geom_point(alpha=0.6) + facet_grid(. ~ replicon, scales = "free_x") +
#   xlab(label="") + scale_color_viridis(discrete=T)

# counting clusters and comparing strains ####
# creating function to parse BAMs
# counting only aligned and
# non supplementary reads
countBamReads = function(x){
  out = system2(command = "samtools",
                args = paste("view -F 0x4 -F 0x800", x, "| wc -l"),
                stdout = T) %>% 
    as.numeric()
  
  return(out)
}

# creating tibble to store results
# and performing counting
bamFiles = paste0("../20190826_sniffles_ref/bam/",
                  paste0("barcode0", c(4:6,1:3)),
                  ".bam")
resCounts = tibble()
for(i in bamFiles){
  resCounts = bind_rows(resCounts,
                        tibble(lib = i,
                               counts = countBamReads(i)))
}

# creating a tibble containing the sum of insertions and deletions
# and normalizing it by read depth
# all cases
insDelCounts = df %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(readCount = resCounts$counts) %>% 
  dplyr::mutate(norm = (max(readCount)/readCount)*count) %>% 
  dplyr::mutate(strain = str_replace(strain, " .$", "")) %>% 
  dplyr::group_by(strain) %>% 
  dplyr::summarise(mean = mean(norm),
                   sd = sd(norm),
                   n = n())

insDelCounts$margin = 1.96*(insDelCounts$sd/sqrt(insDelCounts$n))
insDelCounts$lower95ci = insDelCounts$mean - insDelCounts$margin
insDelCounts$upper95ci = insDelCounts$mean + insDelCounts$margin

# plotting and saving
comparisonPlot = insDelCounts %>% 
  ggplot(aes(x = strain, y = mean)) +
  geom_point() +
  geom_pointrange(aes(ymin = lower95ci,
                      ymax = upper95ci)) +
  coord_flip() +
  xlab("Linhagem") +
  ylab("Média de mobilizações") +
  ggtitle("Total") +
  theme(axis.text.y = element_markdown())

# separated by family ####
fams = df$ISFamily %>% 
  unique() %>% 
  as.character() 
  
isfamdf = list()
isfamplot = list()
isfamplot[["Todas"]] = comparisonPlot
for(i in fams){
  isfamdf[[i]] = df %>%
    dplyr::filter(ISFamily == i) %>% 
    dplyr::group_by(strain) %>%
    dplyr::summarise(count = n()) %>% 
    dplyr::mutate(readCount = resCounts$counts) %>% 
    dplyr::mutate(norm = (max(readCount)/readCount)*count) %>% 
    dplyr::mutate(strain = str_replace(strain, " .$", "")) %>% 
    dplyr::group_by(strain) %>% 
    dplyr::summarise(mean = mean(norm),
                     sd = sd(norm),
                     n = n())
  
  isfamdf[[i]]$margin = 1.96*(isfamdf[[i]]$sd/sqrt(isfamdf[[i]]$n))
  isfamdf[[i]]$lower95ci = isfamdf[[i]]$mean - isfamdf[[i]]$margin
  isfamdf[[i]]$upper95ci = isfamdf[[i]]$mean + isfamdf[[i]]$margin
  
  isfamplot[[i]] = isfamdf[[i]] %>% 
    ggplot(aes(x = strain, y = mean)) +
    geom_point() +
    geom_pointrange(aes(ymin = lower95ci,
                        ymax = upper95ci)) +
    coord_flip() +
    xlab("Linhagem") +
    ylab("Média de mobilizações") +
    ggtitle(i) +
    theme(axis.text.y = element_markdown(),
          plot.title = element_markdown())
}

ggsave(filename = "plots/mobilizationComparisonPerFamily.png",
       plot = ggarrange(plotlist = isfamplot,
                        nrow = isfamplot %>% length(),
                        ncol = 1,
                        labels = "AUTO"),
       dpi = 300,
       width = 7,
       height = 8)

# english version of the aforementioned plot
# plotting general
comparisonPlot = insDelCounts %>% 
  ggplot(aes(x = strain, y = mean)) +
  geom_point() +
  geom_pointrange(aes(ymin = lower95ci,
                      ymax = upper95ci)) +
  coord_flip() +
  xlab("Strain") +
  ylab("Average number of observed mobilizations") +
  ggtitle("Total") +
  theme(axis.text.y = element_markdown())

# separated by family ####
fams = df$ISFamily %>% 
  unique() %>% 
  as.character() 

isfamdf = list()
isfamplot = list()
isfamplot[["Todas"]] = comparisonPlot
for(i in fams){
  if(i != "Outras famílias"){plottitle = i}else{plottitle = "Other families"}
  
  isfamdf[[i]] = df %>%
    dplyr::filter(ISFamily == i) %>% 
    dplyr::group_by(strain) %>%
    dplyr::summarise(count = n()) %>% 
    dplyr::mutate(readCount = resCounts$counts) %>% 
    dplyr::mutate(norm = (max(readCount)/readCount)*count) %>% 
    dplyr::mutate(strain = str_replace(strain, " .$", "")) %>% 
    dplyr::group_by(strain) %>% 
    dplyr::summarise(mean = mean(norm),
                     sd = sd(norm),
                     n = n())
  
  isfamdf[[i]]$margin = 1.96*(isfamdf[[i]]$sd/sqrt(isfamdf[[i]]$n))
  isfamdf[[i]]$lower95ci = isfamdf[[i]]$mean - isfamdf[[i]]$margin
  isfamdf[[i]]$upper95ci = isfamdf[[i]]$mean + isfamdf[[i]]$margin
  
  isfamplot[[i]] = isfamdf[[i]] %>% 
    ggplot(aes(x = strain, y = mean)) +
    geom_point() +
    geom_pointrange(aes(ymin = lower95ci,
                        ymax = upper95ci)) +
    coord_flip() +
    xlab("Strain") +
    ylab("Average number of observed mobilizations") +
    ggtitle(plottitle) +
    theme(axis.text.y = element_markdown(),
          plot.title = element_markdown())
}

ggsave(filename = "plots/mobilizationComparisonPerFamily_en.png",
       plot = ggarrange(plotlist = isfamplot,
                        nrow = isfamplot %>% length(),
                        ncol = 1,
                        labels = "AUTO"),
       dpi = 300,
       width = 7,
       height = 8)
