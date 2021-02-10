# alorenzetti 20191004

####setup####
# set up
library(circlize)
library(rtracklayer)
library(scales)
library(viridis)
library(ggthemes)

# save svg?
save="n"
ht=10
wh=10

# min, mean and max gc should be computed for entire genome?
# if "n" the program will compute min, mean and max separately for chr
# and the combination of plasmids
genomewidegc="y"

# default pal to match viridis colors
strainPal=ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[c(3,2,6,1,4,5)]
discretePal=viridis_pal()(6)[c(1,4,6)]
#defaultPal=hue_pal()(6)
#defaultPal=defaultPal[c(1,4,3,2,5,6)]

# setting parameters
circos.clear()
circos.par("start.degree" = 90,
           "track.height" = 0.075,
           "gap.degree" = c(2,15))


# setting replicon names
chr="NC_002607.1"
pnrc100="NC_001869.1"
pnrc200="NC_002608.1"

# creating a df to represent the layout of genome
df=data.frame(acc=c(chr,pnrc100,pnrc200),
               start=c(0,0,0),
               end=c(2014240,191345,365425),
               foo=c(chr,pnrc100,pnrc200),
               bar=c(chr,pnrc100,pnrc200),
               stringsAsFactors = F)

# saving
if(save=="y"){svg("plasmids.svg", height = ht, width = wh)}
par(cex=1.5)

####init plot with the genome layout for the PLASMIDS####

circos.initializeWithIdeogram(df, plotType = c("labels","axis")[2], chromosome.index = c(chr, pnrc100, pnrc200)[2:3])

#####regions plasmids#####
# initializing regions used in this analysis
chr1=c(pnrc100,pnrc200)
start1=c(0,112795)
end1=c(150252,332792)
regions=(data.frame(chr=chr1,start=start1,end=end1))

# plotting the regions
circos.genomicTrackPlotRegion(regions, ylim = c(0,1), track.height=0.1,
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = c("grey"))
                              },
                              bg.border = NA
)

#####links plasmids#####
# initializing and plotting links
chr1=c(pnrc100,pnrc100)
start1=c(0,158783)
end1=c(112795,191345)
chr2=c(pnrc200,pnrc200)
start2=c(0,332792)
end2=c(112795,365424)
bed1=data.frame(chr=chr1,start=start1,end=end1)
bed2=data.frame(chr=chr2,start=start2,end=end2)
col=adjustcolor(discretePal[2],alpha.f = 0.5)
circos.genomicLink(bed1, bed2, col = col)

chr1=c(pnrc100,pnrc200)
start1=c(32041,32041)
end1=c(71211,64674)
chr2=c(pnrc100,pnrc200)
start2=c(191346,365424)
end2=c(150252,332793)
bed1=data.frame(chr=chr1,start=start1,end=end1)
bed2=data.frame(chr=chr2,start=start2,end=end2)
col=adjustcolor(discretePal[1],alpha.f = 0.5)
circos.genomicLink(bed1, bed2, col = col)

####regions IS plasmids#####
is = as.data.frame(rtracklayer::import("~/gdrive/is-derived-RNAs/Hsalinarum/misc/Hsalinarum-ISSaga-checked.gff3"))
is$rpt_family = as.character(sub("\\+.*$", "", is$rpt_family))
is[is$rpt_family != "IS4" & is$rpt_family != "ISH3","rpt_family"] = "Other_Families"
isfamilies=c("IS4", "ISH3", "Other_Families")
bed_list = vector("list", length = length(isfamilies))
n=1
for(i in isfamilies){
  bed_list[[n]] = is[is$rpt_family == i,c("seqnames", "start", "end")]
  n=n+1
}

# plotting annotated IS regions
circos.genomicTrackPlotRegion(bed_list, ylim = c(0,1), track.height=0.05, bg.col = adjustcolor("white", alpha.f = 0.9),
                              panel.fun = function(region, value, ...) {
                                i=getI(...)
                                circos.genomicRect(region, value, col = discretePal[i])
                              }
)

# initializing and plotting arrows
circos.arrow(sector.index = pnrc100, track.index = 2, x1 = 0, x2 = 112795, col = discretePal[2])
circos.arrow(sector.index = pnrc100, track.index = 2, x1 = 158783, x2 = 191345, col = discretePal[2])
circos.arrow(sector.index = pnrc200, track.index = 2, x1 = 0, x2 = 112795, col = discretePal[2])
circos.arrow(sector.index = pnrc200, track.index = 2, x1 = 332793, x2 = 365424, col = discretePal[2])

circos.arrow(sector.index = pnrc100, track.index = 2, x1 = 32041, x2 = 71210, arrow.position = "start", col = discretePal[1])
circos.arrow(sector.index = pnrc100, track.index = 2, x1 = 150252, x2 = 191345, col = discretePal[1])
circos.arrow(sector.index = pnrc200, track.index = 2, x1 = 332793, x2 = 365424, col = discretePal[1])
circos.arrow(sector.index = pnrc200, track.index = 2, x1 = 32041, x2 = 64674, arrow.position = "start", col = discretePal[1])

# plotting approximate regions of origins of replication Coker2009
circos.rect(xleft = 93976, xright = 95443, ybottom = 0, ytop = 1, sector.index = pnrc100, track.index = 2, col="black")
circos.rect(xleft = 93976, xright = 95443, ybottom = 0, ytop = 1, sector.index = pnrc200, track.index = 2, col="black")
circos.rect(xleft = 288858, xright = 290880, ybottom = 0, ytop = 1, sector.index = pnrc200, track.index = 2, col="black")

#####gc content plasmid#####
# reading gccontent and plotting
gc = read.delim(file = "../20190420_computegc_ref/gccontent.bedgraph.gz", header=F)
if(genomewidegc == "y"){
  mingc=min(gc$V4)
  maxgc=max(gc$V4)
  meangc=mean(gc$V4)
}else{
  mingc=min(gc[gc$V1 == pnrc100 | gc$V1 == pnrc200,4])
  maxgc=max(gc[gc$V1 == pnrc100 | gc$V1 == pnrc200,4])
  meangc=mean(gc[gc$V1 == pnrc100 | gc$V1 == pnrc200,4])
}
col=colorRamp2(c(mingc,meangc,maxgc), c(strainPal[6], "black", strainPal[1]))
circos.genomicHeatmap(gc, col = col, numeric.column = 4, connection_height = 0.0001, heatmap_height = 0.075)

#####coverage plasmids#####
# reading coverage and plotting
barcodes=paste0("barcode0", 1:6)
n=1
bed_list=vector("list", length(barcodes))
for(i in barcodes){
  df1 = read.delim(paste0("../20190826_sniffles_ref/bam/", i,"_500bp.bedgraph"), header=F, na.strings = ".")
  df1 = df1[!is.na(df1$V4),]
  df1[df1$V1 == pnrc200,2] = df1[df1$V1 == pnrc200,2] + 112795
  df1[df1$V1 == pnrc200,3] = df1[df1$V1 == pnrc200,3] + 112795
  bed_list[[n]] = df1
  n=n+1
}

color=adjustcolor("white",alpha.f = 0.9)
circos.genomicTrackPlotRegion(bed_list, track.height=0.15, bg.col=color,
                                panel.fun = function(region, value, ...) {
                                  i = getI(...)
                                  color = adjustcolor(strainPal[i],alpha.f = 0.5)
                                  circos.genomicLines(region, value, col=color, lwd=2.5, ...)
                                }
)
circos.yaxis(side = "left", sector.index = pnrc100, track.index = CELL_META$track.index, labels.cex = 0.5)

#####insertions plasmids####
# reading dfins containing information about new insertions
dfins=read.delim("../dfins.txt", header=T)

# unifying IS1595 IS5 and ISH6 in a single family
filter = as.character(dfins$ISFamily) != "IS4" & as.character(dfins$ISFamily) != "ISH3"
dfins[,"ISFamily"] = as.character(dfins[,"ISFamily"])
dfins[filter,"ISFamily"] = "Other_Families"
dfins[,"ISFamily"] = as.factor(dfins[,"ISFamily"])

# adjusting coordinates of pnrc200
dfins[dfins$replicon == pnrc200,"meanStart"] = dfins[dfins$replicon == pnrc200,"meanStart"] + 112795

# adjusting meanLength to numeric
dfins[,"meanLength"] = as.numeric(dfins[,"meanLength"])
# creating a list of dataframes. each data frame an IS family
isfamilies=c("IS4", "ISH3", "Other_Families")
bed_list=vector("list", length(isfamilies))
n=1
for(i in isfamilies){
  bed_list[[n]] = dfins[dfins$ISFamily == i, c("replicon","meanStart","meanStart","meanLength")]
  colnames(bed_list[[n]]) = c("chr","start","end","meanLength")
  n=n+1
}

color=adjustcolor("white",alpha.f = 0.9)
circos.genomicTrackPlotRegion(bed_list, track.height=0.15, bg.col=color,
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                color = adjustcolor(discretePal[i],alpha.f = 1)
                                circos.genomicPoints(region, value, col=color, pch=1, cex=0.3, ...)
                              }
)
circos.yaxis(side = "left", sector.index = pnrc100, track.index = CELL_META$track.index, labels.cex = 0.5)

#####excisions plasmids#####
# reading dfdel containing information about new insertions
dfdel=read.delim("dfdel.txt", header=T)

# unifying IS1595 IS5 and IS66 in a single family
filter = as.character(dfdel$ISFamily) != "IS4" & as.character(dfdel$ISFamily) != "ISH3"
dfdel[,"ISFamily"] = as.character(dfdel[,"ISFamily"])
dfdel[filter,"ISFamily"] = "Other_Families"
dfdel[,"ISFamily"] = as.factor(dfdel[,"ISFamily"])

# adjusting coordinates of pnrc200
dfdel[dfdel$replicon == pnrc200,"meanStart"] = dfdel[dfdel$replicon == pnrc200,"meanStart"] + 112795

# adjusting meanLength to numeric
dfdel[,"meanLength"] = as.numeric(dfdel[,"meanLength"])
# creating a list of dataframes. each data frame an IS family
isfamilies=c("IS4", "ISH3", "Other_Families")
bed_list=vector("list", length(isfamilies))
n=1
for(i in isfamilies){
  bed_list[[n]] = dfdel[dfdel$ISFamily == i, c("replicon","meanStart","meanStart","meanLength")]
  colnames(bed_list[[n]]) = c("chr","start","end","meanLength")
  n=n+1
}

color=adjustcolor("white",alpha.f = 0.9)
circos.genomicTrackPlotRegion(bed_list, track.height=0.15, bg.col=color, ylim = c(600,1600),
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                color = adjustcolor(strainPal[i], alpha.f = 1)
                                circos.genomicPoints(region, value, col=color, pch=1, cex=0.3, ...)
                              }
)
circos.yaxis(side = "left", sector.index = pnrc100, track.index = CELL_META$track.index, labels.col = "#FFFFFF00", labels.cex = 10^-7)

# closing device to save image
if(save=="y"){dev.off()}

####init plot with the genome layout for the CHROMOSOME####

circos.clear()
circos.par("start.degree" = 90,
           "gap.degree" = 15)
par(cex=1.25)

# saving
if(save=="y"){svg("chr.svg", width = wh, height = ht)}
par(cex=1.5)

circos.initializeWithIdeogram(df, plotType = c("labels","axis")[2], chromosome.index = c(chr, pnrc100, pnrc200)[1])

# initializing regions used in this analysis
chr=chr
start=0
end=2014238
regions=(data.frame(chr=chr,start=start,end=end))

####regions chromosome####
# plotting the regions
circos.genomicTrackPlotRegion(regions, ylim = c(0,1), track.height=0.1,
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region, value, col = c("grey"))
                              },
                              bg.border = NA
)

circos.rect(xleft = 37969, xright = 39766, ybottom = 0, ytop = 1, sector.index = chr, track.index = CELL_META$track.index, col="black")
circos.rect(xleft = 846958, xright = 851813, ybottom = 0, ytop = 1, sector.index = chr, track.index = CELL_META$track.index, col="black")
circos.rect(xleft = 1248368, xright = 1253403, ybottom = 0, ytop = 1, sector.index = chr, track.index = CELL_META$track.index, col="black")
circos.rect(xleft = 1798148, xright = 1801564, ybottom = 0, ytop = 1, sector.index = chr, track.index = CELL_META$track.index, col="black")

####regions IS chromosome#####
is = as.data.frame(rtracklayer::import("~/gdrive/is-derived-RNAs/Hsalinarum/misc/Hsalinarum-ISSaga-checked.gff3"))
is$rpt_family = as.character(sub("\\+.*$", "", is$rpt_family))
is[is$rpt_family != "IS4" & is$rpt_family != "ISH3","rpt_family"] = "Other_Families"
isfamilies=c("IS4", "ISH3", "Other_Families")
bed_list = vector("list", length = length(isfamilies))
n=1
for(i in isfamilies){
  bed_list[[n]] = is[is$rpt_family == i,c("seqnames", "start", "end")]
  n=n+1
}

# plotting annotated IS regions
circos.genomicTrackPlotRegion(bed_list, ylim = c(0,1), track.height=0.05, bg.col = adjustcolor("white", alpha.f = 0.9),
                              panel.fun = function(region, value, ...) {
                                i=getI(...)
                                circos.genomicRect(region, value, col = discretePal[i], lwd = 0.15)
                              }
)

####gccontent chromosome####
# plotting gc
if(genomewidegc == "y"){
  mingc=min(gc$V4)
  maxgc=max(gc$V4)
  meangc=mean(gc$V4)
}else{
  mingc=min(gc[gc$V1 == chr,4])
  maxgc=max(gc[gc$V1 == chr,4])
  meangc=mean(gc[gc$V1 == chr,4])
}
col=colorRamp2(c(mingc,meangc,maxgc), c(strainPal[6], "black", strainPal[1]))
circos.genomicHeatmap(gc, col = col, numeric.column = 4, connection_height = 0.0001, heatmap_height = 0.075)

####coverage chromosome####
#reading coverage and plotting
barcodes=paste0("barcode0", 1:6)
n=1
bed_list=vector("list", length(barcodes))
maxchr=0
for(i in barcodes){
  df1 = read.delim(paste0("../20190826_sniffles_ref/bam/", i,"_500bp.bedgraph"), header=F, na.strings = ".")
  df1 = df1[!is.na(df1$V4),]
  bed_list[[n]] = df1
  n=n+1
  if(maxchr < max(df1[df1$V1 == chr,]$V4)){maxchr = max(df1[df1$V1 == chr,]$V4)}
}

circos.genomicTrackPlotRegion(bed_list, track.height=0.15, bg.col="white", ylim=c(0,maxchr),
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                color = adjustcolor(strainPal[i], alpha.f = 0.5)
                                circos.genomicLines(region, value, col=color, lwd=2.5, ...)
                              }
)
circos.yaxis(side = "left", sector.index = chr, track.index = CELL_META$track.index, labels.cex = 0.5)

####insertions chromosome####
# reading dfins containing information about new insertions
dfins=read.delim("../dfins.txt", header=T)

# unifying IS1595 IS5 and ISH6 in a single family
filter = as.character(dfins$ISFamily) != "IS4" & as.character(dfins$ISFamily) != "ISH3"
dfins[,"ISFamily"] = as.character(dfins[,"ISFamily"])
dfins[filter,"ISFamily"] = "Other_Families"
dfins[,"ISFamily"] = as.factor(dfins[,"ISFamily"])

# adjusting coordinates of pnrc200
dfins[dfins$replicon == pnrc200,"meanStart"] = dfins[dfins$replicon == pnrc200,"meanStart"] + 112795

# adjusting meanLength to numeric
dfins[,"meanLength"] = as.numeric(dfins[,"meanLength"])
# creating a list of dataframes. each data frame an IS family
isfamilies=c("IS4", "ISH3", "Other_Families")
bed_list=vector("list", length(isfamilies))
n=1
for(i in isfamilies){
  bed_list[[n]] = dfins[dfins$ISFamily == i, c("replicon","meanStart","meanStart","meanLength")]
  colnames(bed_list[[n]]) = c("chr","start","end","meanLength")
  n=n+1
}

color=adjustcolor("white",alpha.f = 0.9)
circos.genomicTrackPlotRegion(bed_list, track.height=0.15, bg.col=color,
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                color = adjustcolor(discretePal[i],alpha.f = 1)
                                circos.genomicPoints(region, value, col=color, pch=1, cex=0.3, ...)
                              }
)
circos.yaxis(side = "left", sector.index = chr, track.index = CELL_META$track.index, labels.cex = 0.5)

####excisions chromosome####
# reading dfdel containing information about new insertions
dfdel=read.delim("../dfdel.txt", header=T)

# unifying IS1595 IS5 and IS66 in a single family
filter = as.character(dfdel$ISFamily) != "IS4" & as.character(dfdel$ISFamily) != "ISH3"
dfdel[,"ISFamily"] = as.character(dfdel[,"ISFamily"])
dfdel[filter,"ISFamily"] = "Other_Families"
dfdel[,"ISFamily"] = as.factor(dfdel[,"ISFamily"])

# adjusting coordinates of pnrc200
dfdel[dfdel$replicon == pnrc200,"meanStart"] = dfdel[dfdel$replicon == pnrc200,"meanStart"] + 112795

# adjusting meanLength to numeric
dfdel[,"meanLength"] = as.numeric(dfdel[,"meanLength"])
# creating a list of dataframes. each data frame an IS family
isfamilies=c("IS4", "ISH3", "Other_Families")
bed_list=vector("list", length(isfamilies))
n=1
for(i in isfamilies){
  bed_list[[n]] = dfdel[dfdel$ISFamily == i, c("replicon","meanStart","meanStart","meanLength")]
  colnames(bed_list[[n]]) = c("chr","start","end","meanLength")
  n=n+1
}

color=adjustcolor("white",alpha.f = 0.9)
circos.genomicTrackPlotRegion(bed_list, track.height=0.15, bg.col=color,
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                color = adjustcolor(discretePal[i], alpha.f = 1)
                                circos.genomicPoints(region, value, col=color, pch=1, cex=0.3, ...)
                              }
)
circos.yaxis(side = "left", sector.index = chr, track.index = CELL_META$track.index, labels.col = "#FFFFFF00", labels.cex = 10^-7)

# closing device to save img
if(save=="y"){dev.off()}
