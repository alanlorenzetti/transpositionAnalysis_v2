#/usr/bin/R

# alorenzetti 20190727
# set up
library(dplyr)
library(stringr)

# getting args
args=commandArgs(trailingOnly = T)

# setting parameters
dist=50
filename=args[1]
prefix=sub(".txt","",filename)

# reading dataframe
df = read.delim(filename, sep="\t", header=F, stringsAsFactors = F)

# getting all possible types of IS to iterate
istypes=unique(df[,3])

# parse line function
parseLine = function(x){
  cap=sub("^(.*?):(.*?)-.*$","\\1:\\2",x)
  cap=unlist(strsplit(cap,split=":"))
  acc=cap[1]
  start=cap[2]
  return(c(acc, start))
}

# starting iteration to find clusters
table=NULL
for(i in istypes){
  subdf=df[df[,3] == i,]
  n=1
  flag=0
  if(nrow(subdf) == 1){
    subdf[1,5] = paste0(i, ":", n)
  }else{
    for(j in 2:nrow(subdf)){
      if(flag == 0){
        capNow=parseLine(subdf[j,1])
        accNow=capNow[1] ; startNow=as.integer(capNow[2])
        
        capPast=parseLine(subdf[j-1,1])
        accPast=capPast[1] ; startPast=as.integer(capPast[2])
      }else{
        capNow=parseLine(subdf[j,1])
        accNow=capNow[1] ; startNow=as.integer(capNow[2])
      }
      
      if(accNow == accPast){
        diff=startNow-startPast
        if(diff <= dist){
          subdf[j,5] = paste0(i, ":", n)
          subdf[j-1,5] = paste0(i, ":", n)
          startPast=mean(c(startNow,startPast))
          flag=1
        }else{
          if(j == 2){
            subdf[j-1,5] = paste0(i, ":", n)
            n=n+1
            subdf[j,5] = paste0(i, ":", n)
            flag=0
          }else{
            n=n+1
            subdf[j,5] = paste0(i, ":", n)
            flag=0
          }
        }
      }else{
        if(j == 2){
          subdf[j-1,5] = paste0(i, ":", n)
          n=n+1
          subdf[j,5] = paste0(i, ":", n)
          flag=0
        }else{
          n=n+1
          subdf[j,5] = paste0(i, ":", n)
          flag=0
        }
      }
    }
  }
  table=rbind(table,subdf)
}
table=as.data.frame(table, stringsAsFactors=F)

# parsing generated table
regex="^(.*):(.*)-(.*):(.*):(.*):(.*):(.*)$"
df2 = data.frame(barcode=table$V2,
                 acc=sub(regex,"\\1",table$V1),
                 start=as.integer(sub(regex,"\\2",table$V1)),
                 length=as.integer(sub(regex,"\\5",table$V1)),
                 id=sub(regex,"\\4",table$V1),
                 isName=table$V3,
                 isFamily=table$V4,
                 cluster=table$V5,
                 rnames=sub(regex,"\\7",table$V1),
                 stringsAsFactors = F)

# condensing using dplyr
df2 = df2 %>%
  group_by(barcode,cluster,acc,isName,isFamily) %>%
  summarise(startMean=round(mean(start)),
            startSd=round(sd(start)),
            lengthMean=round(mean(length)),
            lengthSd=round(sd(length)),
            count=length(start),
            rnames=paste(rnames,collapse=","))

# ordering df2 by clusters
df2 = df2[str_order(df2$cluster, numeric=T),]

# writing df2
outputfile=paste0(prefix,"_clusters.txt")
write.table(x = df2, file = outputfile, row.names = F, col.names = F, sep = "\t", na = "NA",quote = F)
