library(GenomicRanges)
#Fijense si est? el de la cabra, si no indican los largos a mano
library(BSgenome.Btaurus.UCSC.bosTau8) 
library(stringr)
library(data.table)
library(gtools)
library(reshape2)
library(ggplot2)
library(plyr)

setwd("C:/Users/H/Desktop/UCDavis 2018/Runs of Homozygosity/cgaTOH_cabrasbien")

##Importar salidas del cgaTOH, filtrar datos y generar lista de listas por cromosoma y muestra

ff <- list.files(
           full.names=TRUE, 
           pattern = ".homozygousruns")

# Elimino los de 0.5 Mb
ff = ff[!ff %in% grep('500000', ff, value=T)]

myfilelist <- lapply(ff, read.table, stringsAsFactors = F, sep = ",", header = T)

names(myfilelist) = unlist(lapply(ff, function(x) str_match(x,"[0-9]{1,2}_[0-9]{2}_[0-9]{7,8}")))

vector=rep(1:29, each=5)

crom = split(myfilelist,vector)
rm(vector)
crom_bind = lapply(crom, function(x) rbind(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]]))
remove = c("RTT_11132","RTT_11133","RTT_11151","RTT_11161","RTT_11168","RTT_11172","RTT_11174") #Eliminadas por an?malas
crom_bind = lapply(crom_bind, function(x) x[! x$Label %in% remove, ])

#mismopa = c("RTT_11148","RTT_11134","RTT_11150","RTT_11129")
#crom_bind = lapply(crom_bind, function(x) x[x$Label %in% mismopa, ])

crom_sample = lapply(crom_bind, function(x) split( x , f = x[[1]]))
crom_sample_IR = lapply(crom_sample, function(x) lapply(x, function(x) IRanges(x[[9]],x[[10]])))
crom_sample_Red = lapply(crom_sample_IR, function(x) lapply(x, reduce))
crom_sample_dfH = lapply(crom_sample_Red, function(x) lapply(x, data.frame))

##################

#Para FROHs parciales

#crom_sample_df12H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width < 2000000,]))
#crom_sample_df24H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 2000000 & x$width < 4000000, ]))
#crom_sample_df48H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 4000000 & x$width < 8000000, ]))
#crom_sample_df816H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 8000000 & x$width < 16000000, ]))
#crom_sample_df16H = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width >= 16000000, ]))


crom_sample_sum = lapply(crom_sample_dfH, function(x) lapply(x, function(x) sum(x[[3]])))
#En caso de querer estimar un FROH parcial XX hay que reemplazar "crom_sample_df" por "crom_sample_dfXX"
sumcrom = rbindlist(crom_sample_sum, fill=TRUE)
sumcrom[is.na(sumcrom)] <- 0
sumcrom = as.data.frame(sumcrom)
rownames(sumcrom) = c("01","10","11","12","13","14","15","16","17","18","19","02","20","21","22","23","24","25","26","27","28","29","03","04","05","06","07","08","09")
sumcrom = sumcrom[mixedsort(rownames(sumcrom)),mixedsort(colnames(sumcrom))]

sumcromH = sumcrom

#Repetir para los L
ff <- list.files(path="C:/Users/H/Desktop/RttPEDMAP/Seba_Espana/ROHs_N1_L",
                 full.names=TRUE, pattern = ".homozygousruns")
ff = ff[!ff %in% grep('500000', ff, value=T)] # Elimino los de 0.5 Mb
myfilelist <- lapply(ff, read.table, stringsAsFactors = F, sep = ",", header = T)
names(myfilelist) = unlist(lapply(ff, function(x) str_match(x,"[0-9]{1,2}_[0-9]{2}_[0-9]{7,8}")))
vector=rep(1:29, each=5)
crom = split(myfilelist,vector)
rm(vector)
crom_bind = lapply(crom, function(x) rbind(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]]))
remove = c("RTT_11132","RTT_11133","RTT_11151","RTT_11161","RTT_11168","RTT_11172","RTT_11174") #Eliminadas por an?malas
crom_bind = lapply(crom_bind, function(x) x[! x$Label %in% remove, ])

#mismopa = c("RTT_11148","RTT_11134","RTT_11150","RTT_11129")
#crom_bind = lapply(crom_bind, function(x) x[x$Label %in% mismopa, ])

crom_sample = lapply(crom_bind, function(x) split( x , f = x[[1]]))
crom_sample_IR = lapply(crom_sample, function(x) lapply(x, function(x) IRanges(x[[9]],x[[10]])))
crom_sample_Red = lapply(crom_sample_IR, function(x) lapply(x, reduce))
crom_sample_dfL = lapply(crom_sample_Red, function(x) lapply(x, data.frame))


##################

#Para FROHs parciales

#crom_sample_df12 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width < 2000000,]))
#crom_sample_df24 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 2000000 & x$width < 4000000, ]))
#crom_sample_df48 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 4000000 & x$width < 8000000, ]))
#crom_sample_df816 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 8000000 & x$width < 16000000, ]))
#crom_sample_df16 = lapply(crom_sample_df, function(x) lapply(x, function(x) x[x$width >= 16000000, ]))


crom_sample_sum = lapply(crom_sample_dfL, function(x) lapply(x, function(x) sum(x[[3]])))
#En caso de querer estimar un FROH parcial XX hay que reemplazar "crom_sample_df" por "crom_sample_dfXX"
sumcrom = rbindlist(crom_sample_sum, fill=TRUE)
sumcrom[is.na(sumcrom)] <- 0
sumcrom = as.data.frame(sumcrom)
rownames(sumcrom) = c("01","10","11","12","13","14","15","16","17","18","19","02","20","21","22","23","24","25","26","27","28","29","03","04","05","06","07","08","09")
sumcrom = sumcrom[mixedsort(rownames(sumcrom)),mixedsort(colnames(sumcrom))]

sumcromL = sumcrom

sumcrom = cbind(sumcromH, sumcromL)
#largos = seqlengths(BSgenome.Btaurus.UCSC.bosTau8)
seba = read.table("C:/Users/H/Desktop/UCDavis 2018/Runs of Homozygosity/Largos.txt", stringsAsFactors = F, header = F)
largos = seba$V3

FROHSxChr = sapply(sumcrom, function(x) x/largos[1:29])
FROHSxChr = t(FROHSxChr)
FROHSxChr = as.data.frame(FROHSxChr)
colnames(FROHSxChr) = c(1:29)
FROHSxChr$Group = c(rep("HI",ncol(sumcromH)),rep("LI",ncol(sumcromL)))
FROHSxChr$Sample = rownames(FROHSxChr)

FROHSxChr = melt(FROHSxChr, id.vars=c("Group", "Sample"))
colnames(FROHSxChr) = c("Group","Sample","Chromosome","FROH")


#Graficar FROH por cromosoma
ggplot(FROHSxChr, aes(x=Chromosome, y=FROH, fill=Group)) + theme_classic() + 
  geom_boxplot() + ylab(expression(paste("F"[ROH]))) + 
  stat_summary(fun.y="mean", geom="point", size=1.5,
               position=position_dodge(width=0.75), color="white") + 
  scale_x_discrete(name ="Chromosome")

#Para obtener los FROH de la longitud que hayamos especificado en crom_sample_df"XX"
sums = colSums(sumcrom)
largo = sum(as.numeric(largos [1:29]))
FROHs = sums/largo
FROHs = FROHs[order(names(FROHs))]
FROHs

IDs = read.table("C:/Users/H/Desktop/RttPEDMAP/Seba_Espana/IDs_CRUDOS.txt", header = T, stringsAsFactors = F)
FPeds = read.table("C:/Users/H/Desktop/RttPEDMAP/Seba_Espana/Fped.txt", header = T, stringsAsFactors = F)
IDPeds = merge(IDs, FPeds, by.x="Lab", by.y="Sample")
IDPeds = merge(as.data.frame(FROHs),IDPeds, by.x="row.names", by.y="Sample")
cor.test(IDPeds$FROHs, IDPeds$Fped, method="spearman", exact = F)


#Estadisticas de largos por muestra
keys <- unique(unlist(lapply(crom_sample_dfH, names)))
samkeys = setNames(do.call(mapply, c(FUN=c, lapply(crom_sample_dfH, `[`, keys))), keys)
samvec = unlist(unlist(samkeys, recursive = FALSE))
samvec = samvec[grep("width", names(samvec))]
names(samvec) = substr(names(samvec),1,9) #quedarse con los primeros 9 caracteres
sampROH = split(samvec,names(samvec)) #separo por nombres en lista

statslenH = rbind(Mean=lapply(sampROH, function(x) mean(x)),S.d. = lapply(sampROH, function(x) sd(x)),
                  Minimum = lapply(sampROH, function(x) min(x)), Maximum = lapply(sampROH, function(x) max(x)))
colnames(statslenH) = paste("H",c(1:ncol(statslenH)), sep = "")
statslenH
#write.table(statslenH, "C:/Users/H/Desktop/RttPEDMAP/Regreso/Figuras/Confirmadas/Stats_Length_H.txt", quote = F, col.names = T, row.names = T, sep = "\t")

melt=melt(sampROH)
melt = melt[with(melt, order(L1)), ]
SamplesH = unique(melt$L1)
ReempH = paste("H",c(1:length(SamplesH)),sep = "")
map = setNames(ReempH, SamplesH)
melt[,2] <- map[unlist(melt$L1)]

ggplot(melt, aes(x=L1, y=value)) + theme_bw() + 
  geom_boxplot(color="darkblue") + scale_x_discrete(name ="Sample (Group LI)", limits = mixedsort(unique(melt$L1))) +
  scale_y_continuous(name = "ROH Length", breaks=c(10000000,20000000,30000000,40000000,50000000,60000000,70000000),
                     labels=c("10 Mb","20 Mb", "30 Mb","40 Mb","50 Mb","60 Mb","70 Mb")) +
  stat_summary(fun.y="mean", geom="point", size=1.5,
               position=position_dodge(width=0.75), color="black")


#Estad?sticas de largo por cromosoma
dfcromH = lapply(crom_sample_dfH, function(x) do.call("rbind",x))
names(dfcromH) = c("01","10","11","12","13","14","15","16","17","18","19","02","20","21","22","23","24","25","26","27","28","29","03","04","05","06","07","08","09")
statschrH = rbind(Mean=lapply(dfcromH, function(x) mean(x[[3]])),S.D.=lapply(dfcromH, function(x) sd(x[[3]])), Maximum =lapply(dfcromH, function(x) max(x[[3]])), Minimum=lapply(dfcromH, function(x) min(x[[3]])))
statschrH = statschrH[,mixedsort(colnames(statschrH))]
statschrH
#write.table(statschrH, "C:/Users/H/Desktop/RttPEDMAP/Regreso/Figuras/Confirmadas/Stats_Chr_H.txt", quote = F, col.names = T, row.names = T, sep = "\t")


dfcromH = do.call("rbind", dfcromH)
dfcromH$Chr = as.character(sapply(strsplit(rownames(dfcromH), "\\."), "[", 1))#Reemplazar por lo que haya antes del primer punto

ggplot(dfcromH, aes(x=Chr, y=width)) + theme_bw() + 
  geom_boxplot(color="brown") + scale_x_discrete(name ="Chromosome (Group HI)", limits = mixedsort(unique(dfcromH$Chr))) +
  scale_y_continuous(name = "ROH Length", breaks=c(10000000,20000000,30000000,40000000,50000000,60000000,70000000),
                     labels=c("10 Mb", "20 Mb", "30 Mb","40 Mb","50 Mb","60 Mb","70 Mb")) +
  stat_summary(fun.y="mean", geom="point", size=1.5,
               position=position_dodge(width=0.75), color="black")




## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summaH = summarySE(melt, measurevar = "value", groupvars = "L1")


d <- density(summaH$value)#Promedios
plot(d, main=NA, xlab= "ROH Length (Group HI)", las = 1,xaxt = "n", cex.axis = 0.8)
polygon(d, col="grey", border="black")
axis(1, at=c(0,2000000,4000000,6000000,8000000), labels=c("0 Mb", "2 Mb", "4 Mb", "6 Mb", "8 Mb"))

d <- density(melt$value)#Todos los ROH juntos
plot(d, main=NA, xlab= "ROH Length (Group HI)", las = 1,xaxt = "n", cex.axis = 0.8, xlim=c(0, 8000000))
polygon(d, col="grey", border="black")
axis(1, at=c(0,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000, 16000000), labels=c("0 Mb", "1 Mb", "2 Mb", "3 Mb","4 Mb","5 Mb", "6 Mb","7 Mb", "8 Mb", "16 Mb"))


#Graficar overlaps en un cromosoma

crom_sample24 = crom_sample_dfH$`17`#El cromosoma 24 aparece en el orden 17
setdiff(names(crom_sample$`1`),names(crom_sample$`17`))

#.Agrego las muestras que no tienen ROH
crom_sample24$RTT_11137 = data.frame(start=NA, end =NA, width=NA)
crom_sample24$RTT_11138 = data.frame(start=NA, end =NA, width=NA)
crom_sample24$RTT_11144 = data.frame(start=NA, end =NA, width=NA)
crom_sample24$RTT_11147 = data.frame(start=NA, end =NA, width=NA)
crom_sample24$RTT_11156 = data.frame(start=NA, end =NA, width=NA)


intervalos1 = do.call(rbind.data.frame, crom_sample24)
intervalos1$sample = rownames(intervalos1)
intervalos1$sample = sub("\\.[0-9]{1,2}","",intervalos1$sample)
intervalos1 = intervalos1[with(intervalos1, order(sample)), ]
intervalos1 = data.frame(intervalos1$sample,intervalos1$start,intervalos1$end,intervalos1$width)
colnames(intervalos1) = sub("intervalos1.","",colnames(intervalos1))
intervalos1$bin = rep(NA,nrow(intervalos1))
veces = unname(table(intervalos1$sample))
intervalos1$bin = rep(1:length(veces),veces)

#Graficar overlap
ggplot(intervalos1) + scale_x_continuous(breaks=c(10000000,20000000,30000000,40000000,50000000,60000000),
                                         labels=c("10 Mb", "20 Mb", "30 Mb","40 Mb","50 Mb","60 Mb")) + 
  scale_y_continuous(breaks = c(5,10,15,20,25,30,35)) +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = bin, ymax = bin + 0.9), color="black",fill = "grey") + ylab("Sample (Group HI)") + xlab("Position (Mb)") +
  theme_bw()

#### Graficar FROH por categoria de largo


FH = read.table("C:/Users/H/Desktop/RttPEDMAP/Regreso/FROHsParciales_H.txt", stringsAsFactors = F, header = T)
FH$Group = c("H") #Tabla con los FROH parciales un animal por linea
FL = read.table("C:/Users/H/Desktop/RttPEDMAP/Regreso/FROHsParciales_L.txt", stringsAsFactors = F, header = T)
FL$Group = c("L")
FJ = rbind(FH,FL)#Junto las dos tablas
write.table(FJ,"C:/Users/H/Desktop/RttPEDMAP/Regreso/FROHsParciales.txt", col.names = T, row.names = F, quote = F)

FJ = melt(FJ)

FJ$variable <- mapvalues(FJ$variable, from=c("FROH1.2", "FROH2.4", "FROH4.8","FROH8.16","FROH16"), to=c("[1-2 Mb]","[2-4 Mb]","[4-8 Mb]","[8-16 Mb]","[>16 Mb]"))
FJ$Group <- mapvalues(FJ$Group, from=c("H","L"), to=c("HI","LI"))
ggplot(FJ, aes(x=variable, y=value, fill=Group)) + theme_classic() + 
  geom_boxplot() + ylab(expression(paste("F"[ROH]))) + 
  stat_summary(fun.y="mean", geom="point", size=1.5,
               position=position_dodge(width=0.75), color="white") + 
  scale_x_discrete(name ="Length", limits=c("[1-2 Mb]","[2-4 Mb]","[4-8 Mb]","[8-16 Mb]","[>16 Mb]"))

##########

crom_dfH = lapply(crom_sample_dfH, function(x) rbindlist(x))
crom_dfH = lapply(crom_dfH, function(x) { #para completar con 0 los elementos de la lista vacios
  if (nrow(x) == 0) {
    x = data.frame(start=c(0), end=c(0), width=c(0))
  }
  else (x = x)
})

crom_gr = lapply(crom_dfH, function(x) IRanges(x[[1]],x[[2]]))

crom_sample_df8 = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width > 8000000,]))
crom_df8 = lapply(crom_sample_df8, function(x) rbindlist(x))
crom_df8 = lapply(crom_df8, function(x) { #para completar con 0 los elementos de la lista vacios
  if (nrow(x) == 0) {
    x = data.frame(start=c(0), end=c(0))
  }
  else (x = x)
})

crom_gr8 = lapply(crom_df8, function(x) IRanges(start=x[[1]],end=x[[2]]))

crom_sample_df16 = lapply(crom_sample_dfH, function(x) lapply(x, function(x) x[x$width > 16000000,]))
crom_df16 = lapply(crom_sample_df16, function(x) rbindlist(x))
crom_df16 = lapply(crom_df16, function(x) { #para completar con 0 los elementos de la lista vacios
  if (nrow(x) == 0) {
    x = data.frame(start=c(0), end=c(0))
  }
  else (x = x)
})
crom_gr16 = lapply(crom_df16, function(x) IRanges(start=x[[1]],end=x[[2]]))


ffmaps <- list.files(path="C:/Users/H/Desktop/RttPEDMAP/cgaTOH_IGEVET/PEDMAPs-by-Chr/RttPEDMAP-v6-H-by-Chr",
                 full.names=TRUE, pattern = ".map")
maplist = lapply(ffmaps, read.table, stringsAsFactors = F, sep = "\t", header = F)
names(maplist) = c("01","10","11","12","13","14","15","16","17","18",
                      "19","02","20","21","22","23","24","25","26","27",
                      "28","29","03","04","05","06","07","08","09")
colnames = c("Chr","SNP","Zero","Pos")
maplist = lapply(maplist, setNames, colnames)
maplist$`01`$Chr = sub("^","0",maplist$`01`$Chr)
maplist$`02`$Chr = sub("^","0",maplist$`02`$Chr)
maplist$`03`$Chr = sub("^","0",maplist$`03`$Chr)
maplist$`04`$Chr = sub("^","0",maplist$`04`$Chr)
maplist$`05`$Chr = sub("^","0",maplist$`05`$Chr)
maplist$`06`$Chr = sub("^","0",maplist$`06`$Chr)
maplist$`07`$Chr = sub("^","0",maplist$`07`$Chr)
maplist$`08`$Chr = sub("^","0",maplist$`08`$Chr)
maplist$`09`$Chr = sub("^","0",maplist$`09`$Chr)

mapirl = lapply(maplist, function(x) IRanges(x[[4]],x[[4]]))

y = names(mapirl)

maptodos = do.call(rbind.data.frame, maplist)


##################################################
#Overlaps entre H y maps

LH = list()
for (i in 1:29){
  b = countOverlaps(mapirl[[i]], crom_gr[[i]])
  LH[[i]] = b
}

names(LH) = c("01","10","11","12","13","14","15","16","17","18",
              "19","02","20","21","22","23","24","25","26","27",
              "28","29","03","04","05","06","07","08","09")

LH = lapply(LH, function(x) as.data.frame(x))

L2 = do.call(rbind.data.frame, LH)
mapL = merge(maptodos,L2, by = "row.names")
mapL = mapL[with(mapL,order(Chr, Pos)),]
mapL$Row.names = NULL
rownames(mapL) = NULL

ggplot(data=mapL) + geom_bar(aes(x=Pos, y=x,colour = factor(mapL$Chr)), stat="identity") + 
  facet_grid(. ~ Chr, scales = "free") +labs(x="SNP position by chromosome",y="Number of ROH")+
  theme(legend.position="none", axis.text.x = element_blank(),axis.ticks = element_blank()) +
  scale_y_continuous()
