#Figure 1 (see file edgecase-trees.pdf)
# made in DensiTree v2.2.5 and Adobe Illustrator, uploaded in cophylo/data
#Figure 2 (see file edgecases-ai_modified.png)
setwd('cophylo/data')
df <- read.table("oldaverageforRscaledL.csv", header=T, sep=',')

nmetric <- nlevels(df$Metric)
df$Distance[grepl("^k", df$Metric)] <- 1 - df$Distance[grepl("^k", df$Metric)]

xrange <- range(df$L)
yrange <- range(df$Distance)

par(mar=c(5,5,1,0), mfrow=c(1,3), cex=1)

plot(xrange, yrange, type="n", 
     xlab="Coalescence rate (lineage pair/Ma)", ylab="Scaled distance", 
     log='x', cex.lab=1.2, cex.axis=1.1)
colors <- rainbow(nmetric, v=0.8)
for (i in 1:nmetric) {
  m <- levels(df$Metric)[i]
  metric <- subset(df, Metric==m)
  lines(metric$L, metric$Distance, type="l", lwd=1.5,
        col=colors[i])
  points(metric$L, metric$Distance, col='white', cex=2, pch=20)
  text(x=metric$L, y=metric$Distance, label=m, cex=0.6, col=colors[i])
}

df <- read.table("oldaverageforRscaledM.csv", header=T, sep=',')
nmetric <- nlevels(df$Metric)
df$Distance[grepl("^k", df$Metric)] <- 1 - df$Distance[grepl("^k", df$Metric)]
# omit M=0
df <- df[df$M>0,]
par(mar=c(5,1,1,1))
plot(range(df$M), range(df$Distance), type="n", 
     xlab="Migration rate (lineage/Ma)", ylab="", 
     log='x', cex.lab=1.2, cex.axis=1.1, yaxt='n')
colors <- rainbow(nmetric, v=0.8)
for (i in 1:nmetric) {
  m <- levels(df$Metric)[i]
  metric <- subset(df, Metric==m)
  lines(metric$M, metric$Distance, type="l", lwd=1.5,
        col=colors[i])
  points(metric$M, metric$Distance, col='white', cex=2, pch=20)
  text(x=metric$M, y=metric$Distance, label=m, cex=0.6, col=colors[i])
}


df <- read.table("oldaverageforRscaledP.csv", header=T, sep=',')
nmetric <- nlevels(df$Metric)
df$Distance[grepl("^k", df$Metric)] <- 1 - df$Distance[grepl("^k", df$Metric)]

plot(range(df$P), range(df$Distance), type="n", 
     xlab="Cospeciation probability", ylab="Scaled distance", 
     cex.lab=1.2, cex.axis=1.1)
colors <- rainbow(nmetric, v=0.8)
for (i in 1:nmetric) {
  m <- levels(df$Metric)[i]
  metric <- subset(df, Metric==m)
  lines(metric$P, metric$Distance, type="l", lwd=1.5,
        col=colors[i])
  points(metric$P, metric$Distance, col='white', cex=2, pch=20)
  text(x=metric$P, y=metric$Distance, label=m, cex=0.6, col=colors[i])
}
#file saved in pdf as cophylo/data/edgecases.pdf, modified in Adobe Illustrator as cophylo/data/edgecases-ai.pdf and redited in GIMP Image Editor as cophylo/data/edgecases-ai_modified.png

#Figure 3
require('corrplot')

avg <- read.table('cophylo/data/averagefinal.csv', sep=',', header=T)
avg$kUn <- 1-avg$kUn
avg$kU <- 1-avg$kU
avg$kLn <- 1-avg$kLn
avg$kL <- 1-avg$kL

# correlation of statistics
#pairs(avg[,4:ncol(avg)], gap=0, cex=0.5, pch=16)

corr <- cor(avg[,4:ncol(avg)], method='spearman')

adjm <- corr>0.95
require(igraph)
g <- graph_from_adjacency_matrix(adjm)
pdf(file='corrplot.pdf', width=5.5, height=5.2)
par(mar=c(5,5,2,2))
corrplot(corr, type='lower', diag=F, method='ellipse', order='FPC', cl.pos='n', tl.col='black')

#Figure 4 

final <- read.csv('cophylo/data/Final.csv', sep=',', header=T)

require(entropy)

# example
d2d <- discretize2d(final$P, final$kUn, numBins1=10, numBins2=10)
plot(d2d)
mi.empirical(d2d)
entropy.empirical(final$kUn)

# generate data frame
temp <- sapply(4:ncol(final), function(i) {
  d2d <- discretize2d(log(final$L), final[,i], numBins1=10,
                      numBins2=10)
  mi.empirical(d2d)
})
mi <- data.frame(metric=names(final)[4:ncol(final)], L=temp)

mi$M <- sapply(4:ncol(final), function(i) {
  d2d <- discretize2d(log(final$M), final[,i], numBins1=10,
                      numBins2=10)
  mi.empirical(d2d)
})
mi$P <- sapply(4:ncol(final), function(i) {
  d2d <- discretize2d(final$P, final[,i], numBins1=10, numBins2=10)
  mi.empirical(d2d)
})


# display the results
require(RColorBrewer)
pal <- brewer.pal(3, 'Set3')
par(mar=c(5,5,1,1))
barplot(t(mi[,2:4]), beside=T, names.arg=mi$metric, 
        las=2, col=pal, ylab='Mutual information')

temp <- final[final$P>0.8,]
mi$M2 <- sapply(4:ncol(temp), function(i) {
  d2d <- discretize2d(log(temp$M), temp[,i], numBins1=10, numBins2=10)
  mi.empirical(d2d)
})

temp <- final[final$M < 1e-4,]
mi$L2 <- sapply(4:ncol(temp), function(i) {
  d2d <- discretize2d(log(temp$L), temp[,i], numBins1=10, numBins2=10)
  mi.empirical(d2d)
})
mi$P2 <- sapply(4:ncol(temp), function(i) {
  d2d <- discretize2d(temp$P, temp[,i], numBins1=10, numBins2=10)
  mi.empirical(d2d)
})

pal <- brewer.pal(6, 'Paired')

par(mfrow=c(3,1), xpd=NA, family='sans')
barplot(rbind(mi$L, mi$L2), beside=T, col=pal[1:2], ylim=c(0, 0.8),
        names.arg=mi$metric, las=2, ylab='Mutual information')
text(x=2, y=0.8, label=expression(Lambda), cex=2)
barplot(rbind(mi$M, mi$M2), beside=T, col=pal[3:4], ylim=c(0, 0.8),
        names.arg=mi$metric, las=2, ylab='Mutual information')
text(x=2, y=0.8, label='M', cex=2)
barplot(rbind(mi$P, mi$P2), beside=T, col=pal[5:6], ylim=c(0, 0.8),
        names.arg=mi$metric, las=2, ylab='Mutual information')
text(x=2, y=0.8, label='P', cex=2)

#Figure 5

# make some contour plots

norm <- function(x) {
  (x-min(x)) / (max(x)-min(x))
}

avg <- read.csv('cophylo/data/average.csv', sep='\t', header=T)


par(mar=c(5,5,1,1), mfrow=c(1,2))

# contour plot of unlabeled kernel on M and P
lo <- loess(kU ~ log10(M) + P, data=avg)
gr <- expand.grid(list(M=10^(seq(-6.9,-0.6,0.1)), P=seq(0.01,0.99,0.01)))
z <- predict(lo, newdata=gr)
#filled.contour(z, col=terrain.colors(20))
#image(z, col=terrain.colors(10,alpha=0.8), breaks=seq(0,1,0.1), xlab='M', ylab='P', xaxt='n')
plot(norm(log10(avg$M)), norm(avg$P), cex=2*sqrt(norm(avg$kU))+0.25, bg=rgb(norm(avg$kU), 0, 1-norm(avg$kU), alpha=0.5), pch=21, col='white', xaxt='n', xlab=expression('log'[10]*' M'), ylab='P', cex.lab=1.2)
title(main='kU', adj=0)
contour(z, add=T, levels=seq(0,1,0.1), labcex=0.8, vfont=NULL)
axis(side=1, at=seq(0,1,0.2), labels=round(seq(-6.9, -0.6, length.out=6),1))


# contour plot of Robinson-Foulds on M and P
lo <- loess(RF ~ log10(M) + P, data=avg)
z <- predict(lo, newdata=gr)
plot(norm(log10(avg$M)), norm(avg$P), cex=2*sqrt(1-norm(avg$RF))+0.25, bg=rgb(1-norm(avg$RF), 0, norm(avg$RF), alpha=0.5), pch=21, col='white', xaxt='n', xlab=expression('log'[10]*' M'), ylab='P', cex.lab=1.2)
title('RF', adj=0)
contour(z, add=T, labcex=0.8, vfont=NULL)
axis(side=1, at=seq(0,1,0.2), labels=round(seq(-6.9, -0.6, length.out=6),1))

#Figure 6
# make some contour plots

norm <- function(x) {
  (x-min(x)) / (max(x)-min(x))
}

avg <- read.csv('cophylo/data/average.csv', sep='\t', header=T)

# contour plot of Sim on M and L
par(mfrow=c(1,2), mar=c(5,5,1,1))

lo <- loess(Sim ~ log10(M) + log10(L), data=avg)
gr <- expand.grid(list(M=10^(seq(-6.9,-0.6,0.1)), L=10^(seq(-6.9,-0.6,0.1))))
z <- predict(lo, newdata=gr)
plot(norm(log10(avg$M)), norm(log10(avg$L)), 
     cex=2*sqrt(1-norm(avg$Sim))+0.25,
     bg=rgb(1-norm(avg$Sim), 0, norm(avg$Sim), alpha=0.5),
     pch=21, col='white', 
     xaxt='n', yaxt='n',
     xlab=expression('log'[10]*' M'), ylab=expression('log'[10]*' '*Lambda), cex.lab=1.2)
title(main='Sim', adj=0)
axis(side=1, at=seq(0,1,0.2), labels=round(seq(-6.9, -0.6, length.out=6),1))
axis(side=2, at=seq(0,1,0.2), labels=round(seq(-6.9, -0.6, length.out=6),1))
contour(z, add=T, labcex=0.8, vfont=NULL)

# same plot with normalized unlabeled kernel
lo <- loess(kUn ~ log10(M) + log10(L), data=avg)
z <- predict(lo, newdata=gr)
plot(norm(log10(avg$M)), norm(log10(avg$L)), 
     cex=2*sqrt(norm(avg$kUn))+0.25,
     bg=rgb(norm(avg$kUn), 0, 1-norm(avg$kUn), alpha=0.5),
     pch=21, col='white', xaxt='n', yaxt='n',
     xlab=expression('log'[10]*' M'), ylab=expression('log'[10]*' '*Lambda), cex.lab=1.2)
title(main='kUn', adj=0)
axis(side=1, at=seq(0,1,0.2), labels=round(seq(-6.9, -0.6, length.out=6),1))
axis(side=2, at=seq(0,1,0.2), labels=round(seq(-6.9, -0.6, length.out=6),1))
contour(z, add=T, labcex=0.8, vfont=NULL)

#Figure 7
setwd('cophylo')

# gives the percentile for 18 pairs in the General data collection
total <- read.csv('data/TotalandKernel.txt', header=T, row.names=1)

# exclude "norm" suffix columns
#total <- total[ , !grepl("norm$", names(total))]

#concordance <- as.factor(c(
 # 'hi',  # 1-2 
#  'hi',  # 3-4
#  'hi',  # 5-6
 # 'hi',  # 7-8
#  'hi',  # 9-10
#  'hi',  # 11-12
#  'lo',  # 15-16
#  'lo',  # 17-18
#  'hi',  # 19-20
 # 'hi',  # 21-22
 # 'lo',  # 23-24
 # 'lo',  # 25-26
 # 'hi',  # 27-28
 # 'lo',  # 29-30
 # 'lo',  # 31-32
 # 'hi',  # 33-34
 # 'hi',  # 35-36
 # 'hi'  # 37-38
#))

# RF dists not normalized in original file
require(ape)
require(phangorn)

total$RF <- sapply(c(seq(1,12,2), seq(15,38,2)), function(i) {
  t1 <- read.tree(paste('data/', i, '.nwk', sep=''))
  t2 <- read.tree(paste('data/', i+1, '.nwk', sep=''))
  RF.dist(t1, t2, normalize=T, rooted=T, check.labels=F)
})

# normalize
normalize <- function(x) {
  (x-min(x)) / (max(x)-min(x))
}
total <- as.data.frame(apply(total, 2, normalize))
total$kU <- 1-total$kU
total$kUn <- 1 - total$kUn
total$kL <- 1-total$kL
total$kLn <- 1-total$kLn

#z <- rep(concordance=='hi', times=ncol(total))


#fit <- glm((concordance=='hi') ~ kUn+kLn+kU+kL+Sim+Node+nPH85, 
           data=as.data.frame(total), family='binomial')
#summary(fit)
#p <- princomp(total)
#biplot(p)


require(e1071)
temp <- as.data.frame(total)

#biplot(dist(temp[,-14]), temp[,14])
#p <- prcomp(temp[,-14])
p <- prcomp(temp[1:18])
#require(devtools)
#install_github('ggbiplot', 'vqv')
#require(ggbiplot)
require(ggplot2)
source('coplylo/scripts/ggbiplot.R')

g <- ggbiplot(p, groups=temp$Group, 
              labels=rownames(temp), labels.size=3, 
              var.col=rgb(0,0,0,0.4))
g <- g + scale_color_manual(name="Group", values=c('firebrick', 'cadetblue'))
g <- g + theme(legend.position='none')
print(g)

#points(p$x[,1]/p$sdev[1], p$x[,2]/p$sdev[2])

#Figure S1

setwd('cophylo')

# gives the percentile for 18 pairs in the General data collection
total <- read.csv('data/TotalandKernelS1.txt', header=T, row.names=1)

# exclude "norm" suffix columns
#total <- total[ , !grepl("norm$", names(total))]

concordance <- as.factor(c(
  'hi',  # 1-2 
  'hi',  # 3-4
  'hi',  # 5-6
  'hi',  # 7-8
  'hi',  # 9-10
  'hi',  # 11-12
  'lo',  # 15-16
  'lo',  # 17-18
  'hi',  # 19-20
  'hi',  # 21-22
  'lo',  # 23-24
  'lo',  # 25-26
  'hi',  # 27-28
  'lo',  # 29-30
  'lo',  # 31-32
  'hi',  # 33-34
  'hi',  # 35-36
  'hi'  # 37-38
))

# RF dists not normalized in original file
require(ape)
require(phangorn)

total$RF <- sapply(c(seq(1,12,2), seq(15,38,2)), function(i) {
  t1 <- read.tree(paste('cophylo/data/', i, '.nwk', sep=''))
  t2 <- read.tree(paste('cophylo/data/', i+1, '.nwk', sep=''))
  RF.dist(t1, t2, normalize=T, rooted=T, check.labels=F)
})

# normalize
normalize <- function(x) {
  (x-min(x)) / (max(x)-min(x))
}
total <- as.data.frame(apply(total, 2, normalize))
total$kU <- 1-total$kU
total$kUn <- 1 - total$kUn
total$kL <- 1-total$kL
total$kLn <- 1-total$kLn

z <- rep(concordance=='hi', times=ncol(total))

#par(mar=c(6,5,1,1))
plot(x=jitter(rep(1:ncol(total), each=nrow(total))), 
     y=unlist(total), 
     bg=ifelse(z, 'salmon', 'dodgerblue'),
     pch=ifelse(z, 21, 4), xaxt='n', xlab='',
     ylab='Normalized distance', cex.lab=1.5)
axis(side=1, at=1:ncol(total), label=colnames(total), las=2)

for (i in seq(2, 18, 2)) {
  rect(i-.5, -0.2, i+.5, 1.2, border=NA, col=rgb(0,0,1,0.1))
}


#Figure S2
library('kernlab')
setwd('cophylo/data/')
m <- read.csv('kUgeneral.csv', header=F, row.names=1)
m1 <- as.matrix(m[,4:39])
km <- as.kernelMatrix(as.matrix(m1))
kp <- kpca(km)
eig(kp) / sum(eig(kp))
plot(rotated(kp), col=(rainbow(2, v=0.8, alpha=0.5)[m$V3]), cex=1.8, pch=c(16,17)[m$V4], xlab='1st principal component', ylab='2nd principal component', cex.lab=1.2, main="")
arrows(1.67293432,0.49707664, -1.74064841,0.56594978,length=0)
arrows(-0.58154352,0.2675036, -0.03083918,0.35049616,length=0)
arrows(-1.35137908,-0.01964648, -1.33928223,0.58599312,length=0)
arrows(-0.66116896,-1.45550015, -0.49407892,-2.02572709,length=0)
arrows(2.32269969,0.57800922, 2.27908574,-0.08308743,length=0)
arrows(-1.29032551,0.5136863, 0.52819849,0.53309431,length=0)
arrows(-0.34562907,-0.81074566, -0.65943773,-1.85093988,length=0)
arrows(-1.67724325,0.33061813, -2.30397687,0.41437632,length=0)
arrows(0.01352481,0.15357856, -0.12274984,-0.56948255,length=0)
arrows(-3.20327918,0.27875119, 2.42196233,0.09256549,length=0)
arrows(-0.41535254,0.13378274, -3.01457914,0.29607081,length=0)
arrows(-1.95673496,0.01149768, -0.84771361,0.51135405,length=0)
arrows(2.38342947,0.26048135, 2.3686219,0.60501541,length=0)
arrows(1.81169161,-2.80157999, 1.34572019,0.07412895,length=0)
arrows(1.49261528,0.5597438, 1.486601,0.98881483,length=0)
arrows(-0.36903701,0.29207771, -0.26587612,0.72319974,length=0)
arrows(-0.69669524,-0.86314656, 1.23204257,0.41086029,length=0)
arrows(2.53856211,0.13050438, -0.53011913,0.32062523,length=0)
text(rotated(kp), label=row.names(m), cex=0.6, pos=1, offset=0.3)
legend("bottomright",legend=c('Host','Pathogen'), pch=c(16,17), pt.cex=2,cex=1,bty='n')

#Figure s3
setwd('cophylo/data/')
library('kernlab')
par(mfrow=c(1,2))
m <- read.csv('PruningOutgroupmidpointtree.csv', header=F, row.names=1)
m1 <- as.matrix(m[,4:41])
km <- as.kernelMatrix(as.matrix(m1))
kp <- kpca(km)
eig(kp) / sum(eig(kp))
plot(rotated(kp), col=(rainbow(2, v=0.8, alpha=0.5)[m$V3]), cex=1.8, pch=c(16,17)[m$V4], xlab='1st principal component', ylab='2nd principal component', cex.lab=1.2, main="")
text(rotated(kp), label=row.names(m), cex=0.6, pos=1, offset=0.3)
arrows(2.019781195,0.178678011,0.7866806,-1.32511831,length=0)
arrows(0.43227247,0.324543497,-1.234696958,0.424404183,length=0)
arrows(1.5599124,2.419102605,-0.437039603,0.892329585,length=0)
arrows(2.208940779,0.159944534,-2.336593081,0.017426068,length=0)
arrows(1.148661862,0.813102006,-2.495362562,0.352393114,length=0)
arrows(0.518067033,-0.598928628,-2.263649855,-1.852131411,length=0)
arrows(1.389366077,-1.717390492,0.557261147,-0.099808197,length=0)
arrows(-1.056203149,-1.601413381,-3.448215013,0.863251236,length=0)
arrows(1.357952222,-0.898278936,2.07195028,-1.318524457,length=0)
arrows(1.868160775,0.552472365,-1.710579243,-1.004708449,length=0)
arrows(1.75398008,1.046695519,-0.913392495,0.09125463,length=0)
arrows(1.75398008,1.046695519,0.200125107,-1.362946069,length=0)
arrows(1.176328865,1.060860211,0.009351812,-0.329388272,length=0)
arrows(3.211175347,0.004074051,-0.926452283,0.953297388,length=0)
arrows(0.582331356,-0.124255464,-2.164851977,1.463439408,length=0)
arrows(0.429695382,1.281753082,-1.795800733,0.300570819,length=0)
arrows(0.713359976,-0.413385806,-0.4854832,-2.061446102,length=0)
arrows(0.511555904,-0.271364612,-1.420498479,0.963127702,length=0)
arrows(0.36853834,-0.512255161,-3.94061046,0.281928217,length=0)
legend("topleft",legend=c('Host','Virus'), pch=c(16,17), pt.cex=2,cex=1,bty='n')
m <- read.csv('PruningOutgroupmidpointtreesigma0.csv', header=F, row.names=1)
m1 <- as.matrix(m[,4:41])
km <- as.kernelMatrix(as.matrix(m1))
kp <- kpca(km)
eig(kp) / sum(eig(kp))
plot(rotated(kp), col=(rainbow(2, v=0.8, alpha=0.5)[m$V3]), cex=1.8, pch=c(16,17)[m$V4], xlab='1st principal component', ylab='2nd principal component', cex.lab=1.2, main="")
text(rotated(kp), label=row.names(m), cex=0.6, pos=1, offset=0.3)
arrows(1.25539257,-1.55780028,0.7916768,-0.48341399,length=0)
arrows(-0.74850331,0.42914732,-1.67480501,2.25613687,length=0)
arrows(1.38222169,-1.64085549,1.39884159,-0.06927178,length=0)
arrows(-1.8219006,-0.15028689,-1.94644061,1.15137236,length=0)
arrows(0.80719723,0.31257261,1.89463019,1.03374031,length=0)
arrows(-3.53444232,-1.85091539,-3.42557183,-1.44179782,length=0)
arrows(-1.43332247,-0.28736728,1.72754098,0.10531758,length=0)
arrows(-1.37013232,1.60181172,-1.55670581,1.52531783,length=0)
arrows(1.59026978,0.92661944,1.32374315,0.95305973,length=0)
arrows(1.03376234,-0.20506895,0.4881895,0.71848067,length=0)
arrows(-0.27360906,-4.2983023,1.48030028,-0.20458272,length=0)
arrows(1.00060224,0.40319618,2.54817113,0.38097715,length=0)
arrows(-1.53326548,1.34124147,-1.53326548,1.34124147,length=0)
arrows(1.22908843,-0.06273605,0.55661017,0.20655741,length=0)
arrows(0.16794474,-0.99201196,-1.2982071,1.49993375,length=0)
arrows(-1.36507648,-1.66487993,2.73911555,0.10069072,length=0)
arrows(-1.27066349,-0.24245418,-1.90821321,1.07296312,length=0)
arrows(1.15223058,1.28233822,0.08839839,1.44269162,length=0)
arrows(-1.53326548,1.34124147,-1.53326548,1.34124147,length=0)
legend("bottomleft",legend=c('Host','Virus'), pch=c(16,17), pt.cex=2,cex=1,bty='n')

#Figure s4
df <- read.table("oldaveragecoeffvarL.csv", header=T, sep=',',na.strings = "NA")
nmetric <- nlevels(df$Metric)


xrange <- range(df$L)
yrange <- range(df$Distance,na.rm = TRUE)

par(mar=c(5,5,1,0), mfrow=c(1,3), cex=1)

plot(xrange, yrange, type="n", 
     xlab="Coalescence rate (lineage pair/Ma)", ylab="Coeffcient of variance", 
     log='x', cex.lab=0.70, cex.axis=1.1)
colors <- rainbow(nmetric, v=0.8)
for (i in 1:nmetric) {
    m <- levels(df$Metric)[i]
    metric <- subset(df, Metric==m)
    lines(metric$L, metric$Distance, type="l", lwd=1.5,
          col=colors[i])
    points(metric$L, metric$Distance, col='white', cex=2, pch=20)
    text(x=metric$L, y=metric$Distance, label=m, cex=0.6, col=colors[i])
}

df <- read.table("oldaveragecoeffvarM.csv", header=T, sep=',',na.strings = "NA")
nmetric <- nlevels(df$Metric)

# omit M=0
df <- df[df$M>0,]
par(mar=c(5,1,1,1))
plot(range(df$M), range(df$Distance,na.rm = TRUE), type="n", 
     xlab="Migration rate (lineage/Ma)", ylab="", 
     log='x', cex.lab=0.70, cex.axis=1.1, yaxt='n')
colors <- rainbow(nmetric, v=0.8)
for (i in 1:nmetric) {
    m <- levels(df$Metric)[i]
    metric <- subset(df, Metric==m)
    lines(metric$M, metric$Distance, type="l", lwd=1.5,
          col=colors[i])
    points(metric$M, metric$Distance, col='white', cex=2, pch=20)
    text(x=metric$M, y=metric$Distance, label=m, cex=0.6, col=colors[i])
}

df <- read.table("oldaveragecoeffvarP.csv", header=T, sep=',',na.strings = "NA")
nmetric <- nlevels(df$Metric)


plot(range(df$P), range(df$Distance,na.rm = TRUE), type="n", 
     xlab="Cospeciation probability", ylab="Coeffcient of variance", 
     cex.lab=0.70, cex.axis=1.1)
colors <- rainbow(nmetric, v=0.8)
for (i in 1:nmetric) {
    m <- levels(df$Metric)[i]
    metric <- subset(df, Metric==m)
    lines(metric$P, metric$Distance, type="l", lwd=1.5,
          col=colors[i])
    points(metric$P, metric$Distance, col='white', cex=2, pch=20)
    text(x=metric$P, y=metric$Distance, label=m, cex=0.6, col=colors[i])
}
