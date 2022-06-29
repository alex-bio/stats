## k-means clustering 
## 22-04-20 AP

#### Assign data, fake example below ####
## want to maximize between.ss and minimize within.ss, so let's look at the ratio of (between.ss)/(total.ss). This should approach 1 as number of clusters approaches the sample size of the data.
## first make some data you know clusters to test this. These are means and sd for a Normal distribution, randomly chosen values for each "sample", defined as x,y,z,p,q,r,i,j,k. Play around with the means and sd to see how the distributions overlap and how this affects the clusters.
mean1 <- 0
mean2 <- 5
mean3 <- 10
sd1 <- 1.5
sd2 <- 1.5
sd3 <- 1.5
fake.matrix <- data.frame(x = rnorm(100, mean=mean1, sd=sd1), y= rnorm(100, mean=mean2, sd=sd2), z =rnorm(100, mean=mean3, sd=sd3), p = rnorm(100, mean = mean1, sd=sd1), q = rnorm(100, mean=mean2, sd=sd2), r = rnorm(100,mean=20, sd=sd3), i = rnorm(100, mean=mean1, sd=sd1), j= rnorm(100, mean=mean2, sd=sd2), k =rnorm(100, mean=20, sd=sd3))

## just looking at this data, you should see that the means of x and p are very close, means of y and q are close, and means of z and r are close. In theory, we should see ~3 clusters.
## just a quick sanity check to see what this data looks like. I suggest running one line at a time to see the graphs pile up.
hist(fake.matrix$x, col='gray', xlim=c(-5, 25), ylim=c(0,50), alpha=0.5)
hist(fake.matrix$y, col='darkorange', add=TRUE, alpha=0.5)
hist(fake.matrix$z, col='lightblue',add=T, alpha=0.5)
hist(fake.matrix$p, col='honeydew4',add=T, alpha=0.5)
hist(fake.matrix$q, col='darkorange3',add=T, alpha=0.5)
hist(fake.matrix$r, col='steelblue1',add=T, alpha=0.5)
hist(fake.matrix$i, col='azure2', add=T, alpha=0.5)
hist(fake.matrix$j, col='coral1', add=TRUE, alpha=0.5)
hist(fake.matrix$k, col='turquoise2',add=T, alpha=0.5)
#add legend
legend('topright', c('x', 'y', 'z','p','q','r','i','j','k'), fill=c('gray', 'darkorange','lightblue','honeydew4','darkorange3','steelblue1','azure2', 'coral1','turquoise2'), title='samples')

## how many clusters for this fake dataset? Arbitrary starting point. 3 because did 3 different rnorms.
my.clusters <- 3
## assign your data
## data has to be transposed? depends. see here: https://www.biostars.org/p/406259/
## The function kmeans clusters the ROWS. So if you need to tanspose your dataset, can do:
## This transposed step will cluster the "samples", x, y, z, p, q, r. BUT you must change the max number of clusters, otherwise the script will poop out when you have # clusters > # samples. 
k.means.dat <- t(fake.matrix)
# k.means.dat <- fake.matrix
## make new df for kmeans only, if don't have matrix of log.fold.changes/averages. Don't want k-means on the replicates...too messy. Do "profiles" or averages of replicates.
# k.means.dat <- t(fake.matrix)

#### Make a ratio(between.ss/total.ss) vs. # of clutsers graph ####
## this function takes your number of clusters set and runs k-means clustering for 10,000 iterations. It returns the ratio of (k$betweens/k$totss), as we mentioned above. This value will NOT be the same each time you run the clustering algorithm because it picks randomly 3 initial points.

k.means.opt <- function(my.clusters) {
  k <- kmeans(k.means.dat, centers=my.clusters, iter.max = 10000)
  return(k$betweenss/k$totss)
}
## how long does it take?
t0 <- Sys.time()
k.means.opt(my.clusters)
t <- Sys.time()
t-t0

## now that we have this function that returns the between.ss/tot.ss ratio, we can replicate that 1000 times for each cluster number and plot. As we increase the number of clusters, of course more variability will be explained. However, it's a log drop off, so you do want to qualitatively pick which number of clusters where, beyond that cluster number, the small increment of additional (between.ss) is not worth adding another cluster. 
## set the cluster maximum, probably no more than 20... Consider what you are clustering! 20,000 genes? or 6 or 18 samples? have to change it such that # clusters is not > # things to cluster.
cluster.max <- 9
## make empty list to put the data
k.means.opt.list <- vector("list", length=cluster.max) 
## replicate function where, for each cluster from 1:cluster.max, replicates 1000 times and returns the (between.ss/tot.ss) ratio. Each of the [[elements]] of the list is 1000 ratios for each number of clusters. Also I want to time this step.

## this takes 1 second on your desktop.
t0 <- Sys.time()
for (i in 1:cluster.max) {
  k.means.opt.list[[i]] <- replicate(1000, k.means.opt(i))
}
## now make a graph of all these replicates
## assign your x axis, number of clusters
clusters.x <- c(1:cluster.max)
## unlist your k.means.opt.list object to extract Ratio.BSS.TSS. Take the maximum Ratio.BSS.TSS of each 1000 kmeans replicates for each clutser, i (in this case, 1000 kmeans replicates for each cluster from 1 cluster (should be 0%) to 20 clusters (should be close to 100%)). If the number of clusters = number of data points, the "Ratio" will be undefined, which will show on the graph as 0% even though it's undefined.
k.means.data <- data.frame(Number.of.clusters= clusters.x, Ratio.BSS.TSS=unlist(lapply(k.means.opt.list, max)))
## now plot these two things
g.kmeans.opt <- ggplot(data = k.means.data, aes(x = Number.of.clusters, y = Ratio.BSS.TSS)) +
  geom_point() +
  theme_light()
g.kmeans.opt
t1 <- Sys.time()
t1-t0


#### Make a quantitative-based decision about number of clusters ####
k.centers <- 3
##looks so much better with outliers removed. Much more "loggy"
set.seed(22)
k.lumos <- kmeans(k.means.dat, centers=k.centers, iter.max = 10000)
k.lumos$betweenss/k.lumos$totss
## 93.6%! Seems crazy. Feel free to play around with the fake.matrix to see how changing the data will show a drop in the clustering. Again, the samples are clustering into 3 clear groups. Be careful to ask if you are clustering genes or samples. Kmeans is almost never this clean. 