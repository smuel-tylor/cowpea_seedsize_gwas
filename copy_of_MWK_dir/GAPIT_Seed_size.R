#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")
library(GAPIT3)

# source("http://zzlab.net/GAPIT/GAPIT.library.R")
# source("https://zzlab.net/GAPIT/previous/gapit_functions20220122.txt")
# 
# source("http://zzlab.net/GAPIT/GAPIT.library.R")
# source("http://zzlab.net/GAPIT/gapit_functions.txt")

# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)



# import dataGAPIT3
myG <- read.table("Genotype_seed_size.hmp.txt", header = F, sep = "\t")
myY <- read.table("phenotype_Seed_size.txt", header = T, sep = "\t")


myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  model=c("Blink"),# choose model
  PCA.total=3,
  Random.model=FALSE,
  #NJtree.group=4,                                       # set the number of clusting group in Njtree plot
  #QTN.position=mysimulation$QTN.position,
  #Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  #PCA.3d=TRUE,                                          # plot 3d interactive PCA
  Phenotype.View= FALSE, 
  file.output=T
)
