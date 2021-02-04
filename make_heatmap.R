rm(list=ls())
setwd("/Users/emmacorcoran/Documents/Rotations/Gendron-Rotation/Bioinformatics/miP3-master/exampleData")

#Load libraries
library("Biostrings")
library(d3heatmap) #makes interactive heatmap!
library(shiny)
library(scales)

#Load fasta files
targets = readAAStringSet("CelegansFBX_SRdomains.fasta")
target_names = names(targets) #names of fbox proteins
for (i in 1:length(target_names)) {
  target_names[i] = unlist(strsplit(target_names[i], split=' OS='))[1]
}
setwd("/Users/emmacorcoran/Documents/Rotations/Gendron-Rotation/Bioinformatics/miP3-master/output_results")
fbox_decoys = readAAStringSet("CelegansFbox_1e-20_all.fasta")
decoy_names = names(fbox_decoys) #names of decoy proteins
for (i in 1:length(decoy_names)) {
  decoy_names[i] = unlist(strsplit(decoy_names[i], split=' OS='))[1]
}
#Load eval file
setwd("/Users/emmacorcoran/Documents/Rotations/Gendron-Rotation/Bioinformatics/miP3-master/heatmap_results/C.elegans")
evals = data.matrix(read.csv("CelegansFbox1e-20_evals.csv", header=F))
evals = t(evals)
colnames(evals) = target_names
rownames(evals) = decoy_names
toremove = grep("FBX", decoy_names)
toremove = c(toremove, grep("F-box", decoy_names))
new_evals = evals[-toremove,]
evals = new_evals

##Analysis of decoys
number_targets = matrix(list(), nrow=nrow(evals), ncol=3)
#sort the decoys into different categories
for (i in 1:nrow(evals)) {
  nonNa_vals = evals[i,which((evals[i,]!="NA")&(evals[i,]>=0))]
  one_hits = nonNa_vals[which(nonNa_vals<=1e-20)]
  number_targets[i,1] = rownames(evals)[i]
  number_targets[i,2] = length(one_hits)
  number_targets[i,3] = list(names(one_hits))
}

##Analysis of targets
number_decoys = matrix(list(), nrow=ncol(evals), ncol=3)
nodecoys = c()
toremove = c()
#make a matrix for each fbox of number of decoys and name of decoys
for (i in 1:ncol(evals)) {
  nonNa_vals = evals[which((evals[,i]!="NA")&(evals[,i]>=0)),i]
  significant_vals = nonNa_vals[which(nonNa_vals<=1e-20)]
  number_decoys[i,1] = colnames(evals)[i]
  number_decoys[i,2] = length(significant_vals)
  if (length(significant_vals) == 0) {
    nodecoys = c(nodecoys, colnames(evals)[i])
    toremove = c(toremove, i)
    number_decoys[i,3] = "NA"
  }
  else {
    number_decoys[i,3] = list(names(significant_vals))
  }
}
write.table(nodecoys, file="CelegansFboxall_nodecoys.csv", quote=F, sep=",", row.names=F, col.names=F)

##remove targets that don't have a decoy then sort by number of targets
evals_sorted = evals[,-toremove] #remove fboxes with no decoys
colnames(number_decoys) = c("Fbox protein", "Number decoys", "Name of decoys")
number_decoys = number_decoys[-toremove,] #remove fboxes with no decoys

####sort fboxes by number of decoys and shared decoys
evals_sorted = evals_sorted[order(unlist(number_targets[,2])) ,order(unlist(number_decoys[,2]))]
number_decoys = number_decoys[order(unlist(number_decoys[,2])), ]
number_targets = number_targets[order(unlist(number_targets[,2])), ]

newdata = matrix(list(), nrow=0, ncol=3) #rows are decoys and columns are fbox
colnames(newdata) = c("Fbox protein", "Number decoys", "Name of decoys")
indexorder = c()
for (i in 1:max(unlist(number_decoys[,2]))) {
  currdata = matrix(number_decoys[which(unlist(number_decoys[,2])==i),],ncol=3) #all fboxes that have i number of decoys
  ###sort by number of intersections
  j=1
  while (j <= nrow(currdata)) {
    first_targets = unlist(currdata[j,3])
    newdata = rbind(newdata, currdata[j,])
    indexorder = c(indexorder,  which(colnames(evals_sorted) == currdata[j,1]))
    currdata = matrix(currdata[-j,], ncol=3)
    k=1
    while (k <= (nrow(currdata)-1)) {
      second_targets = unlist(currdata[k+1,3])
      if (length(intersect(first_targets, second_targets)) >= i/2) {
        newdata = rbind(newdata, currdata[(k+1), ])
        indexorder = c(indexorder, which(colnames(evals_sorted) == currdata[k+1,1]))
        currdata = matrix(currdata[-(k+1),], ncol=3)
      }
      else {
        k = k+1
      }
    }
  }
}
write.table(newdata[,1:2], file="sortedFBox_1e-20.csv", quote=F, sep=",", row.names=F) ##make table of newdata (sorted fboxes with their number of decoys)

#####sort decoys
decoydata = matrix(list(), nrow=0, ncol=3) #rows are decoys and columns are fbox
colnames(decoydata) = c("Decoy", "Number targets", "Name of targets")
indexorder1 = c()
for (i in 1:max(unlist(number_targets[,2]))) {
  currdata = matrix(number_targets[which(unlist(number_targets[,2])==i),],ncol=3) #all decoys that have i number of targets
  ###sort by number of intersections
  j=1
  while (j <= nrow(currdata)) {
    first_targets = unlist(currdata[j,3])
    decoydata = rbind(decoydata, currdata[j,])
    indexorder1 = c(indexorder1,  which(rownames(evals_sorted) == currdata[j,1]))
    currdata = matrix(currdata[-j,], ncol=3)
    k=1
    while (k <= (nrow(currdata)-1)) {
      second_targets = unlist(currdata[k+1,3])
      if (length(intersect(first_targets, second_targets)) >= i/2) {
        decoydata = rbind(decoydata, currdata[(k+1), ])
        indexorder1 = c(indexorder1, which(rownames(evals_sorted) == currdata[k+1,1]))
        currdata = matrix(currdata[-(k+1),], ncol=3)
      }
      else {
        k = k+1
      }
    }
  }
}
write.table(decoydata[,1:2], file="sorteddecoy_1e-20.csv", quote=F, sep=",", row.names=F) ##make table of newdata (sorted fboxes with their number of decoys)
evals_fboxsort = evals_sorted[indexorder1,indexorder]
write.table(evals_fboxsort, file="tograph.csv", quote=F, sep=",", row.names=T, col.names=T)

#Make heatmap
myCol <- c("red", "Dark Orange 1", "Orange 1", "white")
myBreaks <- c(0, 1e-60, 1e-40, 1e-20, 10)
colorFunc <- col_bin(myCol, bins = rescale(myBreaks))
d3heatmap(evals_fboxsort, dendrogram="none", cexRow=.1, cexCol=.1, Rowv=F, Colv=F, show_grid=T, anim_duration=0, colors=colorFunc, digits=200)
