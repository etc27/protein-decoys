#Load libraries
library("Biostrings")

#Read in fasta and annotation data
fboxes = readAAStringSet("yeast_allfbox.fasta")
annotations = read.delim("yeast_allfbox.tab", as.is=T)

#Remove F-box domain
for (i in 1:length(fboxes)) {
  myelement = fboxes[[i]]
  fbox_coords = annotations[i,"Domain..FT."]
  s1 = unlist(strsplit(fbox_coords, split='DOMAIN '))[2]
  s1 = unlist(strsplit(s1, split=' F-box'))[1]
  real_coords = as.numeric(unlist(strsplit(s1, split=' ')))
  fboxes[[i]] = myelement[(real_coords[2]+1):length(myelement)]
}

#Write to file
writeXStringSet(fboxes, file="FBXproteins_SRdomains.fasta")

