#Load libraries
library("Biostrings")

#Read in fasta and annotation data
uboxes = readAAStringSet("Arabidopsis_Ubox.fasta")
annotations = read.delim("Arabidopsis_Ubox.tab", as.is=T)

#Remove U-box domain
for (i in 1:length(uboxes)) {
  myelement = uboxes[[i]]
  ubox_coords = annotations[i,"Domain..FT."]
  s1 = unlist(strsplit(ubox_coords, split='DOMAIN '))[2]
  s1 = unlist(strsplit(s1, split=' U-box'))[1]
  real_coords = as.numeric(unlist(strsplit(s1, split=' ')))
  uboxes[[i]] = myelement[-(real_coords[1]:real_coords[2])]
}

#Write to file
writeXStringSet(uboxes, file="ArabidopsisUbox_SRdomains.fasta")
