library(stringr)
combocsv <- read.csv("temp_combo_reads.csv",sep=",",header=FALSE);
probs_for_locs <- read.csv("temp_cpg_locs_probs_from_R.csv",header=TRUE);
colnames(probs_for_locs) <- c("V1","V2","V3")
lapply(1:nrow(combocsv), function(x) {
  gregexpr("CG",combocsv[x,2])[[1]] -> matches
  base <- strsplit(strsplit(as.character(combocsv[x,1]),':')[[1]][2], '-')[[1]][1]
  if (length(matches) == 1){
    if (matches == -1) return(NA)
  }
  (as.numeric(base) + as.numeric(matches) - 1)
}) -> cpglocs_in_read;

lapply(cpglocs_in_read, function(x) {
  if (length(x)==1) { if(is.na(x))return(NA)}
  unlist(lapply(x, function(y){ probs_for_locs[probs_for_locs[,2]==y,3] }))
}) -> probs_for_cpgs_in_read;

lapply(probs_for_cpgs_in_read, function(x) {
  if (length(x) == 1) { if (is.na(x)) return(NA)}
  runif(n=length(x),0,1) -> flip;
  flip < x
}) -> Methylate

unlist(lapply(1:nrow(combocsv), function(x) {
  base <- strsplit(strsplit(as.character(combocsv[x,1]),':')[[1]][2], '-')[[1]][1]
  cpglocs_in_read[[x]] -> cpglocs; 
 # if (length(cpglocs) == 1) {if (is.na(cpglocs)) return(combocsv[x,2])}
  methylatecur <- Methylate[[x]]; 
  cpglocs[methylatecur] -> toChange;
#  if (!any(methylatecur)) {return(combocsv[x,2])}
  stridx <- toChange - as.numeric(base)+1; 
  rdstring <- as.character(combocsv[x,2])
  for (idx in stridx){
  substr(rdstring,idx,idx) <- 'T'
  }
  return(rdstring)
})) -> MethylatedReads

write.csv(as.data.frame(cbind(combocsv[,1],MethylatedReads)),file="Methylated_reads_out_of_R.csv")
