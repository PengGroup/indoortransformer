#' Identify desired functional groups in a molecule (in SMILES format) and return a vector of functional group strings named with atom indices (Structure-Data Format)
#'
#' @param smiles A SMILES-formatted molecule string.
#' @param grps A predefined character vector of functional groups.
#'
#' @return A named character vector.
#' @export
#'
#' @examples
#' x <- "CCC=CCCO[P](OCCC)(=O)OCCCC"
#' y <- c("C=C", "PR3", "P=S", "ester_PO")
#' smartBonds(x, y)
smartBonds<-function(smiles, grps){

  if(length(smiles) < 1){return(NULL)}
  if(is.na(smiles)){return(NULL)}

  smart_grps<-dplyr::case_when(
    grps == "ester_PO" ~ "[P]([OX2H0])(=O)",
    grps == "P=S" ~ "[P]=[S]",
    grps == "PR3" ~ "[PX3]",
    grps == "C=C" ~ "[C]=[C]",
    TRUE ~ as.character(grps)
  )
  molecule<-rcdk::parse.smiles(smiles)
  sites<-NULL
  for(i in 1:length(smart_grps)){
    sites_temp<-rcdk::matches(smart_grps[i], molecule, return.matches = TRUE)
    numSites<-length(sites_temp[[1]]$mapping)
    sitesAdd<-NULL
    for(j in 1:numSites){
      sitesAdd[j]<-grps[i]
      atomIndices<-sites_temp[[1]]$mapping[[j]]
      atomIndices<-atomIndices + 1 #Required because the atom mapping in rdck is zero-indexed, but SDF atom indices start at 1
      names(sitesAdd)[j]<-stringr::str_flatten(atomIndices, collapse = " ") #Save the atom indices for the functional group as vector names
    }
    if(numSites > 0){sites<-c(sites, sitesAdd)}else{next}
  }
  return(sites)
}
