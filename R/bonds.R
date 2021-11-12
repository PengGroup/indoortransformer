#' @title Quickly detects specific functional groups in a molecule
#'
#' @description Detects specific functional groups in a molecule (SMILES format) and returns a vector of functional group strings named with atom indices (Structure-Data Format). Faster and more lightweight than findBonds().
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



#' @title Detects specific functional groups in a molecule.
#'
#' @description Detects specific functional groups in a molecule (SMILES format) and returns a table containing detailed SDF information. Slower than smartBonds(), but provides more detailed information.
#'
#'
#' @param smiles A SMILES-formatted molecule string.
#' @param targetGroup A predefined character vector of functional groups.
#'
#' @return A tibble containing bond block and atom block information for each functional group detected.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_detect str_split str_extract
#' @importFrom tibble tibble as_tibble rownames_to_column
#' @import dplyr
findBonds<-function(smiles, targetGroup) { #Function to find atom and bond row IDs containing desired functional group

  sdfList<-ChemmineR::smiles2sdf(smiles)
  atomRows<-names(targetGroup)[1] %>%
    str_split(pattern = " ") %>%
    unlist() %>%
    as.integer()

  atomTable<-ChemmineR::atomblock(sdfList[[1]])
  atomTibble<-as.data.frame(atomTable) %>%  #"as.data.frame()" is required because the atomtable list items are stored as arrays
    rownames_to_column(var="atoms") %>% #Extracts the atom names from the SDF Atoms block table row names and puts them into their own column (presently in "C_1" format)
    as_tibble() %>%
    transmute(atoms = str_extract(.data$atoms,"[:alpha:]"), atomID = c(1:nrow(.data))) %>% #"[:alpha:]" is a REGEX command which extracts only letters from a given string, in this case atom names ("C", "P", etc.)
    filter(.data$atomID %in% atomRows)

  bondTibble<-ChemmineR::bondblock(sdfList[[1]]) %>%
    as.data.frame() %>% #"as.data.frame()" is required because the bondtable list items are stored as arrays
    as_tibble() %>%
    select(c(1:3)) %>% #Extract only the first 3 columns from the SDF Bonds block table.
    stats::setNames(c("atom1", "atom2", "bondType")) %>% #These are what the values in the first 3 columns correspond to: the two atoms involved in the bond, and the type (1 = single, 2 = double, 3 = triple)
    mutate(bondID = c(1:nrow(.data))) #Stores the SDF bond number as its own column. This is important because certain rows will be filtered out below, but we want to be able to easily find important bonds during the transformation prediction stage. So the bond ID# is stored and returned as its own column in the returned variable "bondTable_result".
  if(nrow(atomTibble) > 1){
    bondTibble<-bondTibble %>%
      filter(.data$atom1 %in% atomTibble$atomID) %>%
      filter(.data$atom2 %in% atomTibble$atomID)
  }else{
    bondTibble<-filter(bondTibble, .data$atom1 == atomTibble$atomID | .data$atom2 == atomTibble$atomID)
  }


  if(targetGroup == "ester_PO"){
    bondTibble<-filter(bondTibble, .data$bondType == 1)
    atomTibble<-filter(atomTibble, .data$atomID %in% bondTibble$atom1 | .data$atomID %in% bondTibble$atom2)
    bondTable_result<-tibble(SMILES = smiles, funGroup = targetGroup, bondRow = bondTibble$bondID, Patom = atomTibble$atomID[str_detect(atomTibble$atoms, "P")], Oatom = atomTibble$atomID[str_detect(atomTibble$atoms, "O")], Catom = NA, Satom = NA)
    #CmpdID = ID # of the input compound, SMILES = input compound structure in SMILES format, funGroup = string representing the functional group described by the remaining columns, bondRow = SDF Bonds block row index containing the atoms which are part of this functional group (stored as a list of values), Patom/Oatom/Catom/Satom = SDF Atoms block row index for the phosphorus/oxygen/carbon/sulphur atoms involved in the bond (stored as a list of values)
  }else if(targetGroup == "PR3"){
    bondTable_result<-tibble(SMILES = smiles, funGroup = targetGroup, bondRow = bondTibble$bondID, Patom = atomTibble$atomID[str_detect(atomTibble$atoms, "P")], Oatom = NA, Catom = NA, Satom = NA)
  }else if(targetGroup == "P=S"){
    bondTable_result<-tibble(SMILES = smiles, funGroup = targetGroup, bondRow = bondTibble$bondID, Patom = atomTibble$atomID[str_detect(atomTibble$atoms, "P")], Oatom = NA, Catom = NA, Satom = atomTibble$atomID[str_detect(atomTibble$atoms, "S")])
  }else if(targetGroup == "C=C"){
    bondTable_result<-tibble(SMILES = smiles, funGroup = targetGroup, bondRow = bondTibble$bondID, Patom = NA, Oatom = NA, Catom = list(atomTibble$atomID), Satom = NA)
    #For alkene groups, both atom IDs must be stored.
  }else{bondTable_result<-NULL}

  return(bondTable_result)
}
