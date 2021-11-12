#' Construct a new SDFset object for a single molecule from modified atom block and bond block tables (??"SDFset-class" for details).
#'
#' @param prodTables A list of two matrices: an atoms table and a bonds table
#' @param sdf The unmodified SDF object from which the original atom block and bond block tables were obtained
#'
#' @return A new SDFset object containing one molecule.
#' @export
makeSDF<-function(prodTables, sdf){

  atomTable_temp<-prodTables[[1]]
  bondTable_temp<-prodTables[[2]]
  #  countsProd<-c(nrow(atomTable_temp),nrow(bondTable_temp),"0  0  0  0  0  0  0  0999 V2000")

  atomNum<-as.character(nrow(atomTable_temp))
  atomNum_length<-stringr::str_length(atomNum)
  atomNum<-paste0(paste0(rep(" ", 3-atomNum_length), collapse=""), atomNum, collapse="")

  bondNum<-as.character(nrow(bondTable_temp))
  bondNum_length<-stringr::str_length(bondNum)
  bondNum<-paste0(paste0(rep(" ", 3-bondNum_length), collapse=""), bondNum, collapse="")

  countsLine<-ChemmineR::header(sdf[[1]])[[4]]
  stringr::str_sub(countsLine, 1, 6)<-paste0(atomNum, bondNum, collapse="")

  prodSDF<-sdf[[1]]
  prodSDF@atomblock<-atomTable_temp
  prodSDF@bondblock<-bondTable_temp
  prodSDF@header[[4]]<-countsLine
  prodSDF<-ChemmineR::SDFset(list(prodSDF)) #Needs to be returned as "SDFset" because that is the only format accepted by ChemmineR::sdf2smiles().
  prodSDF<-ChemmineR::canonicalize(prodSDF)

  return(prodSDF)
}
