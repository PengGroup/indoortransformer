#' @title Predict indoor transformation products of organophosphorus compounds (OPCs).
#'
#' @param smilesList Input character vector of organophosphorus compounds. Must be in SMILES (Simplified Molecular Input Line Entry Specification) format.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr distinct mutate bind_rows select rowwise left_join row_number across
#' @importFrom stats setNames
#'
#' @return A tibble containing precursor SMILES, product SMILES, list of sequential transformations, product ID#, exact product mass, calculated logP values, and predicted positive/negative fragment m/z for organophosphate ester products.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c("CO[P](=S)(OC)Oc1ccc(cc1)[S](=O)(=O)N(C)C","CO[P](=O)(OC)OC=C(Cl)Cl")
#' results <- trans.Products(x)
#' }

trans.Products<-function(smilesList){
  prodSMILES <- prodID <- precID <- Reactions <- mz_pos <- mz_neg <- NULL #Global variables declaration
  funGroups<-c("ester_PO", "P=S", "PR3", "C=C")
  rxns<-c("Oxi", "Hyd", "HOCl", "Ozo")

  prodTable_output<-tibble(precSMILES=NULL, Reactions=NULL, prodSMILES=NULL, skipCheck=NULL) #Creates final output table to store predicted transformation product information. This table will be returned at the end of this function.
  #precSMILES = SMILES of the input compound, Reactions = character vector representing the sequence of reaction types resulting in this product, prodSMILES = product structure in SMILES format

  for(i in 1:length(smilesList)){

    prods_temp<-tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL) #Creates final output table to store predicted transformation product information. This table will be returned at the end of this function.
    productRaw<-tibble(precSMILES = smilesList[i], Reactions = list("Raw"), prodSMILES = smilesList[i])

    rxnSave<-NULL
    savePoint<-FALSE
    skipCheck<-NULL
    productsNew<-recursive.Predict(smilesList[i], prods_temp, funGroups, rxns, rxnSave, savePoint, skipCheck)

    if(!is.null(productsNew)){

      productsNew<-productsNew %>%
        distinct(prodSMILES, .keep_all = TRUE) %>%
        mutate(precSMILES = smilesList[i])

      print(paste(i,
                  "out of",
                  length(smilesList),
                  "seed compounds complete.  There were",
                  nrow(productsNew),
                  "products predicted for this compound.", sep = " "))

      prodTable_output<-prodTable_output %>%
        bind_rows(productRaw, productsNew) %>%
        select(-c(savePoint, skipCheck))
    }
  }

  if(nrow(prodTable_output) < 1){
    print("There were no products predicted for the input compounds.")
  }else{
    prodTable_output <- prodTable_output %>%
      mutate(prodID = row_number()) %>%
      rowwise() #%>%
#      mutate(MW = rcdk::get.exact.mass(rcdk::parse.smiles(prodSMILES)[[1]]), xlogP = rcdk::get.xlogp(rcdk::parse.smiles(prodSMILES)[[1]]))
  }

  fragments<-OPCproducts.Frag(prodTable_output$prodSMILES)
  prodTable_output <- prodTable_output %>%
    left_join(fragments, by = "prodSMILES") %>%
    rowwise() %>%
    mutate(across(c(Reactions, mz_pos, mz_neg), ~paste(unlist(.), collapse = ','))) #%>%
#    setNames(c("Precursor", "Reactions", "Product", "Product_ID", "MW", "xlogP", "Pos_Frag_mz", "Neg_Frag_mz"))

  return(prodTable_output)
}
