#' @title Predict common fragments of OPCs (organophosphorus compounds).
#'
#' @param smilesList Input character vector of organophosphorus compounds. Must be in SMILES (Simplified Molecular Input Line Entry Specification) format.
#'
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr select bind_rows distinct rowwise mutate case_when
#' @importFrom rcdk get.mol2formula parse.smiles
#'
#' @return Tibble containing input SMILES and predicted positive/negative fragment m/z values for each.
OPCproducts.Frag<-function(smilesList){

  funGroups <- "ester_PO"
  rxns <- "Hyd"

  fragTable_output<-tibble(prodSMILES = NULL, mz_pos = NULL, mz_neg = NULL) #Creates final output table to store predicted fragment m/z values. This table will be returned at the end of this function.

  for(i in 1:length(smilesList)){

    frags_temp<-tibble(Reactions = NULL, fragSMILES = NULL, skipCheck = NULL)
    unfragmented<-tibble(fragSMILES = smilesList[i], fragSide = "unfrag")

    rxnSave<-NULL
    skipCheck<-NULL
    fragNew<-recursive.Frag(smilesList[i], frags_temp, funGroups, rxns, rxnSave, skipCheck)

    if(!is.null(fragNew)){

      Omass<-iso_list$mass[7]
      Hmass<-iso_list$mass[1]
      Pmass<-iso_list$mass[30]
      emass<-0.000549

      fragNew <- fragNew %>%
        select(-c(Reactions, skipCheck)) %>%
        bind_rows(unfragmented) %>%
        distinct(fragSMILES, .keep_all = TRUE) %>%
        rowwise() %>%
        mutate(fragMass = get.mol2formula(parse.smiles(fragSMILES)[[1]],charge=0)@mass) %>%
        mutate(
          mz_pos = case_when(
            fragSide == "Rside" ~ list(fragMass-Omass-Hmass-emass),
            fragSide == "Pside" ~ list(c(fragMass+Hmass-emass, fragMass-Omass-Hmass-emass)),
            fragSide == "unfrag" ~ list(fragMass+Hmass-emass),
            TRUE ~ list(NA)
          ),
          mz_neg = case_when(
            fragSide == "Rside" ~ list(c(fragMass-Hmass+emass, fragMass-(2*Hmass)+(2*Omass)+Pmass+emass)),
            fragSide == "Pside" ~ list(NA),
            fragSide == "unfrag" ~ list(NA),
            TRUE ~ list(NA)
          )
        )

      fragTable_temp<-tibble(prodSMILES = smilesList[i], mz_pos = list(unlist(fragNew$mz_pos[!is.na(fragNew$mz_pos)])), mz_neg = list(unlist(fragNew$mz_neg[!is.na(fragNew$mz_neg)])))

      print(paste0("Fragments predicted for ",
                   i,
                   " out of ",
                   length(smilesList),
                   " transformation products."))

      fragTable_output<-bind_rows(fragTable_output, fragTable_temp)
    }
  }

  if(nrow(fragTable_output) < 1){
    print("There were no products predicted for the input compounds.")
  }
  return(fragTable_output)
}



#' @title Recursive prediction of transformation products.
#'
#' @param SMI Character vector of length one containing a SMILES string.
#' @param prods Empty three-column tibble in which to store results.
#' @param funGroups Character vector containing functional groups at which transformation reactions can occur.
#' @param rxnList Character vector containing types of reactions to use for prediction.
#' @param rxnSave Empty character vector for storing sequential reaction information.
#' @param skipCheck Empty character vector for storing SMILES strings of predicted products. Used to improve efficiency by not performing the same prediction multiple times.
#'
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stats na.exclude
#' @importFrom dplyr bind_rows filter
#'
#' @return Tibble containing fragment SMILES, list of sequential fragmentation reactions, and saved skipCheck strings.
recursive.Frag<-function(SMI, prods, funGroups, rxnList, rxnSave, skipCheck){

  rxnGroups<-fragSites(SMI, funGroups)

  if(length(rxnGroups)<1){
    return(NULL)
  }else{
    prods_temp<-tibble(Reactions=NULL, fragSMILES=NULL, skipCheck=NULL) #Creates final output table to store predicted fragment information. This table will be returned at the end of this function.
    prods_Cleave2 <- tibble(Reactions=NULL, fragSMILES=NULL, skipCheck=NULL)
    prods_Cleave1 <- tibble(Reactions=NULL, fragSMILES=NULL, skipCheck=NULL)

    for(j in 1:length(rxnGroups)){

      prod_new<-predict.Frag(SMI, rxnGroups[j], rxnList, rxnSave)
      if(is.null(prod_new)){next}

      if(any(prod_new$fragSMILES %in% skipCheck)){prod_new$fragSMILES[which(prod_new$fragSMILES %in% skipCheck)]<-NA} #If any of the predicted fragments are duplicates of previously predicted fragments, remove them from 'prod_new'

      skipCheck<-skipCheck %>%
        c(prod_new$fragSMILES) %>%
        na.exclude() %>%
        unique()


      prod_new2 <- prod_new[2,] #Second cleavage product
      prod_new1 <- prod_new[1,] #First cleavage product

      rxnSave_current<-prod_new2$Reactions

      prods_Cleave2 <- tibble(Reactions=NULL, fragSMILES=NULL, skipCheck=NULL)
      prods_Cleave2 <- bind_rows(prods_Cleave2, prod_new2, recursive.Frag(prod_new2$fragSMILES[1], prods, funGroups, rxnList, rxnSave_current, skipCheck))
      prods_Cleave2 <- filter(prods_Cleave2, !is.na(fragSMILES))
      skipCheck<-c(skipCheck, unlist(last(prods_Cleave2$skipCheck))) #Using the last row to extract SMILES to skip because the first row(s) will be NULL
      skipCheck<-unique(skipCheck)

      rxnSave_current<-prod_new1$Reactions

      prods_Cleave1 <- tibble(Reactions=NULL, fragSMILES=NULL, skipCheck=NULL)
      prods_Cleave1 <- bind_rows(prods_Cleave1, prod_new1, recursive.Frag(prod_new1$fragSMILES[1], prods, funGroups, rxnList, rxnSave_current, skipCheck))
      prods_Cleave1 <- filter(prods_Cleave1, !is.na(fragSMILES))
      skipCheck<-c(skipCheck, unlist(last(prods_Cleave1$skipCheck)))
      skipCheck<-unique(skipCheck)


      prods_temp<-bind_rows(prods_temp, prods_Cleave1, prods_Cleave2)
      prods_temp$skipCheck<-list(skipCheck)
      prods_temp<-filter(prods_temp, !is.na(fragSMILES))
    }

    return(prods_temp)
  }
}



#' @title Quickly detects specific functional groups in a molecule
#'
#' @param smiles A SMILES-formatted molecule string.
#' @param grps A predefined character vector of functional groups.
#'
#' @importFrom rcdk parse.smiles
#' @importFrom stringr str_flatten
#'
#' @return A named character vector.
fragSites<-function(smiles, grps){

  if(length(smiles) < 1){return(NULL)}
  if(is.na(smiles)){return(NULL)}

  if("ester_PO" %in% grps){
    smart_grps<-"[P]([OX2H0])(=O)"
  }else{return(NULL)}

  molecule<-parse.smiles(smiles)
  sites<-NULL
  for(i in 1:length(smart_grps)){
    sites_temp<-rcdk::matches(smart_grps[i], molecule, return.matches = TRUE)
    numSites<-length(sites_temp[[1]]$mapping)
    sitesAdd<-NULL
    for(j in 1:numSites){
      sitesAdd[j]<-grps[i]
      atomIndices<-sites_temp[[1]]$mapping[[j]]
      atomIndices<-atomIndices + 1 #Required because the atom mapping in rdck is zero-indexed, but SDF atom indices start at 1
      names(sitesAdd)[j]<-str_flatten(atomIndices, collapse = " ") #Save the atom indices for the functional group as vector names
    }
    if(numSites > 0){sites<-c(sites, sitesAdd)}else{next}
  }
  return(sites)
}



#' @title Predict fragments from specified reaction events on a single molecule.
#'
#' @param cmpd SMILES string representing a single molecule.
#' @param grp Single-element named character vector containing a functional group. Name represents atom indices. Output of smartBonds().
#' @param rxnList Character vector of allowed reaction types.
#' @param rxnSave Character vector containing saved reactions preceding this one in sequence.
#'
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#'
#' @return Returns a one-line tibble containing: reaction sequence up to this point, product molecule (in SMILES format), and the current value of savePoint (TRUE or FALSE).
predict.Frag<-function(cmpd, grp, rxnList, rxnSave){
  rxnSite_IDs<-findBonds(cmpd, grp)
  product<-NULL
  fragSide<-NULL

  if(grp == "ester_PO"){
    if(any(str_detect(rxnList, "Hyd"))){
      product<-predict.Hydrolysis(rxnSite_IDs)
      rxnSave_temp<-"Hyd"
      fragSide <- c("Rside", "Pside")
    }
  }

  if(length(product) > 0){
    rxnSave_current<-unlist(rxnSave)
    rxnSave<-list(c(rxnSave_current, rxnSave_temp))
    productReturn<-tibble(Reactions = rxnSave, fragSMILES = product, fragSide = fragSide)
  }else{
    productReturn<-NULL
  }

  return(productReturn)
}
