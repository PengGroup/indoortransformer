#' @title Predict reaction products of specified reactions with a single molecule.
#'
#' @param cmpd SMILES string representing a single molecule.
#' @param grp Single-element named character vector containing a functional group. Name represents atom indices. Output of smartBonds().
#' @param rxnList Character vector of allowed reaction types.
#' @param rxnSave Character vector containing saved reactions preceding this one in sequence.
#' @param savePoint Single-element logical vector. Allows for different reactions at the same reaction site within recursivePredict().
#'
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#'
#' @return Returns a one-line tibble containing: reaction sequence up to this point, product molecule (in SMILES format), and the current value of savePoint (TRUE or FALSE).
predict.Product<-function(cmpd, grp, rxnList, rxnSave, savePoint){
  rxnSite_IDs<-findBonds(cmpd, grp)
  product<-NULL
  savecheck<-FALSE

  if(grp == "P=S" | grp == "PR3"){
    if(any(str_detect(rxnList, "Oxi"))){
      product<-predict.Oxi(rxnSite_IDs)
      rxnSave_temp<-"Oxi"
    }
  }

  if(grp == "ester_PO"){
    if(any(str_detect(rxnList, "Hyd"))){
      product<-predict.Hydrolysis(rxnSite_IDs)
      rxnSave_temp<-"Hyd"
    }
  }

  if(grp == "C=C"){
    if(savePoint == TRUE){
      if(any(str_detect(rxnList, "Ozo"))){
        product_Ozonolysis<-predict.Ozonolysis(rxnSite_IDs)
        product_SecOzon<-predict.SecOzon(rxnSite_IDs)
        product<-c(product_Ozonolysis, product_SecOzon)
        rxnSave_temp<-"Ozo"
      }
    }
    if(savePoint == FALSE){
      if(any(str_detect(rxnList, "HOCl"))){
        product<-predict.HOClAdd(rxnSite_IDs)
        rxnSave_temp<-"HOCl"
        savecheck<-TRUE
      }
    }
  }


  if(length(product) > 0){
    rxnSave_current<-unlist(rxnSave)
    rxnSave<-list(c(rxnSave_current, rxnSave_temp))
    productReturn<-tibble(Reactions = rxnSave, prodSMILES = product, savePoint = savecheck)
  }else{
    productReturn<-NULL
  }

  return(productReturn)
}



#' Predicts the product of oxidation of a phosphine or thiophosphate.
#'
#' @param rxnSite One-row tibble containing: molecule (as SMILES string), reaction site type (string), SDF bondblock row number, and SDF atomblock row numbers for phosphorus, oxygen, sulfur, and/or carbon. See findBonds() for more details.
#'
#' @importFrom magrittr %>%
#' @importFrom ChemmineR smiles2sdf atomblock bondblock sdf2smiles
#' @importFrom tidyr unnest
#' @importFrom dplyr distinct group_by filter ungroup row_number
#' @importFrom stringr str_replace
#'
#' @return A character vector of length one containing a SMILES string.
predict.Oxi<-function(rxnSite){
  bondRow <- Patom <- Oatom <- Satom <- funGroup <- NULL

  currentSDF<-smiles2sdf(rxnSite$SMILES[1])
  rxnSite<-rxnSite %>% #Extract the rows corresponding to the current molecule and assign reaction types based on functional group
    unnest(c(bondRow, Patom, Oatom, Satom)) %>%
    distinct() %>%
    group_by(funGroup) %>%
    filter(row_number()==1 | funGroup != "PR3") %>% #Added to prevent the same PR3 group from being oxidized 3 times. May need to remove this line if other reactions are added which apply to PR3 groups.
    ungroup()

  if(rxnSite$funGroup[1] == "P=S"){
    atomTable_prod<-atomblock(currentSDF[[1]])
    bondTable_prod<-bondblock(currentSDF[[1]])

    atomReplace<-rxnSite$Satom
    atomID_prod<-rownames(atomTable_prod)[atomReplace]
    rownames(atomTable_prod)[atomReplace]<-str_replace(atomID_prod, "S", "O")

    prodTable_list<-list(atomTable_prod, bondTable_prod)
    currentSDF<-makeSDF(prodTable_list, currentSDF)
    prodSMI<-sdf2smiles(currentSDF)
    return(prodSMI@smilist[[1]])
  }

  if(rxnSite$funGroup[1] == "PR3"){
    atomTable_prod<-atomblock(currentSDF[[1]])
    bondTable_prod<-bondblock(currentSDF[[1]])

    rxnSite_temp<-filter(rxnSite, funGroup=="PR3")

    atomTable_prod<-rbind(atomTable_prod, c(rep(0, ncol(atomTable_prod))))
    newAtom_ID<-nrow(atomTable_prod)
    rownames(atomTable_prod)[newAtom_ID]<-paste0("O_", newAtom_ID)

    atomSite<-rxnSite_temp$Patom
    bondTable_prod<-rbind(bondTable_prod, c(atomSite, newAtom_ID, 2, rep(0, ncol(bondTable_prod)-3)))
    rownames(bondTable_prod)[nrow(bondTable_prod)]<-as.character(nrow(bondTable_prod))

    prodTable_list<-list(atomTable_prod, bondTable_prod)
    currentSDF<-makeSDF(prodTable_list, currentSDF) #Saves the product SDF for use by reaction blocks later in the sequence
    prodSMI<-sdf2smiles(currentSDF)
    return(prodSMI@smilist[[1]])
  }
}


#' @title Predicts the products of hydrolysis at a phosphate ester linkage.
#'
#' @param rxnSite One-row tibble containing: molecule (as SMILES string), reaction site type (string), SDF bondblock row number, and SDF atomblock row numbers for phosphorus, oxygen, sulfur, and/or carbon. See findBonds() for more details.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom ChemmineR smiles2sdf atomblock bondblock sdf2smiles
#' @importFrom tidyr unnest
#' @importFrom dplyr distinct filter bind_rows mutate if_else rowwise select lag n last
#' @importFrom tibble rownames_to_column rowid_to_column column_to_rownames
#' @importFrom stringr str_extract
#'
#' @return A character vector of length two containing SMILES strings, one for each hydrolysis product.
predict.Hydrolysis<-function(rxnSite){

  currentSDF<-smiles2sdf(rxnSite$SMILES[1])
  rxnSite<-rxnSite %>% #Extract the rows corresponding to the current molecule and assign reaction types based on functional group
    unnest(c(.data$bondRow, .data$Patom, .data$Oatom, .data$Catom, .data$Satom)) %>%
    distinct()

  atoms<-as.data.frame(atomblock(currentSDF[[1]]))
  bonds<-as.data.frame(bondblock(currentSDF[[1]]))

  check <- TRUE
  atomTest <- rxnSite$Oatom
  atomNeg <- rxnSite$Patom

  cleaveProd1_bonds<-NULL
  while(check == TRUE){
    keepBonds<-bonds %>%
      filter(.data$C1 %in% atomTest | .data$C2 %in% atomTest) %>%
      filter(!(.data$C1 %in% atomNeg | .data$C2 %in% atomNeg))

    if(is.null(cleaveProd1_bonds)==FALSE){
      keepBonds<-keepBonds %>%
        filter(!(.data$C1 %in% c(rxnSite$Oatom, rxnSite$Patom) & .data$C2 %in% c(rxnSite$Oatom, rxnSite$Patom))) #Do not store the ester bond itself. (Required for alkenes which are inside ring groups. Otherwise the 'while loop' will continue indefinitely.)
    }

    if(nrow(keepBonds) < 1){
      check == FALSE
      break
    }
    atomNeg <- atomTest
    atomCheck<-keepBonds[,1:2]

    atomTest_temp<-unlist(atomCheck)
    atomTest<-atomTest_temp[-which(unlist(keepBonds[,1:2]) %in% atomTest)]

    cleaveProd1_bonds<-bind_rows(cleaveProd1_bonds, keepBonds)
  }

  if(is.null(cleaveProd1_bonds) == FALSE){

    newAtom_empty<-setNames(c(rep(0, ncol(atoms))), colnames(atoms))
    cleaveProd1_atoms<-atoms %>%
      rownames_to_column() %>%
      rowid_to_column() %>%
      filter(.data$rowid %in% cleaveProd1_bonds$C1 | .data$rowid %in% cleaveProd1_bonds$C2)

    if(length(which(unlist(cleaveProd1_bonds[,c(1:2)]) %in% rxnSite$Patom)) == 3){ #Accounts for the possibility of cyclic phosphate compounds. Without this check, the product might end up with a P-H group.
      cleaveProd1_atoms<-cleaveProd1_atoms %>%
        bind_rows(newAtom_empty) %>%
        mutate(rowid_old = if_else(is.na(.data$rowid), as.integer(lag(.data$rowid)+1), .data$rowid))  %>%
        mutate(rowid = seq(from=1, to=n(), by=1)) %>%
        mutate(rowname = if_else(is.na(.data$rowname), paste0("O_", .data$rowid), .data$rowname)) %>%
        mutate(rowname = paste0(str_extract(.data$rowname, "[:alpha:]"), "_", .data$rowid))
    }else{
      cleaveProd1_atoms<-cleaveProd1_atoms %>%
        mutate(rowid_old = if_else(is.na(.data$rowid), as.integer(lag(.data$rowid)+1), .data$rowid))  %>%
        mutate(rowid = seq(from=1, to=n(), by=1)) %>%
        mutate(rowname = paste0(str_extract(.data$rowname, "[:alpha:]"), "_", .data$rowid))
    }


    if(length(which(unlist(cleaveProd1_bonds[,c(1:2)]) %in% rxnSite$Patom)) == 3){ #Accounts for the possibility of cyclic phosphate compounds. Without this check, the product might end up with a P-H group.
      newBond_empty<-setNames(c(rep(0, ncol(bonds))), colnames(bonds))
      cleaveProd1_bonds<-cleaveProd1_bonds %>%
        rowid_to_column %>%
        rowwise() %>%
        mutate(C1 = if_else(.data$C1 %in% cleaveProd1_atoms$rowid_old, cleaveProd1_atoms$rowid[which(.data$C1 == cleaveProd1_atoms$rowid_old)], 0)) %>%
        mutate(C2 = if_else(.data$C2 %in% cleaveProd1_atoms$rowid_old, cleaveProd1_atoms$rowid[which(.data$C2 == cleaveProd1_atoms$rowid_old)], 0)) %>%
        column_to_rownames(var = "rowid") %>%
        bind_rows(newBond_empty) %>%
        mutate(
          C1 = replace(.data$C1, list = which(.data$C1 == 0), values = cleaveProd1_atoms$rowid[which(cleaveProd1_atoms$rowid_old == rxnSite$Patom)]),
          C2 = replace(.data$C2, list = which(.data$C2 == 0), values = last(cleaveProd1_atoms$rowid)),
          C3 = replace(.data$C3, list = which(.data$C3 == 0), values = 1)
        ) %>%
        as.matrix()
    }else{
      cleaveProd1_bonds<-cleaveProd1_bonds %>%
        rowid_to_column %>%
        rowwise() %>%
        mutate(C1 = if_else(.data$C1 %in% cleaveProd1_atoms$rowid_old, cleaveProd1_atoms$rowid[which(.data$C1 == cleaveProd1_atoms$rowid_old)], 0)) %>%
        mutate(C2 = if_else(.data$C2 %in% cleaveProd1_atoms$rowid_old, cleaveProd1_atoms$rowid[which(.data$C2 == cleaveProd1_atoms$rowid_old)], 0)) %>%
        column_to_rownames(var = "rowid") %>%
        as.matrix()
    }

    cleaveProd1_atoms<-cleaveProd1_atoms %>%
      select(-c(.data$rowid, .data$rowid_old)) %>%
      column_to_rownames(var = "rowname") %>%
      as.matrix()


    cleaveProd1<-list(cleaveProd1_atoms, cleaveProd1_bonds)
    prod1_SDF<-makeSDF(cleaveProd1, currentSDF)
    prod1_SMI<-sdf2smiles(prod1_SDF)
    prod1_SMI<-prod1_SMI@smilist[[1]]

  }else{prod1_SMI<-NULL}



  check <- TRUE
  atomTest <- rxnSite$Patom
  atomNeg <- rxnSite$Oatom
  cleaveProd2_bonds<-NULL
  while(check == TRUE){
    keepBonds<-bonds %>%
      filter(.data$C1 %in% atomTest | .data$C2 %in% atomTest) %>%
      filter(!(.data$C1 %in% atomNeg | .data$C2 %in% atomNeg))

    if(is.null(cleaveProd1_bonds)==FALSE){
      keepBonds<-keepBonds %>%
        filter(!(.data$C1 %in% c(rxnSite$Oatom, rxnSite$Patom) & .data$C2 %in% c(rxnSite$Oatom, rxnSite$Patom))) #Do not store the ester bond itself. (Required for alkenes which are inside ring groups. Otherwise the 'while loop' will continue indefinitely.)
    }

    if(nrow(keepBonds) < 1){
      check == FALSE
      break
    }
    atomNeg <- atomTest
    atomCheck<-keepBonds[,1:2]

    atomTest_temp<-unlist(atomCheck)
    atomTest<-atomTest_temp[-which(unlist(keepBonds[,1:2]) %in% atomTest)]

    cleaveProd2_bonds<-bind_rows(cleaveProd2_bonds, keepBonds)
  }

  if(is.null(cleaveProd2_bonds) == FALSE){

    newAtom_empty<-setNames(c(rep(0, ncol(atoms))), colnames(atoms))
    cleaveProd2_atoms<-atoms %>%
      rownames_to_column() %>%
      rowid_to_column() %>%
      filter(.data$rowid %in% cleaveProd2_bonds$C1 | .data$rowid %in% cleaveProd2_bonds$C2) %>%
      bind_rows(newAtom_empty) %>%
      mutate(rowid_old = if_else(is.na(.data$rowid), as.integer(lag(.data$rowid)+1), .data$rowid))  %>%
      mutate(rowid = seq(from=1, to=n(), by=1)) %>%
      mutate(rowname = if_else(is.na(.data$rowname), paste0("O_", .data$rowid), .data$rowname)) %>%
      mutate(rowname = paste0(str_extract(.data$rowname, "[:alpha:]"), "_", .data$rowid))

    newBond_empty<-setNames(c(rep(0, ncol(bonds))), colnames(bonds))
    cleaveProd2_bonds<-cleaveProd2_bonds %>%
      rowid_to_column %>%
      rowwise() %>%
      mutate(C1 = if_else(.data$C1 %in% cleaveProd2_atoms$rowid_old, cleaveProd2_atoms$rowid[which(.data$C1 == cleaveProd2_atoms$rowid_old)], 0)) %>%
      mutate(C2 = if_else(.data$C2 %in% cleaveProd2_atoms$rowid_old, cleaveProd2_atoms$rowid[which(.data$C2 == cleaveProd2_atoms$rowid_old)], 0)) %>%
      column_to_rownames(var = "rowid") %>%
      bind_rows(newBond_empty) %>%
      mutate(
        C1 = replace(.data$C1, list = which(.data$C1 == 0), values = cleaveProd2_atoms$rowid[which(cleaveProd2_atoms$rowid_old == rxnSite$Patom)]),
        C2 = replace(.data$C2, list = which(.data$C2 == 0), values = last(cleaveProd2_atoms$rowid)),
        C3 = replace(.data$C3, list = which(.data$C3 == 0), values = 1)
      ) %>%
      as.matrix()

    cleaveProd2_atoms<-cleaveProd2_atoms %>%
      select(-c(.data$rowid, .data$rowid_old)) %>%
      column_to_rownames(var = "rowname") %>%
      as.matrix()



    cleaveProd2<-list(cleaveProd2_atoms, cleaveProd2_bonds)
    prod2_SDF<-makeSDF(cleaveProd2, currentSDF)
    prod2_SMI<-sdf2smiles(prod2_SDF)
    prod2_SMI<-prod2_SMI@smilist[[1]]

  }else(prod2_SMI<-NULL)



  prodSMI<-c(prod1_SMI, prod2_SMI)

  return(prodSMI)

}


#' @title Predicts the products of ozonolysis at a single alkene group.
#'
#' @param rxnSite One-row tibble containing: molecule (as SMILES string), reaction site type (string), SDF bondblock row number, and SDF atomblock row numbers for phosphorus, oxygen, sulfur, and/or carbon. See findBonds() for more details.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom ChemmineR smiles2sdf atomblock bondblock sdf2smiles
#' @importFrom tidyr unnest
#' @importFrom dplyr distinct filter bind_rows mutate if_else rowwise select lag n last
#' @importFrom tibble rownames_to_column rowid_to_column column_to_rownames
#' @importFrom stringr str_extract
#'
#' @return A character vector of length two containing SMILES strings, one for each ozonolysis product.
predict.Ozonolysis<-function(rxnSite){

  currentSDF<-smiles2sdf(rxnSite$SMILES[1])
  rxnSite<-rxnSite %>% #Extract the rows corresponding to the current molecule and assign reaction types based on functional group
    unnest(c(.data$bondRow, .data$Patom, .data$Oatom, .data$Catom, .data$Satom)) %>%
    distinct()

  atoms<-as.data.frame(atomblock(currentSDF[[1]]))
  bonds<-as.data.frame(bondblock(currentSDF[[1]]))

  check <- TRUE
  atomTest <- rxnSite$Catom[1]
  atomNeg <- rxnSite$Catom[2]

  cleaveProd1_bonds<-NULL
  while(check == TRUE){
    keepBonds<-bonds %>%
      filter(.data$C1 %in% atomTest | .data$C2 %in% atomTest) %>% #Keep any bonds involving the currently targeted atom
      filter(!(.data$C1 %in% atomNeg | .data$C2 %in% atomNeg)) #Remove any bonds involving the currently rejected atom

    if(is.null(cleaveProd1_bonds)==FALSE){
      keepBonds<-keepBonds %>%
        filter(!(.data$C1 %in% c(rxnSite$Catom[1], rxnSite$Catom[2]) & .data$C2 %in% c(rxnSite$Catom[1], rxnSite$Catom[2]))) #Do not store the double bond itself. (Required for alkenes which are inside ring groups. Otherwise the 'while loop' will continue indefinitely.)
    }

    if(nrow(keepBonds) < 1){
      check == FALSE
      break
    }

    atomNeg <- atomTest
    atomCheck<-keepBonds[,1:2]

    atomTest_temp<-unlist(atomCheck)
    atomTest<-atomTest_temp[-which(unlist(keepBonds[,1:2]) %in% atomTest)]

    cleaveProd1_bonds<-bind_rows(cleaveProd1_bonds, keepBonds)
  }

  if(is.null(cleaveProd1_bonds)==FALSE){
    newAtom_empty<-setNames(c(rep(0, ncol(atoms))), colnames(atoms))
    cleaveProd1_atoms<-atoms %>%
      rownames_to_column() %>%
      rowid_to_column() %>%
      filter(.data$rowid %in% cleaveProd1_bonds$C1 | .data$rowid %in% cleaveProd1_bonds$C2) %>%
      bind_rows(newAtom_empty) %>%
      mutate(rowid_old = if_else(is.na(.data$rowid), as.integer(lag(.data$rowid)+1), .data$rowid))  %>%
      mutate(rowid = seq(from=1, to=n(), by=1)) %>%
      mutate(rowname = if_else(is.na(.data$rowname), paste0("O_", .data$rowid), .data$rowname)) %>%
      mutate(rowname = paste0(str_extract(.data$rowname, "[:alpha:]"), "_", .data$rowid))


    newBond_empty<-setNames(c(rep(0, ncol(bonds))), colnames(bonds))
    cleaveProd1_bonds<-cleaveProd1_bonds %>%
      rowid_to_column %>%
      rowwise() %>%
      mutate(C1 = if_else(.data$C1 %in% cleaveProd1_atoms$rowid_old, cleaveProd1_atoms$rowid[which(.data$C1 == cleaveProd1_atoms$rowid_old)], 0)) %>%
      mutate(C2 = if_else(.data$C2 %in% cleaveProd1_atoms$rowid_old, cleaveProd1_atoms$rowid[which(.data$C2 == cleaveProd1_atoms$rowid_old)], 0)) %>%
      column_to_rownames(var = "rowid") %>%
      bind_rows(newBond_empty) %>%
      mutate(
        C1 = replace(.data$C1, list = which(.data$C1 == 0), values = cleaveProd1_atoms$rowid[which(cleaveProd1_atoms$rowid_old == rxnSite$Catom[1])]),
        C2 = replace(.data$C2, list = which(.data$C2 == 0), values = last(cleaveProd1_atoms$rowid)),
        C3 = replace(.data$C3, list = which(.data$C3 == 0), values = 2)
      ) %>%
      as.matrix()

    cleaveProd1_atoms<-cleaveProd1_atoms %>%
      select(-c(.data$rowid, .data$rowid_old)) %>%
      column_to_rownames(var = "rowname") %>%
      as.matrix()


    cleaveProd1<-list(cleaveProd1_atoms, cleaveProd1_bonds)
    prod1_SDF<-makeSDF(cleaveProd1, currentSDF)
    prod1_SMI<-sdf2smiles(prod1_SDF)
    prod1_SMI<-prod1_SMI@smilist[[1]]
  }else{prod1_SMI<-NULL}




  ##Now do the same thing for the other side of the alkene bond
  check <- TRUE
  atomTest <- rxnSite$Catom[2]
  atomNeg <- rxnSite$Catom[1]

  cleaveProd2_bonds<-NULL
  while(check == TRUE){
    keepBonds<-bonds %>%
      filter(.data$C1 %in% atomTest | .data$C2 %in% atomTest) %>% #Keep any bonds involving the currently targeted atom
      filter(!(.data$C1 %in% atomNeg | .data$C2 %in% atomNeg)) #Remove any bonds involving the currently rejected atom

    if(is.null(cleaveProd2_bonds)==FALSE){
      keepBonds<-keepBonds %>%
        filter(!(.data$C1 %in% c(rxnSite$Catom[1], rxnSite$Catom[2]) & .data$C2 %in% c(rxnSite$Catom[1], rxnSite$Catom[2]))) #Do not store the double bond itself. (Required for alkenes which are inside ring groups. Otherwise the 'while loop' will continue indefinitely.)
    }

    if(nrow(keepBonds) < 1){
      check == FALSE
      break
    }
    atomNeg <- atomTest
    atomCheck<-keepBonds[,1:2]

    atomTest_temp<-unlist(atomCheck)
    atomTest<-atomTest_temp[-which(unlist(keepBonds[,1:2]) %in% atomTest)]

    cleaveProd2_bonds<-bind_rows(cleaveProd2_bonds, keepBonds)
  }

  if(is.null(cleaveProd2_bonds)==FALSE){
    newAtom_empty<-setNames(c(rep(0, ncol(atoms))), colnames(atoms))
    cleaveProd2_atoms<-atoms %>%
      rownames_to_column() %>%
      rowid_to_column() %>%
      filter(.data$rowid %in% cleaveProd2_bonds$C1 | .data$rowid %in% cleaveProd2_bonds$C2) %>%
      bind_rows(newAtom_empty) %>%
      mutate(rowid_old = if_else(is.na(.data$rowid), as.integer(lag(.data$rowid)+1), .data$rowid))  %>%
      mutate(rowid = seq(from=1, to=n(), by=1)) %>%
      mutate(rowname = if_else(is.na(.data$rowname), paste0("O_", .data$rowid), .data$rowname)) %>%
      mutate(rowname = paste0(str_extract(.data$rowname, "[:alpha:]"), "_", .data$rowid))

    newBond_empty<-setNames(c(rep(0, ncol(bonds))), colnames(bonds))
    cleaveProd2_bonds<-cleaveProd2_bonds %>%
      rowid_to_column %>%
      rowwise() %>%
      mutate(C1 = if_else(.data$C1 %in% cleaveProd2_atoms$rowid_old, cleaveProd2_atoms$rowid[which(.data$C1 == cleaveProd2_atoms$rowid_old)], 0)) %>%
      mutate(C2 = if_else(.data$C2 %in% cleaveProd2_atoms$rowid_old, cleaveProd2_atoms$rowid[which(.data$C2 == cleaveProd2_atoms$rowid_old)], 0)) %>%
      column_to_rownames(var = "rowid") %>%
      bind_rows(newBond_empty) %>%
      mutate(
        C1 = replace(.data$C1, list = which(.data$C1 == 0), values = cleaveProd2_atoms$rowid[which(cleaveProd2_atoms$rowid_old == rxnSite$Catom[2])]),
        C2 = replace(.data$C2, list = which(.data$C2 == 0), values = last(cleaveProd2_atoms$rowid)),
        C3 = replace(.data$C3, list = which(.data$C3 == 0), values = 2)
      ) %>%
      as.data.frame() %>%
      as.matrix()

    cleaveProd2_atoms<-cleaveProd2_atoms %>%
      select(-c(.data$rowid, .data$rowid_old)) %>%
      column_to_rownames(var = "rowname") %>%
      as.matrix()


    cleaveProd2<-list(cleaveProd2_atoms, cleaveProd2_bonds)
    prod2_SDF<-makeSDF(cleaveProd2, currentSDF)
    prod2_SMI<-sdf2smiles(prod2_SDF)
    prod2_SMI<-prod2_SMI@smilist[[1]]

  }else{prod2_SMI<-NULL}


  prodSMI<-c(prod1_SMI, prod2_SMI)

  return(prodSMI)
}


#' @title Predicts the secondary ozonide product of ozonolysis at a single alkene group.
#'
#' @param rxnSite One-row tibble containing: molecule (as SMILES string), reaction site type (string), SDF bondblock row number, and SDF atomblock row numbers for phosphorus, oxygen, sulfur, and/or carbon. See findBonds() for more details.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom ChemmineR smiles2sdf atomblock bondblock sdf2smiles
#' @importFrom tidyr unnest
#' @importFrom dplyr distinct bind_rows rowwise mutate if_else select slice
#' @importFrom stats setNames
#' @importFrom tibble rownames_to_column rowid_to_column column_to_rownames
#'
#' @return A character vector of length one containing a SMILES string.
predict.SecOzon<-function(rxnSite){ #Formation of a secondary ozonide from an alkene.
  currentSDF<-smiles2sdf(rxnSite$SMILES[1])
  rxnSite<-rxnSite %>% #Extract the rows corresponding to the current molecule and assign reaction types based on functional group
    unnest(c(.data$bondRow, .data$Patom, .data$Oatom, .data$Catom, .data$Satom)) %>%
    distinct()

  atoms<-as.data.frame(atomblock(currentSDF[[1]]))
  bonds<-as.data.frame(bondblock(currentSDF[[1]]))

  newAtom_empty<-setNames(c(rep(0, ncol(atoms))), colnames(atoms)) #Presently ignoring 3d positioning and treating the new oxygens as all being on the same plane.

  atomTable_temp <- atoms %>% #Need to add the 3 new oxygen atoms
    rownames_to_column() %>%
    bind_rows(newAtom_empty, newAtom_empty, newAtom_empty) %>%
    rowid_to_column() %>%
    rowwise() %>%
    mutate(rowname = if_else(is.na(.data$rowname), paste0("O_", .data$rowid), .data$rowname)) %>%
    select(-.data$rowid) %>%
    column_to_rownames(var = "rowname") %>%
    as.matrix()


  newAtoms_ID<-c(nrow(atomTable_temp)-2, nrow(atomTable_temp)-1, nrow(atomTable_temp))
  newBond_empty<-setNames(c(rep(0, ncol(bonds))), colnames(bonds))

  newBonds_C1<-c(newAtoms_ID[1], newAtoms_ID[1], newAtoms_ID[2], newAtoms_ID[2], newAtoms_ID[3]) #O1, O1, O2, O2, O3 (new oxygen atoms).
  newBonds_C2<-as.numeric(c(bonds[unique(rxnSite$bondRow), 1], bonds[unique(rxnSite$bondRow), 2], bonds[unique(rxnSite$bondRow), 1], newAtoms_ID[3], bonds[unique(rxnSite$bondRow), 2])) #C1, C2, C1, O3, C2.
  #When added to the Bonds table below, should result in the following (ozonide) bonds: O1-C1, O1-C2, O2-C1, O2-O3, O3-C2.

  bondTable_temp<-bonds %>% #Need to create 5 new bonds (C-O x4, O-O x1) and delete the original C=C bond
    bind_rows(newBond_empty, newBond_empty, newBond_empty, newBond_empty, newBond_empty) %>% #Add 5 new (empty) bonds
    mutate(
      C1 = replace(.data$C1, list = which(.data$C1 == 0), values = newBonds_C1),
      C2 = replace(.data$C2, list = which(.data$C2 == 0), values = newBonds_C2),
      C3 = replace(.data$C3, list = which(.data$C3 == 0), values = 1)
    ) %>%
    slice(-rxnSite$bondRow) %>%
    rowid_to_column() %>%
    column_to_rownames(var = "rowid") %>%
    as.matrix()

  atoms_bonds<-list(atomTable_temp, bondTable_temp)
  currentSDF<-makeSDF(atoms_bonds, currentSDF)
  if(!is.null(currentSDF)){
    prodSMI<-sdf2smiles(currentSDF)
    prodSMI<-prodSMI@smilist[[1]]
  }else{prodSMI<-NULL}


  return(prodSMI)
}



#' @title Predicts the product of addition of HOCl to a single alkene group.
#'
#' @param rxnSite One-row tibble containing: molecule (as SMILES string), reaction site type (string), SDF bondblock row number, and SDF atomblock row numbers for phosphorus, oxygen, sulfur, and/or carbon. See findBonds() for more details.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom ChemmineR smiles2sdf atomblock bondblock bonds sdf2smiles
#' @importFrom tidyr unnest
#' @importFrom dplyr distinct
#'
#' @return A character vector of length one containing a SMILES string.
predict.HOClAdd<-function(rxnSite){

  currentSDF<-smiles2sdf(rxnSite$SMILES[1])
  rxnSite<-rxnSite %>% #Extract the rows corresponding to the current molecule and assign reaction types based on functional group
    unnest(c(.data$bondRow, .data$Catom)) %>%
    distinct()

  atomTable_temp<-atomblock(currentSDF[[1]])
  bondTable_temp<-bondblock(currentSDF[[1]])

  allBonds_count<-bonds(currentSDF[[1]]) #The bonds() function returns a simplified bond table containing the total number of bonds for each atom in the molecule.

  Catoms<-unlist(rxnSite$Catom)
  CBonds_count<-allBonds_count[Catoms,]

  if(length(unique(CBonds_count$Nbondcount))>1){ #If one carbon is more highly substituted, the reaction will be regioselective.
    #Cl (electrophile) will bond to the less substituted carbon, while OH (nucleophile) will bond to the more substituted carbon (Markovnikov's Rule).

    Catom_OHtarget<-Catoms[which.max(CBonds_count$Nbondcount)]
    Catom_Cltarget<-Catoms[which.min(CBonds_count$Nbondcount)]

  }else{ #If they are equally substituted, then we don't care too much about the isomer formed (for now)

    Catom_OHtarget<-Catoms[1]
    Catom_Cltarget<-Catoms[2]
  }

  atomTable_temp<-rbind(atomTable_temp, c(rep(0, ncol(atomTable_temp)))) #Adds a new row to the atoms table for the C-OH bond
  newAtom_ID<-nrow(atomTable_temp)
  rownames(atomTable_temp)[newAtom_ID]<-paste0("O_", newAtom_ID) #Labels the row as O_#, where # is the current row index

  bondTable_temp<-rbind(bondTable_temp, c(Catom_OHtarget, newAtom_ID, 1, rep(0, ncol(bondTable_temp)-3))) #Adds a new row to the bonds table between the target C and the O added in the previous step (the "1" means single bond)
  rownames(bondTable_temp)[nrow(bondTable_temp)]<-as.character(nrow(bondTable_temp)) #bonds table rows are just labeled with the row index


  atomTable_temp<-rbind(atomTable_temp, c(rep(0, ncol(atomTable_temp)))) #Adds a new row to the atoms table for the C-Cl bond
  newAtom_ID<-nrow(atomTable_temp)
  rownames(atomTable_temp)[newAtom_ID]<-paste0("Cl_", newAtom_ID) #Labels the row as Cl_#, where # is the current row index

  bondTable_temp<-rbind(bondTable_temp, c(Catom_Cltarget, newAtom_ID, 1, rep(0, ncol(bondTable_temp)-3))) #Adds a new row to the bonds table between the target C and the O added in the previous step (the "1" means single bond)
  rownames(bondTable_temp)[nrow(bondTable_temp)]<-as.character(nrow(bondTable_temp)) #bonds table rows are just labeled with the row index

  bondTable_temp[unique(rxnSite$bondRow),3]<-1 #Changes the old C=C bond to a C-C bond
  atoms_bonds<-list(atomTable_temp, bondTable_temp)
  currentSDF<-makeSDF(atoms_bonds, currentSDF)

  if(!is.null(currentSDF)){
    prodSMI<-sdf2smiles(currentSDF)
    prodSMI<-prodSMI@smilist[[1]]
  }else{prodSMI<-NULL}

  return(prodSMI)
}

