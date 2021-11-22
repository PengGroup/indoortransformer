#' @title Recursive prediction of transformation products.
#'
#' @param SMI Character vector of length one containing a SMILES string.
#' @param prods Empty three-column tibble in which to store results.
#' @param funGroups Character vector containing functional groups at which transformation reactions can occur.
#' @param rxnList Character vector containing types of reactions to use for prediction.
#' @param rxnSave Empty character vector for storing sequential reaction information.
#' @param savePoint Single logical value. Value determines whether to perform hydrolysis or ozonolysis on a given alkene group. Recursion allows for both reactions to be performed on a single alkene group.
#' @param skipCheck Empty character vector for storing SMILES strings of predicted products. Used to improve efficiency by not performing the same prediction multiple times.
#'
#' @importFrom tibble tibble
#' @importFrom stats na.exclude
#' @importFrom dplyr filter bind_rows
#'
#' @return Tibble containing product SMILES, list of sequential reactions for each product, and saved skipCheck strings.
recursive.Predict<-function(SMI, prods, funGroups, rxnList, rxnSave, savePoint, skipCheck){
  prodSMILES<-NULL

  rxnGroups<-smartBonds(SMI, funGroups)

  if(length(rxnGroups)<1){
    return(NULL)
  }else{
    prods_temp<-tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL) #Creates final output table to store predicted transformation product information. This table will be returned at the end of this function.
    prods_temp2<-tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
    prods_Cleave2 <- tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
    prods_Cleave1 <- tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
    prods_Cleave3 <- tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)

    for(j in 1:length(rxnGroups)){

      prod_new<-predict.Product(SMI, rxnGroups[j], rxnList, rxnSave, savePoint) #'predict.Product' should return a multi-row tibble for cleavage reactions, and a one-row tibble for all other reactions
      if(is.null(prod_new)){next}

      if(any(prod_new$prodSMILES %in% skipCheck)){prod_new$prodSMILES[which(prod_new$prodSMILES %in% skipCheck)]<-NA} #If any of the predicted products are duplicates of previously predicted products, remove them from 'prod_new'

      skipCheck<-skipCheck %>%
        c(prod_new$prodSMILES) %>%
        na.exclude() %>%
        unique()

      if(last(unlist(prod_new$Reactions[1])) %in% c("Hyd", "Ozo")){
        prod_new2 <- prod_new[2,] #Second cleavage product
        prod_new1 <- prod_new[1,] #First cleavage product

        rxnSave_current<-prod_new2$Reactions
        savePoint_current<-prod_new2$savePoint

        prods_Cleave2 <- tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
        prods_Cleave2 <- bind_rows(prods_Cleave2, prod_new2, recursive.Predict(prod_new2$prodSMILES[1], prods, funGroups, rxnList, rxnSave_current, savePoint_current, skipCheck))
        prods_Cleave2 <- filter(prods_Cleave2, !is.na(prodSMILES))
        skipCheck<-c(skipCheck, unlist(last(prods_Cleave2$skipCheck))) #Using the last row to extract SMILES to skip because the first row(s) will be NULL
        skipCheck<-unique(skipCheck)

        rxnSave_current<-prod_new1$Reactions
        savePoint_current<-prod_new1$savePoint

        prods_Cleave1 <- tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
        prods_Cleave1 <- bind_rows(prods_Cleave1, prod_new1, recursive.Predict(prod_new1$prodSMILES[1], prods, funGroups, rxnList, rxnSave_current, savePoint_current, skipCheck))
        prods_Cleave1 <- filter(prods_Cleave1, !is.na(prodSMILES))
        skipCheck<-c(skipCheck, unlist(last(prods_Cleave1$skipCheck)))
        skipCheck<-unique(skipCheck)


        prod_new3<-NULL
        if(last(unlist(prod_new$Reactions[1])) %in% c("Ozo")){
          prod_new3<-prod_new[3,]
          rxnSave_current<-prod_new3$Reactions
          savePoint_current<-prod_new3$savePoint

          prods_Cleave3 <- tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
          prods_Cleave3 <- bind_rows(prods_Cleave3, prod_new3, recursive.Predict(prod_new3$prodSMILES[1], prods, funGroups, rxnList, rxnSave_current, savePoint_current, skipCheck))
          prods_Cleave3 <- filter(prods_Cleave3, !is.na(prodSMILES))
          skipCheck<-c(skipCheck, unlist(last(prods_Cleave3$skipCheck)))
          skipCheck<-unique(skipCheck)

        }
      }else{
        rxnSave_current<-prod_new$Reactions[1]
        savePoint_current<-prod_new$savePoint[1]

        prods_temp<-bind_rows(prods_temp, prod_new, recursive.Predict(prod_new$prodSMILES[1], prods, funGroups, rxnList, rxnSave_current, savePoint, skipCheck))
        prods_temp <- filter(prods_temp, !is.na(prodSMILES))
        skipCheck<-c(skipCheck, unlist(last(prods_temp$skipCheck)))
        skipCheck<-unique(skipCheck)

        prods_temp2<-tibble(Reactions=NULL, prodSMILES=NULL, skipCheck=NULL)
        if(savePoint_current == TRUE){
          prods_temp2<-bind_rows(prods_temp2, prod_new, recursive.Predict(SMI, prods, funGroups, rxnList, rxnSave, savePoint_current, skipCheck))
          prods_temp2 <- filter(prods_temp2, !is.na(prodSMILES))
          skipCheck<-c(skipCheck, unlist(last(prods_temp2$skipCheck)))
          skipCheck<-unique(skipCheck)
          savePoint<-FALSE
        }

      }
      prods_temp<-bind_rows(prods_temp, prods_temp2, prods_Cleave1, prods_Cleave2, prods_Cleave3)
      prods_temp$skipCheck<-list(skipCheck)
      prods_temp<-filter(prods_temp, !is.na(prodSMILES))
    }

    return(prods_temp)
  }
}

