#!/usr/bin/env Rscript

### Iterations new method for two-class classification model
### Author: Raquel Bosch Romeu

pacman::p_load(ExPosition, ks)

# load database
wd <- ""
load(paste0(wd,"db_nps.RData"))
db <- db_nps  # Data base used

# Parameters
n_iterations = 1000
n_model <- 270   # number of individuals to build the map
n_test <- nrow(db_nps) - n_model    # number of individuals to validate the map

# Order of the dummy variables
predictors_order           <- c("Sex.Man", "Sex.Woman", "More_than_one_SPN.No",   
                                "More_than_one_SPN.Yes", "Diameter.Small", 
                                "Diameter.Medium", "Diameter.Large", 
                                "NPS_localisation.Upper", 
                                "NPS_localisation.Middle", 
                                "NPS_localisation.Lower",
                                "NPS_border.Smooth", "NPS_border.Spiculation",
                                "NPS_border.Lobulation", "NPS_border.Other", 
                                "Previous_malignancy.No", 
                                "Previous_malignancy.Yes", "Smoker.Never",           
                                "Smoker.Yes", "COPD.No", "COPD.Yes")
response_order <- c("dxlungcancer.0", "dxlungcancer.1")


#' Data preparation
#' @description 
#' Transforms the data.frame with predictor and response 
#' variables into a data frame with dummy variables. 
#' @param db Database predictor and response variables. The last 
#' column must contain the response variable
#' @return The database with the dummy variables. the antepenultimate 
#' and penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test) 
#' to which the individual has been assigned.
#' @import ExPosition
DataPreparation = function(db) {
  # divide the database into training and validation groups
  index_model <- sample(nrow(db), size=n_model)
  index_test <- setdiff(1:nrow(db),index_model)
  group <- as.vector(rep(NA, nrow(db)))
  group[index_model] <- rep("model", length(index_model))
  group[index_test] <- rep("test", length(index_test))
  
  # convert variables into dummy variables
  db_dummy <- makeNominalData(db)
  db_dummy <- db_dummy[, c(predictors_order, response_order)]   
  db_dummy <- data.frame(db_dummy, group)
  colnames(db_dummy)[c(length(db_dummy)-2, length(db_dummy)-1)] <- c("class0",
                                                                     "class1")
  
  return(db_dummy)
}

#' Construction of the base
#' @description 
#' It gets the new base obtained after the singular value decomposition
#' @param db_dummy The database with the dummy variables. The antepenultimate 
#' and penultimate column must correspond to the response variable. 
#' The last column must correspond to the group (training or test) 
#' to which the individual has been assigned. 
#' @param variables_order The order of the dummy variables 
#' @return The base obtained after the singular value decomposition 
#' of the density matrix. Only vectors whose singular values are 
#' greater than 0 are retained.
BaseConstruction= function(db_dummy, variables_order) {
  predictors <- as.matrix(db_dummy[db_dummy$group == "model", 
                                   c(1:(length(db_dummy)-3))])
  response <- as.matrix(db_dummy[db_dummy$group == "model", 
                                 c("class0", "class1")])
  D <- t(response)%*%predictors
  X <- sqrt(D)
  
  F <- t(X)%*%X
  F <- (1/sum(diag(F)))*F  # density matrix
  F <- F[, predictors_order]
  F <- F[predictors_order, ]
  
  svd_result <- svd(F)
  
  base <-  svd(F)$u[, which(svd_result$d > .Machine$double.eps)]
}

#' Calculation of the new coordinates
#' @description 
#' It calculates the coordinates in the new basis of the vectors 
#' of the categorical co-variables
#' @param db_dummy The database with the dummy variables. the antepenultimate 
#' and penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test) 
#' to which the individual has been assigned. 
#' @param base The base obtained after the singular value decomposition 
#' of the density matrix
#' @return The data frame with the coordinates in the new basis. The 
#' penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test).
CalculateNewCoordinates <- function(db_dummy, base) {
  predictors <- as.matrix(db_dummy[, c(1:(length(db_dummy)-3))])
  new_coordinates <- t(apply(predictors, 1, 
                             function(y) t(base)%*%(y/sqrt(sum(y)))))
  new_coordinates <- data.frame(new_coordinates, db_dummy$class1, 
                                db_dummy$group) 
  colnames(new_coordinates) <- c("coor1", "coor2", "response", "group")
  
  return(new_coordinates)
}

#' Calculation of the polar coordinates
#' @description 
#' It calculates the polar coordinates from the data frame containing 
#' the coordinate information in the new base.
#' @param new_coordinates The data frame with the coordinates in the new basis. 
#' The penultimate column must correspond to the response variable. 
#' The last column must correspond to the group (training or test). 
#' @return The data frame with the polar coordinates. The 
#' penultimate column corresponds to the response variable. 
#' The last column corresponds to the group (training or test).
CalculatePolarCoordinates <- function(new_coordinates) {
  theta <- atan(new_coordinates[,"coor2"]/new_coordinates[,"coor1"])
  rho <- sqrt((new_coordinates[,"coor1"]^2) + (new_coordinates[,"coor2"]^2))
  
  polar_coor <- data.frame(theta, rho, db_dummy$class1,  db_dummy$group)
  colnames(polar_coor) <- c("theta", "rho", "response", "group")
  
  return(polar_coor)
} 

#' Calculation of the maximum of the density functions of each class
#' @description 
#' It calculates the maximum of the density functions of the coordinate(s)
#' of the "model" group (training set) of each class. If two-coordinate maxima 
#' are calculated, the Gaussian kernel is used, if one-coordinate maxima 
#' are calculated, the Epanechnikov kernel is used.
#' @param coor The data frame  with the coordinates for which you want to 
#' calculate the density maxima. The penultimate column must correspond to 
#' the response variable. The last column must correspond to the 
#' group (training or test).
#' @param Coordinates Coordinates for which the maxima of density functions
#' are to be calculated. By default, it is calculated for two coordinates, 
#' corresponding to the first two columns of the "coor" data frame. 
#' If you only want to do it for one, you must indicate the name of 
#' the column containing that coordinate.
#' @return Matrix with the maxima of the density functions of each class (rows) 
#' @import ks
MaximumDensityFunctions <- function(coor, coordinates = "two") {
  if (coordinates == "two") {
    den_class0 <- kde(coor[coor$group == "model" & coor$response == "0", 
                           c(1,2)])
    # Position of the maximum density value
    max_pos_class0 <- which(den_class0$estimate == max(den_class0$estimate), 
                            arr.ind = TRUE)[1,]
    # Value of the coordinates at the point of maximum density
    class0 <- c(den_class0$eval.points[[1]][max_pos_class0[1]], 
                den_class0$eval.points[[2]][max_pos_class0[2]])
    
    den_class1 <- kde(coor[coor$group == "model" & coor$response == "1", 
                           c(1,2)])
    max_pos_class1 <- which(den_class1$estimate == max(den_class1$estimate), 
                            arr.ind = TRUE)[1,]
    class1 <- c(den_class1$eval.points[[1]][max_pos_class1[1]], 
                den_class1$eval.points[[2]][max_pos_class1[2]])
    
    maxima <- rbind(class0, class1)
    colnames(maxima) <- colnames(coor[,c(1,2)]) 
    
    return(maxima)
  }
  else {
    den_class0 <- density(coor[coor$group == "model" & coor$response == "0", 
                               c(coordinates)], kernel = "epanechnikov")
    # Position of the maximum density value
    max_pos_class0 <- which.max(den_class0$y)
    # Value of the coordinates at the point of maximum density
    class0 <- den_class0$x[max_pos_class0] 
    
    den_class1 <- density(coor[coor$group == "model" & coor$response == "1", 
                               c(coordinates)], kernel = "epanechnikov")
    max_pos_class1 <- which.max(den_class1$y)
    class1 <- den_class1$x[max_pos_class1] 
    
    maxima <- rbind(class0, class1)
    colnames(maxima) <- coordinates
    
    return(maxima)
  }
}

#' Assignment of the map-predicted class
#' @description 
#' It calculates for each individual in the validation group ("test") the 
#' class predicted by the map. It builds a data frame with the prediction 
#' made, the actual classification, whether these are matched and the 
#' score used for discrimination.
#' @param coor The data frame with coordinates of individuals. The penultimate 
#' column must correspond to the response variable. The last column must 
#' correspond to the group (training or test).
#' @param maxima_points A list from the function "MaximumDensityFunctions" 
#' with the maxima of the density functions of each class (class0 and
#' class1). 
#' @param coordinates Coordinates used in the discrimination. 
#' By default, they are two coordinates, corresponding to the first 
#' two columns of the "coor" data frame. If you only want to do it for one, 
#' you must indicate the name of the column containing that coordinate.
#' @return Data frame with the prediction made, the actual classification, 
#' and the score used for discrimination (the difference between the distance 
#' of the coordinates of each individual with each of the maxima of the 
#' two classes).
ClassAssignment <- function(coor, maxima_points, coordinates = "two") {
  coor_test_set <- coor[coor$group == "test",]
  
  if (coordinates == "two") {
    # Differences with the maximums of each class
    distance_to_class0 <- sqrt((coor_test_set[,1] - maxima_points["class0",1])^2 
                               + (coor_test_set[,2] 
                                  - maxima_points["class0", 2])^2)  
    distance_to_class1 <- sqrt((coor_test_set[,1] - maxima_points["class1",1])^2 
                               + (coor_test_set[,2] 
                                  - maxima_points["class1", 2])^2)  
  }
  else {
    distance_to_class0 <- abs(coor_test_set[,coordinates] 
                              - maxima_points["class0", coordinates]) 
    distance_to_class1 <- abs(coor_test_set[,coordinates] 
                              - maxima_points["class1", coordinates]) 
  }
  difference <- distance_to_class0 - distance_to_class1
  # each individual is assigned to the class with the 
  # maximum density closest to its coordinates.
  prediction <- ifelse(difference < 0, 0, 1)
  prediction_dt <- data.frame(coor_test_set$response, prediction, difference)
  colnames(prediction_dt) <- c("real", "predicted", "difference")
  
  return(prediction_dt)
}

#' Calculation of the predictive parameters
#' @description 
#' It calculates from the data frame constructed by the function 
#' "ClassAssignment "the values of the confusion matrix in relative 
#' frequencies and other predictive parameters (predictive values, accuracy, 
#' F1 score, Cohen's Kappa and Matthews Correlation coefficients).
#' @param prediction_dt Data frame with the prediction made, the actual 
#' classification, and the score used for discrimination (the difference 
#' between the distance of the coordinates of each individual with each 
#' of the maxima of the two classes).
#' @return List of calculated predictive parameters
PredictiveParameters <- function(prediction_dt) {
  TN <- nrow(prediction_dt[prediction_dt$real == 0 
                           & prediction_dt$predicted == 0, ])
  TP <- nrow(prediction_dt[prediction_dt$real == 1 
                           & prediction_dt$predicted == 1, ])
  FN <- nrow(prediction_dt[prediction_dt$real == 1 
                           & prediction_dt$predicted == 0, ])
  FP <- nrow(prediction_dt[prediction_dt$real == 0 
                           & prediction_dt$predicted == 1, ])
  TN_rel <- TN/(TN + FP)
  TP_rel <- TP/(TP + FN)
  FN_rel <- FN/(TP + FN)
  FP_rel <- FP/(TN + FP)
  
  NPV <- TN/(TN + FN)
  PPV <- TP/(TP + FP)
  
  accuracy <- (TN + TP)/(TN + TP + FN + FP)
  f1_score <- (2*TP)/((2*TP) + FN + FP)
  cohens_kappa <- 2*((TP*TN)-(FN*FP))/((TP+FP)*(FP+TN)+(TP+FN)*(FN+TN))
  mcc <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  return(list(confusion_matrix =(c(TN = TN_rel, FP = FP_rel, 
                                   FN = FN_rel, TP = TP_rel)),
              predictive_parameters = c(NPV = NPV, PPV = PPV, 
                                        accuracy = accuracy, 
                                        f1_score = f1_score, 
                                        cohens_kappa = cohens_kappa, 
                                        mcc = mcc)))
}


# MAIN CODE

# Initiation of variables to save results
confusion_matrix_two_coor <- matrix(data = NA, nrow = n_iterations, ncol = 4, 
                                   byrow = FALSE)
colnames(confusion_matrix_two_coor) <- c("TN", "FP", "FN", "TP")
confusion_matrix_theta <- matrix(data = NA, nrow = n_iterations, ncol = 4, 
                                 byrow = FALSE)
colnames(confusion_matrix_theta) <- c("TN", "FP", "FN", "TP")
prediction_parameters_two_coor <- matrix(data = NA, nrow = n_iterations, 
                                         ncol = 6, byrow = FALSE)
colnames(prediction_parameters_two_coor) <- c("NPV", "PPV", "accuracy", 
                                              "f1_score", "cohens_kappa", "mcc")
prediction_parameters_theta <- matrix(data = NA, nrow = n_iterations, ncol = 6, 
                                      byrow = FALSE)
colnames(prediction_parameters_theta) <- c("NPV", "PPV", "accuracy", 
                                           "f1_score", "cohens_kappa", "mcc")
# the difference between the distance to the two maximums is saved.  
# It will be used to construct the ROC and precision-recall curves
dist_diff_two_coor <- list(score = list(), class = list())
dist_diff_theta <- list(score = list(), class = list())


for (i in 1:n_iterations) {
  db_dummy <- DataPreparation(db)
  base <- BaseConstruction(db_dummy, variables_order)
  new_coordinates <- CalculateNewCoordinates(db_dummy, base)
  polar_coor <- CalculatePolarCoordinates(new_coordinates) 

  # Analysis with the two polar coordinates
  maxima_two_coordinates  <- MaximumDensityFunctions(polar_coor)
  prediction_two_coor <- ClassAssignment(polar_coor, maxima_two_coordinates)
  results_two_coor <- PredictiveParameters(prediction_two_coor)
  # Save results 
  confusion_matrix_two_coor[i,] <- results_two_coor$confusion_matrix
  prediction_parameters_two_coor[i, ] <- results_two_coor$predictive_parameters
  # The differences between the distance to the maxima of each class is saved
  dist_diff_two_coor$score[[i]] <- prediction_two_coor$difference
  dist_diff_two_coor$class[[i]] <- prediction_two_coor$real
  
  # Analysis with the angular coordinate  
  maxima_theta  <- MaximumDensityFunctions(polar_coor, "theta")
  prediction_theta <- ClassAssignment(polar_coor, maxima_theta, "theta")
  results_theta <- PredictiveParameters(prediction_theta)
  
  confusion_matrix_theta[i,] <- results_theta$confusion_matrix
  prediction_parameters_theta[i, ] <- results_theta$predictive_parameters
  dist_diff_theta$score[[i]] <- prediction_theta$difference
  dist_diff_theta$class[[i]] <- prediction_theta$real
}

save(confusion_matrix_two_coor, prediction_parameters_two_coor, 
     confusion_matrix_theta, prediction_parameters_theta,
     file = paste0(wd, "results_analysis.RData"))

save(dist_diff_two_coor, dist_diff_theta, 
     file = paste0(wd, "results_dist_diff.RData"))

# The coordinate values of the last run are saved for example graphs
save(polar_coor, file = paste0(wd, "example_coordinates.RData"))