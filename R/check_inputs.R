#    Copyright 2020 Australian Institute of Marine Science
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

#' check_inputs
#'
#' Checks the input data for a jagsNEC model call
#'
#' @param  data a data.frame containing the data to use for the model
#' 
#' @param x.var the column heading indicating the concentration (x) variable
#' 
#' @param y.var the column heading indicating the response (y) variable
#' 
#' @param trials.var the column heading indicating the column for the number of "trials" for binomial response data. 
#' If not supplied, the model may run but will not be the model you intended!  

check_inputs <- function(data,
                         x.var,
                         y.var,
                         trials.var){
  # check the specified columns exist in data
  use.vars <- na.omit(c(y.var=y.var, x.var=x.var, trials.var))
  var.colms <- match(use.vars, colnames(data))
  missing.colms <- data.frame(val=use.vars[which(is.na(var.colms))], stringsAsFactors = FALSE)
  missing.colms$element <- rownames(missing.colms)
  if(length(na.omit(var.colms))<length(use.vars)){
    stop(paste("Your indicated ", paste(paste(missing.colms$element, " '", missing.colms$val,"'", sep=""),
                                        collapse = ", "),
               " is not present in your input data. Has this been mispecified?", sep=""))
  }
  
  # extract the data
  y.dat <- data[, y.var]
  x.dat <- data[, x.var]
  
  # check the x data are numeric
  if(class(x.dat)!="numeric"){
    stop(paste("Your indicated x.var column ", x.var," contains data that is class ", class(x.dat),".
                   The function jagsNEC requires the concentration data (argument x.var) to be numeric.",sep=""))
  }
  
  # check data contains only finite values
  test.x <- mean(x.dat)
  test.y <- mean(y.dat)
  if(is.finite(test.x)!=TRUE){
    stop("Your x.var column contains values that are not finite.")
  }
  if(is.finite(test.y)!=TRUE){
    stop("Your y.var column contains values that are not finite.")
  }
  
  # check the data are lower at high x compared to low x (ie the response variable declines within increase in the x)
  if(mean(y.dat[which(x.dat<mean(x.dat))])< mean(y.dat[which(x.dat>mean(x.dat))]) & model != "NECHormesis"){
    stop("The mean value of the response for the lower half of the
             concentration data are lower than that of the upper half of the concentration data.
             jagsNEC only fits concentration response data where the
             response declines with increasing values of concentration.")
 }
}

