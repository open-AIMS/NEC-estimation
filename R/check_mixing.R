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

#' check.mixing
#'
#' Calculates the within and between chain mixing and tests if the between is greater than within based on the supplied quantile probability value.
#'
#' @param  J a jag model fit as returned by a call to jags
#'
#' @param  prob.val the probability value to pass to quantile to assess the between 
#' chain coefficient of variation against the within chain coefficient of variation. 
#' Please inspect the model code for more details.
#'
#' @export
#' @return A list of within and between coefficient of variations for parameter, 
#' the outcome of the test if between is greater than within, and the ratio.

check.mixing <- function(J, prob.val=0.8){

  x <- J$BUGSoutput$sims.array
  params <- J$parameters.to.save
  num.chains <- ncol(x[,,params[1]])
  within.cv <- lapply(params, FUN=function(i){
     x1          <- as.vector(x[,,i])
     chain.id     <- rep(1:num.chains, each = nrow(x[,,i]))
     unlist(lapply(1:num.chains,FUN=function(k){
        sd(x1[which(chain.id==k)])/mean(x1[which(chain.id==k)])     
      }))       
    }) 
  names(within.cv) <- params
  
  between.cv <- lapply(params, FUN=function(i){
     x1          <- as.vector(x[,,i])
     chain.id     <- rep(1:num.chains, each = nrow(x[,,i]))
     sd(x1)/mean(x1) 
    }) 
  names(between.cv) <- params
  
  cv.test <- lapply(params, FUN=function(i){
    between.cv[[i]]>quantile(within.cv[[i]], probs = prob.val)
  })
  names(cv.test) <- params
    
  cv.ratio <- mean(unlist(between.cv)/unlist(lapply(within.cv, FUN=mean))) 
    
  return(list('within.cv'=within.cv,
              'between.cv'=between.cv,
              'cv.test'=cv.test,
              'cv.ratio'=cv.ratio))  
}

