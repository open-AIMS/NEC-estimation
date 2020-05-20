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

#'modify.jagsMANEC
#'
#' Modifys an existing jagsMANECfit, for example, but adding or removing fitted models.
#'
#' @param  jagsMANECfit a jagsMANECfit output list, as returned by fit.jagsMANEC
#'
#' @param model.set A character vector containing the of names of model types to be included in the modified fit.
#' 
#' @param drop.models A character vector containing the names of model types you which to drop for the modified fit.
#' 
#' @param add.models A character vector containing the names of model types to add to the modified fit.
#'
#' @export
#' @return All successully fitted jagsMANECfit model fit. 
modify.jagsMANEC <- function(jagsMANECfit, model.set=NA, drop.models=NA, add.models=NA){
  
  # if the model set is NA 
  if(is.na(model.set)){
   model.set <- names(jagsMANECfit$mod.fits)
  }
  
  if(model.set=="NEC"){model.set=c("NEC3param", "NEC4param", "NECHormesis", "NECsigmoidal")}
  if(model.set=="ECx"){model.set=c("ECx4param", "ECxWeibull1", "ECxWeibull2")}#, "ECxLinear"
  if(model.set=="all"){model.set=c("NEC3param", "NEC4param", "NECHormesis", "NECsigmoidal", 
                                   #"ECxLinear",
                                   "ECx4param", "ECxWeibull1", "ECxWeibull2"
  )}
  
  # if drop.models is not NA 
  if(is.na(drop.models)==F){
    model.set <- model.set[is.na(match(model.set, drop.models))]
  }
  
  # if add.models is not NA 
  if(is.na(add.models)==F){
    model.set <- unique(c(model.set,add.models))
  }
  
  # Fit each of the models
  mod.fits <- vector(mode = 'list', length = length(model.set))
  names(mod.fits) <- model.set
  
  for(m in 1:length(model.set)){
    model <- model.set[m] 
    mod.m <- NULL
    mod.m <- try(jagsMANECfit$mod.fits[[model]], silent=T)
    if(class(mod.m)!="jagsNECfit"){
      fit.m <- try(
      fit.jagsNEC(data=jagsMANECfit$data,
                  x.var=jagsMANECfit$x.var,
                  y.var=jagsMANECfit$y.var,
                  trials.var = jagsMANECfit$trials.var,
                  x.type = jagsMANECfit$x.type, 
                  y.type = jagsMANECfit$y.type,
                  over.disp=jagsMANECfit$over.disp,
                  model=model,
                  added.model=TRUE), 
      silent = TRUE)
    if (!inherits(fit.m, 'try-error')) {
      mod.fits[[model]] <- fit.m  
    } else {
      mod.fits[[model]] <- NA 
    }   
      
      
      
    }else{
      mod.fits[[m]] <- mod.m
    }
    

    
  }
  
  
  # collate all the elements
  export.list <- extract.modStats(mod.fits)
  # assign a class to the output
  class(export.list) <- "jagsMANECfit"
  
  # return the collated output
  return(export.list) 
  
  
  
}