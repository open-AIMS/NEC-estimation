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

#' @param model.set A character vector containing the of names of model types to be included in the modified fit.
#'
#' @export
#' @return All successully fitted jagsMANECfit model fit. 
modify.jagsMANEC <- function(jagsMANECfit, model.set=NA){
  
  # if the model set is NA 
  if(is.na(model.set)){
   model.set <- names(jagsMANECfit$mod.fits)
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
                  burnin = jagsMANECfit$n.burnin,
                  n.iter = jagsMANECfit$n.iter,
                  over.disp=jagsMANECfit$mod.fits[[1]]$over.disp,
                  model=model), 
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