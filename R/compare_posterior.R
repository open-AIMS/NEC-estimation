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

#' compare_posterior
#'
#' Compares two posterior samples through differencing. 
#' Can be used to compare two NEC values from jagsNEC or jagsMANEC model fits, 
#' two ECx values also from jagsNEC or jagsMANEC model fits, or simply any two posteriors passed as numeric vectors.
#' 
#' @param x a jagsNEC or jagsMANEC model fit as returned by fit.jagsNEC or fit.jagsMANEC, or a numeric vector.
#' 
#' @param y a jagsNEC or jagsMANEC model fit as returned by fit.jagsNEC or fit.jagsMANEC, or a numeric vector.
#' 
#' @param comparison a character vector indicating the type of comparison to make. Takes values of NEC, ECx or custom.
#' 
#' @param ECx.val the desired percentage effect value. This must be a value between 1 and 99 (for type = "relative" 
#' and "absolute"), defaults to 10.
#' 
#' @param type a character vector indicating the type of ECx.val to calculate, taking values of "relative",  "absolute" (the default) or "direct". 
#' Type "relative" is calculated as the percentage decrease from the maximum predicted value of the response (top) to the minimum predicted value 
#' of the response. Type "absolute" (the default) is calculated as the percentage decrease from the maximum value of the response (top) 
#' to 0 (or bot for a 4 parameter model fit). Type "direct" provides a direct estimate of the x value for a given y.
#' Note that for the current version, ECx for an NECHormesis model is estimated at a percent decline from the control
#' 
#' @import ggplot2 
#' @import tidybayes
#' 
#' @export
#' @return a list containing:
#' \itemize{
#'    \item "post_vals1" the extracted posterior sample for x
#'    \item "post_vals2" the extracted posterior sample for y
#'    \item "diff_posterior" the full differenced posterior sample
#'    \item "Difference" median, lower and upper 95 percent credible bounds of the differenced posterior
#'    \item "Prob_diff" the probability that x is great than y
#'    \item "plot1" a ggplot2 plot of the input posteriors
#'    \item "plot2" a ggplot2 plot of the difference of the input posteriors (x-y)
#' }

#'Probability (%) Posterior 1 is greater than posterior 2

compare_posterior = function(x, y, comparison="NEC", ECx.val=10, ECx.type="absolute"){
  
  
  if(comparison=="NEC"){
    if(class(x)=="jagsNECfit" | class(x)=="jagsMANECfit"){
      post.vals1 = x$sims.list$NEC      
    }else{stop("Your supplied x must be a jagsNEC or jagsMANEC model fit for comparison type NEC")}
    if(class(y)=="jagsNECfit" | class(y)=="jagsMANECfit"){
      post.vals2 = y$sims.list$NEC
    }else{stop("Your supplied y must be a jagsNEC or jagsMANEC model fit for comparison type NEC")}  
  }
  
  if(comparison=="ECx"){
    if(class(x)=="jagsNECfit" | class(x)=="jagsMANECfit"){
      post.vals1 = extract_ECx(x, ECx.val = ECx.val, posterior = TRUE, type=ECx.type)    
    }else{stop("Your supplied x must be a jagsNEC or jagsMANEC model fit for comparison type ECx")}
    if(class(y)=="jagsNECfit" | class(y)=="jagsMANECfit"){
      post.vals2 = extract_ECx(y, ECx.val = ECx.val, posterior = TRUE, type=ECx.type)    
    }else{stop("Your supplied p2 must be a jagsNEC or jagsMANEC model fit for comparison type ECx")}  
  }
  
  if(comparison=="custom"){
    if(class(x)!="numeric"){
      stop("Your supplied x must be a numeric vector for comparison type custom")
    }
    if(class(y)!="numeric"){
      stop("Your supplied p2 must be a numeric vector for comparison type custom")
    }
    post.vals1 <- x
    post.vals2 <- y
    
  }
  

    df1.s = data.frame(PS = post.vals1)
    df1.s$curve = rep('blue', nrow(df1.s))
    med1 = median(df1.s$PS)
    lowerCI = 0.025*length(post.vals1)
    upperCI = 0.975*length(post.vals1)
    medCI = 0.5*length(post.vals1)
  
    df2.s = data.frame(PS = post.vals2)
    df2.s$curve = rep('orange', nrow(df2.s))
    d1.max = max(density(df1.s$PS)$y)
    d2.max = max(density(df2.s$PS)$y)
    ttt = if (d1.max > d2.max){
      (d1.max)
    } else {
      (d2.max)
    }
    
    med2 = median(df2.s$PS)
    
    
    df3.s = rbind(df1.s, df2.s)

    p1 = ggplot2::ggplot(df3.s, aes(x=PS))+geom_density(aes(group=curve, color =curve , fill=curve), alpha=0.3) 
    p1 = p1 + tidybayes::stat_pointintervalh(aes(y = 0.00, x = PS, group=curve),.width = c(.66, .95))+#+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
      theme_light()
    p1 = p1+scale_fill_manual( values = c("steelblue4", "orange"))+
      scale_color_manual( values = c("steelblue4","grey", "steelblue1","steelblue4", "grey","grey", "grey","grey"))+theme(legend.position="none")#nice
    p1 = p1 + scale_y_continuous(name ="Posterior probability density") 
    p1 = p1 + coord_cartesian(ylim = c(0, ttt)) 
    p1 = p1 + scale_x_continuous(name ="Posterior values") 
    p1
    
    #Differences posterior
    #med1 - med2
    df4.s  = data.frame(diff= df1.s$PS - df2.s$PS)
    lower1 = sort(df4.s$diff)[lowerCI]
    upper1 = sort(df4.s$diff)[upperCI]
    med1 = sort(df4.s$diff)[medCI]
    diff.df = data.frame(med1,lower1, upper1)
    d2.max = max(density(df4.s$diff)$y)
    # plot(density(df4.s$diff))
    # density(df4.s$diff)
    p2 = ggplot2::ggplot(df4.s, aes(x=df4.s$diff))+geom_density(aes(x=df4.s$diff, fill = 'steelblue4'), alpha=0.3)+ 
      tidybayes::stat_pointintervalh(aes(y = 0.00, x = diff),.width = c(.66, .95))+#+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
      geom_vline(xintercept = 0, color = "red", lty = 2)+ theme_light()
    p2 = p2+scale_fill_manual( values = c("steelblue4", "orange"))+
      scale_color_manual( values = c("steelblue4","grey", "steelblue1","steelblue4", "grey","grey", "grey","grey"))+theme(legend.position="none")#nice
    p2 = p2 + scale_y_continuous(name ="Differences posterior density") 
    p2 = p2 + coord_cartesian(ylim = c(0.0, d2.max)) 
    p2 = p2 + scale_x_continuous(name ="Standardized effect size") 
    p2
    
    prob.out1.greater = length(which (df1.s$PS > df2.s$PS))/length(post.vals1)*100
    
    list.data = list(post_vals1=post.vals1, post_vals2=post.vals2, 
                     diff_posterior=unlist(df4.s), 
                     Difference=diff.df, 
                     Prob_diff=prob.out1.greater, 
                     plot1=p1, plot2=p2)
    return(list.data)

} 