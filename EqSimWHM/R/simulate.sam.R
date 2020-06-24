##' Simulate from a sam object 
##' @method simulate sam 
##' @param object sam fitted object as returned from the \code{\link{sam.fit}} function
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param full.data logical, should each inner list contain a full list of data. Defaults to TRUE  
##' @param ... extra arguments
##' @importFrom stats simulate
##' @details simulates data sets from the model fitted and conditioned on the random effects estimated 
##' @return returns a list of lists. The outer list has length \code{nsim}. Each inner list contains simulated values of \code{logF}, \code{logN}, and \code{obs} with dimensions equal to those parameters.
##' @export
simulate.sam<-function(object, nsim=1, seed=NULL, full.data=TRUE, ...){
  if(!is.null(seed)) set.seed(seed)
  est <- unlist(object$pl)
  if(full.data){
    ret <- replicate(nsim, 
                     c(object$data[names(object$data)!="logobs"],#all the old data
                       object$obj$simulate(est)["logobs"])#simulated observations
                     , simplify=FALSE)
    ret<-lapply(ret, function(x){attr(x,"fleetNames") <- attr(object$data,"fleetNames");x})
  }else{
    ret <- replicate(nsim, object$obj$simulate(est), simplify=FALSE)
  }
  ret
}
