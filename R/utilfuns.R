







#' @title getvectors exacts a chosen column from data.frames held in a list
#' 
#' @description getvectors assumes one has a list of model outputs where each 
#'     model output contains matrices or data.frames of values. The aim is to 
#'     extract particular columns from particular data.frames from each list
#'     member.
#'
#' @param x a list of model outputs where each member of the list contains at
#'     least one named data.frame or matrix with named columns.
#' @param dfname the name of the selected data.frame within each list member
#' @param var the name of the column from each data.frame to be combined.
#'
#' @returns a matrix of the variable var
#' @export
#'
#' @examples
#' models <- vector(mode="list",length=5)
#' names(models) <- 1:5
#' for (i in 1:5) {
#'   outmod <- matrix(rnorm(25,mean=5,sd=1),nrow=5,ncol=5,
#'                    dimnames=list(1:5,c("a","b","c","d","f")))
#'   models[[i]] <- list(outmod=outmod,index=i)
#' }
#' getvectors(models,dfname="outmod",var="d")
getvectors <- function(x,dfname,var) {  #  x=models; matname="fishery"; var="deplete"
  mats <- lapply(x,"[[",dfname)
  nmod <- length(mats)
  res <- NULL
  for (i in 1:nmod) res <- cbind(res,mats[[i]][,var])
  colnames(res) <- names(x)
  return(res)
} # end of getvectors

#' @title propdiff gives the proportional difference of the range of a vector
#' 
#' @description propdiff takes in either a list of objects all of which 
#'     contain named scalars, or a numeric vector of numbers and calculates 
#'     the range of the values and then outputs the proportional 
#'     difference between the smallest and the largest. So the result is the
#'     proportional difference between the largest and smallest divided by the 
#'     largest. This is useful when conducting a jitter analysis of fitting a
#'     model, each of which uses a set of likelihood components.
#'
#' @param x a list containing multiple named scalar values or a numeric vector
#' @param var the name of the scalar to be extracted from each member of the 
#'     list and converted into a vector of values. If the input is a vector,
#'     this argument can be ''.
#' @param giverge default = TRUE, which means the range of values will be 
#'     returned as well as the propdiff. If FALSE, only the propdiff is
#'     returned.
#'
#' @returns either the range of values and the proportional difference between
#'     the largest and smallest, or just the proportional difference of the 
#'     range.
#' @export
#'
#' @examples
#'   x <- vector(mode="list",length=10)
#'   names(x) <- 1:10
#'   for (i in 1:10) x[[i]] <- c(value=rnorm(1,mean=5, sd=1),n=1)
#'   print(x)
#'   propdiff(x,var="value")
#'   x <- rnorm(10,mean=5, sd=1)
#'   propdiff(x,var="",giverge=FALSE)
propdiff <- function(x,var,giverge=TRUE) {  # x=models; var="LL"
  if (inherits(x,"list")) {
    rge <- range(unlist(lapply(x,"[[",var)),na.rm=TRUE)
  } 
  if (inherits(x,c("matrix","data.frame"))) {
    rge <- range(x[,var],na.rm=TRUE)
    }  else {
    rge <- range(x,na.rm=TRUE)
  }
  diff <- abs(rge[1] - rge[2])
  pdiff <- diff/max(rge)
  if (giverge) {
    return(list(rge=rge,diff=diff,propdiff=pdiff))
  } else {
    return(pdiff)
  }
} # end of propdiff



