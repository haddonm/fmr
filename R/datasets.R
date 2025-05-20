

#' @title fish fishery data for the plaice example from Beverton and Holt, 1957
#'
#' @description A dataset containing the fish data.frame, for use in fmr
#'
#' @format A data.frame with four columns
#' \describe{
#'   \item{year}{ 9 fishing years from 1929 to 1937, listed as 29 - 37}
#'   \item{effort}{ the annual effort}
#'   \item{catch}{ the annual catch of plaice}
#'   \item{obsce}{ the cpue calculated by dividing the catch by the effort}
#' }
#' 
#' @references Beverton, R.J.H. and S.J. Holt (1957) \emph{On the dynamics of 
#'     exploited fish populations.} U.K. Ministry of Agriculture and Fisheries, 
#'     Fisheries Investigations (Series 2), \emph{19}: 1-533.
#'     
#' @examples
#' data(fish)
#' fish
"fish"

#' @title francis92 Three data objects suitable for use with an ASPM.
#'
#' @description A dataset containing the fish data.frame, the glb list, and the
#'     props data.frame set up ready for use with an ASPM. These are data 
#'     obtained from Francis 1992 concerning orange roughy from New Zealand's
#'     Chatham Rise.
#'
#' @format A list of three objects
#' \describe{
#'   \item{fish}{ a data.frame containing year, catch (tonnes), index (trawl
#'       survey index), and Coefficient of Variation.}
#'   \item{glb}{ a list of global variables including maxage, M, parameters for
#'       growth, weight-at-age, maturity-at-age, steepness, R0, selectivity,
#'       resilience, number of ages, ages, number of years of catch data,
#'       and the name of the species. }
#'   \item{props}{ a data.frame of age, laa, waa, maa, and sela, that is the
#'       length-, weight, maturity, and selectivity-at-age}
#' }
#' 
#' @references Francis, R.I.C.C. (1992) Use of Risk Analysis to Assess Fishery 
#'     Management Strategies: A Case Study using Orange Roughy 
#'     (\emph{Hoplostethus atlanticus}) on the Chatham Rise, New Zealand.
#'     \emph{Canadian Journal of Fisheries and Aquatic Sciences} \emph{49}:
#'     922-930.
#' 
#' @examples
#' data(francis92)
#' str(francis92)
#' 
"francis92"

#
#' @title orhdat1 Three data objects suitable for use with asmreduct.
#'
#' @description A dataset containing a fish data.frame, the glb list, and 
#'     the props data.frame set up ready for use with asmreduct. 
#'
#' @format A list of three objects
#' \describe{
#'   \item{fish}{ a data.frame containing year, catch}
#'   \item{glb}{ a list of global variables including maxage, M, parameters 
#'       for growth, weight-at-age, maturity-at-age, steepness, R0, 
#'       selectivity, resilience, number of ages, the ages themselves, the
#'       number of years of catch data, and the species name}
#'   \item{props}{ a data.frame of age, laa, waa, maa, sela, and feca}
#' }
"orhdat1"

#' @title ocaa the observed catch-at-age for plaice. 
#'
#' @description Table 13.1 from Beverton & Holt, 1957 containing the relative
#'     catch-at-age for plaice from the north sea. This is a data.frame of nine 
#'     years (29/30 - 37/38) by 9 ages (2 - 10). 
#'
#' @format A 9 x 9 data.frame
#' 
#' @references Beverton, R.J.H. and S.J. Holt (1957) \emph{On the dynamics of 
#'     exploited fish populations.} U.K. Ministry of Agriculture and Fisheries, 
#'     Fisheries Investigations (Series 2), \emph{19}: 1-533.
#'     
#' @examples
#' data(ocaa)
#' ocaa
"ocaa"

#' @title owaa the observed weight-at-age for plaice. 
#'
#' @description Table 16.2 from Beverton & Holt, 1957 containing the relative
#'     weight-at-age for plaice from the north sea. This is a data.frame of nine 
#'     years (29/30 - 37/38) by 9 ages (2 - 10). Older ages are omitted.
#'
#' @format A 9 x 9 data.frame
#' 
#' @references Beverton, R.J.H. and S.J. Holt (1957) \emph{On the dynamics of 
#'     exploited fish populations.} U.K. Ministry of Agriculture and Fisheries, 
#'     Fisheries Investigations (Series 2), \emph{19}: 1-533.
#'     
#' @examples
#' data(owaa)
#' owaa
"owaa"

#' @title param preliminary parameter estimates for a SCAA model of plaice
#' 
#' @description Preliminary parameters used when fitting a statistical catch-at-
#'     age model to the plaice data from Beverton & Holt, 1957. Much of this 
#'     data has also been included in the age-structured model described in 
#'     Haddon, 2011.
#'
#' @format A data.frame of two columns
#' \describe{
#'   \item{var}{ the name of each parameter }
#'   \item{value}{ the natural logarithm transformed value of each parameter }
#' }
#' 
#' @references Beverton, R.J.H. and S.J. Holt (1957) \emph{On the dynamics of 
#'     exploited fish populations.} U.K. Ministry of Agriculture and Fisheries, 
#'     Fisheries Investigations (Series 2), \emph{19}: 1-533.
#'     Haddon, M (2011) \emph{Modelling and Quantitative Methods in Fisheries} 
#'     2ns edition. Chapman & Hall, CRC Press. Boca Raton. 449p.
#' 
#' @examples 
#'  data(plaice)
#'  str(plaice)
#'  print(plaice$fish)
#'  print(plaice$agedata)
"param"

#' @title plaice data derived from Beverton and Holt, 1957 for European Plaice.
#' 
#' @description plaice data including fish, glb, props, agedata, and lendata
#'     for North sea plaice derived from tables and the text of the classical
#'     Beverton and Holt, 1957, book. Includes age data that is useful for 
#'     illustrating the catch curves. Much of this data has also been included
#'     in the age-structured model described in Haddon, 2011.
#'
#' @format A list of five objects with only the first four containing data, the
#'     lendata only contains formatted data for illustrating that format, it is
#'     not real data. The other objects contain real data.
#' \describe{
#'   \item{fish}{ a data.frame containing year, catch, cpue, SE of the cpue }
#'   \item{glb}{biological parameters relating to growth, selectivity, 
#'     weight-at-age, steepness, and resilience and spsname }
#'   \item{props}{ contains six variables ages, laa, waa, maa, sela, and feca,
#'     which are all relative to age.}
#'   \item{agedata}{ a list of 5 objects, yrage - the years in which age data are
#'     available, ages - the observed ages, agemax - the maximum age, nage - 
#'     the number of observed ages, and naa - the numbers-at-age by year}
#'   \item{lendata}{ a list of 5 objects akin to the agedata object but for
#'     length data.}
#' }
#' 
#' @references Beverton, R.J.H. and S.J. Holt (1957) \emph{On the dynamics of 
#'     exploited fish populations.} U.K. Ministry of Agriculture and Fisheries, 
#'     Fisheries Investigations (Series 2), \emph{19}: 1-533.
#' 
#' @examples 
#'  data(plaice)
#'  str(plaice)
#'  print(plaice$fish)
#'  print(plaice$agedata)
"plaice"

#' @title westroughy Three data objects describing western Tasmanian Orange Roughy.
#'
#' @description westroughy a dataset containing the fish data.frame, the glb 
#'     list, and the props data.frame set up ready for analysis using an 
#'     age-structured production model. The analysis of the CPUE, not typically 
#'     available for orange roughy fisheries is described in the reference 
#'     below.
#'
#' @format A list of three objects
#' \describe{
#'   \item{fish}{ a data.frame containing Year, Catch, standardized CPUE, and 
#'       SE, the standard error of the CPUE estimates, if present}
#'   \item{glb}{ a list of global variables including maxage, M, parameters for
#'       growth, weight-at-age, maturity-at-age, steepness, R0, selectivity,
#'       resilience, number of ages, and the ages themselves. }
#'   \item{props}{ a data.frame of age, laa, waa, maa, sela, and feca}
#' }
#' 
#' @references Haddon, M. (2018) Western Orange Roughy. pp806-821 \emph{in} 
#'     Tuck, G.N. (ed.) \emph{Stock Assessment for the Southern and Eastern 
#'     Scalefish and Shark Fishery 2016 and 2017. Part 2, 2017}. Australian 
#'     Fisheries Management Authority and CSIRO Oceans and Atmosphere, Hobart.
#'     837p.
#'     
#' @examples
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' str(glb)
"westroughy"
