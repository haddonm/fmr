
# abdat------------------

#' @title abdat fishery data for Block 12 of Tasmanian blacklip fishery
#'
#' @description A dataset containing columns of year, LML, catches, and cpue 
#'     data for use in fmr when exploring size-based integrated models.
#'
#' @format A 58 row x 4 column matrix
#' \describe{
#'   \item{year}{ 9 fishing years from 1963 to 2020}
#'   \item{lml}{ the legal minimum length in the year of catch}
#'   \item{catch}{ the annual catch of blacklip abalone \emph{Haliotis rubra}}
#'   \item{cpue}{ the standardizaed cpue for each year}
#' }
#' 
#' @references McAllister, J. and C. Mundy (2024) \emph{Tasmanian Abalone 
#'     Assessment}, Institute for Marine & Antarctic Studies, University of 
#'     Tamsania p 1-144.
#'     
#' @examples
#' data(abdat)
#' abdat
"abdat"

# fish --------------------

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

# francis92-----------------------

#' @title francis92 Three data objects suitable for use with an ASPM
#'
#' @description A dataset containing the fish data.frame, the glb list, and the
#'     props data.frame set up ready for use with an ASPM. These are data 
#'     obtained from Francis 1992 concerning orange roughy from New Zealand's
#'     Chatham Rise. Francis' 1992 paper is well known as an important 
#'     discussion of risk analysis when providing fisheries management advice.
#'     It was one of the first published uses of an age-structured production 
#'     model. It is even better known as providing a clear exposition of what 
#'     has become a standard re-parameterization of the Beverton-Holt stock-
#'     recruitment model. Francis demonstrates how the B-H model can be 
#'     re-defined using steepness and unfished recruitment level. The full
#'     derivation is provided in Haddon (2021, p172).
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
#' @references Haddon, M. (2021) \emph{Using R for Modelling and Quantitative 
#'     Methods in Fisheries} CRC Press. Chapman & Hall, Boca Raton. 337p. ISBN:
#'     978-0-367-46988-7 see also https://haddonm.github.io/URMQMF/index.html.
#' 
#' @examples
#' data(francis92)
#' str(francis92)
#' 
"francis92"


# fournarch82----------------

#' @title fournarch82 the observed catch-at-age for plaice 
#'
#' @description Table 1 from Fournier and Archibald, 1982, p1200 containing the 
#'     actual numbers-at-age in a simulated population. 
#'
#' @format A 20 row x 10 column  matrix purporting to be the number-at-age in a 
#'    simulated population across 20 years for 10 ages, 4 - 13.
#' 
#' @references Fournier, D. and C.P. Archibald (1982) A general theory for 
#'    analyzing catch at age data. \emph{Canadian Journal of Fisheries and
#'    Aquatic Sciences} \emph{39}: 1195-1207.
#'     
#' @examples
#' data(fournarch82)
#' fournarch82
#' rowfreqboot(x=fournarch82,n=250)
"fournarch82"


# outlike-------------------------------------------------

#' @title outlike a 51 x 51 matrix of likelihoods for LnR0 and Lnq
#'
#' @description A 51 x 51 data.frame with rows as LnR0 = seq(5.9,8.4,0.05) and
#'     columns of Lnq = seq(-8.7,-6.2,0.05). The values are the -veLL for each
#'     combination when combined with a constant LnsigCE = -1.3. This can also
#'     be used with persp to generate a 3D plot of the likelihood surface. 
#'     Note, in teh example, 
#'
#' @format A 51 x 51 data.frame
#' 
#' @examples
#' data(outlike)
#' p1 <- seq(5.9,8.4,0.05)
#' p2 <- -1.3
#' p3 <- seq(-8.7,-6.2,0.05)
#' z <- outlike
#' x <- p1
#' y <- p3
#' trans <- persp(x,y,z,zlim=c(0,20),xlim=c(5.8,8.5),ylim=c(-8.7,-6.0),
#'                theta=-10,phi=50,shade=1e-02,border=NULL,
#'                lwd=0.25,xlab="LnR0",ylab="Lnq",zlab="-veLL",
#'                scale=TRUE,expand=0.4,axes=TRUE,ticktype="detailed")
"outlike"

# orhdat1-----------------------------------

#' @title orhdat1 Three data objects suitable for use with asmreduct
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

# ocaa----------------------------------

#' @title ocaa the observed catch-at-age for plaice 
#'
#' @description Table 13.1 from Beverton & Holt, 1957, p451 containing the relative
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

# owaa ------------------------------------

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

# param ------------------------

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

# pindat ----------------------------

#' @title pindat preliminary parameter estimates for a SBM model for abalone
#' 
#' @description Preliminary parameters used when fitting a statistical catch-at-
#'     size model to the abalone data from McAllister & Mundy (2024). This 
#'     data has also been included in the size-based MSE, 'aMSE' described at
#'     https://haddonm.github.io/aMSEGuide
#'
#' @format A data.frame of two columns
#' \describe{
#'   \item{var}{ the name of each parameter }
#'   \item{value}{ the natural logarithm transformed value of each parameter }
#' }
#' 
#' @references McAllister, J. and C. Mundy (2024) \emph{Tasmanian Abalone 
#'     Assessment}, Institute for Marine & Antarctic Studies, University of 
#'     Tamsania p 1-144.
#' 
#' @examples 
#'  data(plaice)
#'  str(plaice)
#'  print(plaice$fish)
#'  print(plaice$agedata)
"pindat"


# plaice ----------------------------

#' @title plaice data derived from Beverton and Holt, 1957 for European Plaice
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

# robustresults-------------------

#' @title robustresults are 'robustASPM' outputs western Tasmanian Orange Roughy
#'
#' @description robustresults a data.frame containing the results of running
#'     the robustASPM function using the westroughy data sets and using the 
#'     dynamicsF function that has instantaneous fishing mortality rates. It
#'     constitutes a 200 x 14 data.frame with columns iLnR0, isigmaCE, iavq, 
#'     -veLL, LnR0, LsigCE, Lavq, R0, sigCE, avq, MSY, B0, pardist, and Iters.
#'     See the description in teh Introduction to Age-Structured Models 
#'     chapter in More Fisheries Modelling using R.
#'     
#' @format A data.frame of 14 columns and 200 rows  
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
#' @examples
#' \dontrun{
#'  data(robustresults)
#'  str(robustresults)
#' }
"robustresults"

# sizecomp -----------------------

#' @title sizecomp a 38 x 20 matrix of sizecomposition data for abalone
#'
#' @description sizecomp a a 38 x 20 matrix of 2mm size-classes from 136 - 210
#'     mm and 20 years from 1993 - 2020, with counts in each size-class in each
#'     year. Actual sample sizes vary from 345 - 9960, although effective
#'     sample size, of course, is lower than these limits. 
#'
#' @format A 38 x 20 matrix of integer counts
#' \describe{
#'   \item{rows}{ 38 2 mm size-classes from 136 - 210 mm}
#'   \item{columns}{ 20 years from 1993 - 2020}
#' }
#' 
#' @references Haddon, M. (2025) aMSE: A Software Framework for Abalone 
#'     Management Strategy Evaluation_. R package version 1.0.5,
#'     <https://github.com/haddonm/aMSE>.
#'     
#' @examples
#' data(sizecomp)
#' sizecomp
"sizecomp"


# westroughy-------------------

#' @title westroughy A list of 3 data objects for western Tasmanian Orange Roughy
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


