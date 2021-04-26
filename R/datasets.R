


#' @title dataspm Three data objects suitable for use with datalowSA.
#'
#' @description A dataset containing the fish data.frame, the glb list, and the
#'     props data.frame set up ready for use with datalowSA. In particular it can
#'     be used with the SPM functions, as well as the ASPM functions.
#'
#' @format A list of three objects
#' \describe{
#'   \item{fish}{ a data.frame containing Year, Catch, CPUE, SE, Records, and
#'       GeoM which is the unstandardized geometric mean CPUE }
#'   \item{glb}{ a list of global variables including maxage, M, parameters for
#'       growth, weight-at-age, maturity-at-age, steepness, R0, selectivity,
#'       resilience, number of ages, and the ages themselves. }
#'   \item{props}{ a data.frame of age, laa, waa, maa, sela, and feca}
#' }
"dataspm"

#' @title fishdat Three data objects suitable for use with datalowSA.
#'
#' @description A dataset containing the fish data.frame, the glb list, and the
#'     props data.frame set up ready for use with datalowSA. In particular it can
#'     be used with fitASPM, fitSPM, run_cMSY, and DBSRA.
#'
#' @format A list of three objects
#' \describe{
#'   \item{fish}{ a data.frame containing Year, Catch, CPUE, and SE, the standard
#'       error of the CPUE estimates, if present}
#'   \item{glb}{ a list of global variables including maxage, M, parameters for
#'       growth, weight-at-age, maturity-at-age, steepness, R0, selectivity,
#'       resilience, number of ages, and the ages themselves. }
#'   \item{props}{ a data.frame of age, laa, waa, maa, sela, and feca}
#' }
"fishdat"

#' @title invert data derived from a trawl caught invertebrate fishery.
#'
#' @description A dataset containing the fish data.frame as a 31 x 7 matrix, 
#'     the glb and props data.frames are set to NULL. The fish data.frame has
#'     both the standardized cpue as well as the unstandardized geom, that is
#'     the geometric mean cpue.  This is particularly set up to
#'     be used with the SPM functions but also the Catch-MSY routines.
#'
#' @format A list of three objects only two of which contains data
#' \describe{
#'   \item{fish}{ a data.frame containing year, catch, cpue, SE of the cpue, 
#'       geom, which is the unstandardized geometric mean CPUE, vessel, which
#'       is the number of active vessels reporting catches, and records, which is
#'       the number of cpue records reported each year }
#'   \item{glb}{ contains the resilience and spsname }
#'   \item{props}{ set to NULL}
#' }
#' @examples 
#'  \dontrun{
#'  data(invert)
#'  str(invert)
#'  print(invert$fish)
#' }
"invert"

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

#' @title plaice data derived from Beverton and Holt, 1957 for European Plaice.
#' 
#' @description plaice data including fish, glb, props, agedata, and lendata
#'     for North sea plaice dervied from tables and the text of the classical
#'     Beverton and Holt, 1957, book. Includes age data that is useful for 
#'     illustratung the catch curves. Much of this data has also been included
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
#' @examples 
#'  \dontrun{
#'  data(plaice)
#'  str(plaice)
#'  print(plaice$fish)
#'  print(plaice$agedata)
#' }
"plaice"
