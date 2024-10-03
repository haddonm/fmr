



library(microbenchmark)

cppFunction('double matchCC(double f, double M, double cyr, double Byr) {
  return abs((Byr * (1 - exp(-(M + f))) * f/(M + f)) - cyr);
}')


grsearch <- function(f, interval, M, cyr, Byr, tol = 1e-09) {#.Machine$double.eps^0.25) {
  # Define the golden ratio
  golden_ratio <- (sqrt(5) - 1) / 2
  lower=interval[1]
  upper=interval[2]
  # Calculate initial points
  x1 <- upper - golden_ratio * (upper - lower)
  x2 <- lower + golden_ratio * (upper - lower)
  # Evaluate the function at the initial points
  f1 <- f(x1,M, cyr, Byr)
  f2 <- f(x2,M, cyr, Byr)
  # Iteratively narrow the interval
  while (abs(upper - lower) > tol) {
    if (f1 < f2) {
      upper <- x2
      x2 <- x1
      f2 <- f1
      x1 <- upper - golden_ratio * (upper - lower)
      f1 <- f(x1,M, cyr, Byr)
    } else {
      lower <- x1
      x1 <- x2
      f1 <- f2
      x2 <- lower + golden_ratio * (upper - lower)
      f2 <- f(x2,M, cyr, Byr)
    }
  }
  # Return the midpoint of the final interval as the estimated minimum
  return((lower + upper) / 2)
} # end of grsearch




result <- grsearch(matchC, interval=c(0,1.25), M=M,cyr=catch[2],Byr=expB[2])
print(result)


