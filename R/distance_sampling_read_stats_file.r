
# Two functions that extracts data from the stats file

split.line <- function(line) {
# Be careful with this function! Parsing might go wrong.
# It still gives some errors,
# I should probably update the function to store the number of the row
# where parsing fails, so that I could go back and check what went wrong.
# Another option is to edit the file in Excel first, safer but time consuming

# Line is an indexing vector, called from the read.stats.file function

# Was changed slightly (substring values differed)
# takes each line and returns the different values for stratum, samp, estimator, ...
  if(nchar(line) < 6) stop ("This isn't a stats file!")
  stratum <- as.integer(substr(line,1,6))
  if(is.na(stratum)) {
    #this means that 
    ok=FALSE
    samp <- NULL
    estimator <- NULL
    module <- NULL
    statistic <- NULL
    value <- NULL  
    cv <- NULL
    lcl <- NULL
    ucl <- NULL
    degrees.freedom <- NULL
  } else {
    ok=TRUE
#       Parse the rest of the line when ok is true
    samp <- as.integer(substr(line, 8, 13))
    estimator <- as.integer(substr(line,14,15))
    module <- as.integer(substr(line, 16,17))
    statistic <- as.integer(substr(line,18,20))
    value <- as.numeric(substr(line, 21,36))
    # Index values were not correct from here:
    cv <- as.numeric(substr(line, 37,52)) # was 37,45
    lcl <- as.numeric(substr(line, 53,67)) # was 46,60
    ucl <- as.numeric(substr(line, 68,82)) # was 61,75
    degrees.freedom <- as.numeric(substr(line, 83,95)) # was 76,91 - also changed from as.integer to as.numeric
   }
  return(list(ok=ok, stratum=stratum, samp=samp,
              estimator=estimator, module=module, statistic=statistic,
              value=value, cv=cv, lcl=lcl, ucl=ucl, degrees.freedom=degrees.freedom))
}
################################################################################

read.stats.file <- function(stat.file.name) {
#   Purpose: extracts results statistics from the MCDS stat file

#Input:
#  stat.file.name - name of file to look in

#Returns list:
#   Dhat.ind    -   density of individuals
#   Dhat.ind.se
#   Dhat.ind.df
#   Es          -   mean cluster size
#   Es.se
#   f0          -   pooled f0
#   f0.se 
#   n           -   number of detected animals
#   n.se
#   L           -   effort
#   nL          -   encounter rate
#   nL.se

#read the file in and test that it has something in it
  lines.v <- readLines(stat.file.name)
  n.lines <- length(lines.v)
  if(n.lines==0) {
    stop ("Nothing to read!")
  }

# go through each line, looking for the results we want and storing them in the appropriate place

  Dhat.ind <- NULL; Dhat.ind.cv <- NULL; Dhat.ind.df <- NULL; Es <- NULL; Es.cv <- NULL; f0 <- NULL;
  f0.cv <- NULL; n <- NULL; n.cv <- NULL; L <- NULL; nL <- NULL; nL.cv <- NULL

  for (line in 1:n.lines) {
    parsed.line <- split.line(lines.v[line])
    if(parsed.line$ok) {
       # Dhat.ind
       if(parsed.line$module==4 & parsed.line$statistic==2){# density of individuals line
          Dhat.ind[parsed.line$stratum] <- parsed.line$value
          Dhat.ind.cv[parsed.line$stratum] <- parsed.line$cv
          Dhat.ind.df[parsed.line$stratum] <- parsed.line$degrees.freedom}
       # Es
       if(parsed.line$module==3 & parsed.line$statistic==4){# size bias adjusted cluster size 
          Es[parsed.line$stratum] <- parsed.line$value
          Es.cv[parsed.line$stratum] <- parsed.line$cv}
       # f0
       if(parsed.line$module==2 & parsed.line$statistic==4){# f0
          f0 <- parsed.line$value
          f0.cv <- parsed.line$cv}
       # n
       if(parsed.line$module==1 & parsed.line$statistic==1){# n
          n[parsed.line$stratum] <- parsed.line$value}
       # L
       if(parsed.line$module==1 & parsed.line$statistic==3){# L
          L[parsed.line$stratum] <- parsed.line$value}
       # nL
       if(parsed.line$module==1 & parsed.line$statistic==4){# nL
          nL[parsed.line$stratum] <- parsed.line$value
          nL.cv[parsed.line$stratum] <- parsed.line$cv}         
    }
  }
  return(list(
          Dhat.ind=Dhat.ind, Dhat.ind.se=Dhat.ind*Dhat.ind.cv, Dhat.ind.df=Dhat.ind.df, 
          Es=Es, Es.se=Es*Es.cv, f0=f0, f0.se=f0*f0.cv,
          n=n, L=L, nL=nL, nL.se=nL*nL.cv
          ))
 }
