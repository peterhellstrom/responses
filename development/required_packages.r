# Check that required packages are installed ----

is.installed <- function(mypkg) mypkg %in% installed.packages()[,1]

install.missing <- function(mypkg, 
                            repos = "https://ftp.acc.umu.se/mirror/CRAN/", 
                            dependencies = TRUE) {
  for (i in 1:length(mypkg)) {
    if (is.installed(mypkg[i]) == FALSE) {
      install.packages(
        pkgs = mypkg[i], 
        lib =.Library, 
        repos = repos, 
        dependencies = dependencies)
    }
  }
}

# Exclude glmmadmb from list of packages (it's not available from CRAN For R 3.x)
#install.packages("R2admb")
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")), type="source")

pkgs <- c(
	"bbmle",
	"MASS",
	"emdbook",
	"nlstools",
	"robustbase",
	"nnet",
	"mgcv",
	"boot",
	"simpleboot",
	"lme4",
	"lmerTest",
	"MCMCglmm",
	"xlsx",
	"nlme",
	"coefplot2",
	"pscl",
	"RLRsim",
	"gstat",
	"sp",
	"plotrix",
	"lattice",
	"RODBC",
	"ltm",
	"fitdistrplus",
	"fdrtool",
	"R.utils",
	"gamlss",
	"gamlss.tr",
	"VGAM",
	"glmmADMB")

is.installed(pkgs)
install.missing(pkgs)

# Load required packages ----
# pkgs <- as.list(c(pkgs))
# lapply(pkgs, require, character.only = TRUE)
