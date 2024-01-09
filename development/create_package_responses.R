# https://r-pkgs.org

library(devtools)
library(usethis)

p <- "W:/projects/R/responses"
#usethis::create_package(p, check_name = FALSE)

devtools::load_all()

# Must run document() to add export functions to NAMESPACE
devtools::document()
devtools::install()

devtools::test()

# Document data:
# https://r-pkgs.org/data.html

usethis::use_mit_license()

use_git_config(user.name = "peterhellstrom", user.email = "peter.hellstrom@nrm.se")
usethis::use_git()
usethis::use_github()
# GitHub API error (401): Bad credentials

usethis::create_github_token()

use_readme_rmd()
build_readme()

devtools::load_all()

# Must run document() to add export functions to NAMESPACE
devtools::document()
devtools::install()

devtools::test()

# Ignore ----
use_build_ignore(c("backup", "data-raw", "development", "examples"))

## Data sets ----
usethis::use_data_raw()

# Document data:
# https://r-pkgs.org/data.html

install_github("peterhellstrom/responses")


## Load package ----
library(eagles)

## Data sets ----
storrutor
ekorutor
fastighetsblad
wms_layers_data
tms_layers_data

st_as_sf(storrutor, coords = c("easting", "northing"), crs = 3021)
