# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}




usethis::use_package("dplyr")
usethis::use_package("haven")
usethis::use_package("sf")
usethis::use_package("surveyPrev")
usethis::use_package("SUMMER")
usethis::use_package("stringr")
usethis::use_package("rlang")

usethis::use_package("ggplot2")
usethis::use_package("patchwork")
