##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  packageStartupMessage(get(".mmes_gpu_message", mode="function")(
    get(".mmes_gpu_status", mode="function")()))
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(magenta(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]  Solving Mixed Model Equations in R (sommer) ", desc$Version," (", desc$Date, ") []",sep="")),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]  ------------- Multivariate Linear Mixed Models --------------   []")),appendLF=TRUE)
    packageStartupMessage(paste0(magenta("[]  Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen
                                                                                           (white(" ")), bgWhite(magenta("M")), bgRed(white(" ")),"  ", bgRed(bold(yellow(" (") )),bgRed(bold(white("W"))), bgRed(bold(yellow(") "))) ) ,"                []")),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Published: PLoS ONE 2016, 11(6):1-15                            []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Dedicated to the University of Chapingo and UW-Madison          []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Type 'vignette('sommer.qg')' for a short tutorial               []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Type 'citation('sommer')' to know how to cite sommer            []"),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(magenta("sommer is updated on CRAN every 3-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(magenta("Current source is available at https://github.com/covaruber/sommer"),appendLF=TRUE)
    packageStartupMessage(magenta("If needed, install as: remotes::install_github('covaruber/sommer')"),appendLF=TRUE)
    packageStartupMessage(magenta("Visit https://covaruber.github.io/sommer/ for more information"),appendLF=TRUE)
  }
  invisible()
}

.onLoad = function(libname, pkgname)
{
  if(requireNamespace("emmeans", quietly=TRUE))
    emmeans::.emm_register("mmes", pkgname)
}
