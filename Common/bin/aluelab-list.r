#!/usr/bin/env Rscript

## ////////////////////////////////////////////////////////////////////////////
##
## aluelab-list.r
##
## - elaborate alucomb.r fit results
##
## ////////////////////////////////////////////////////////////////////////////

require(methods)

args <- commandArgs(TRUE)
if (length(args) > 0) {
  file.name = args[1]
} else {
  file.name = "average_alucomb.rdata"
}

##--- get alucomb results and data
load(file.name)

quant.val = quant
rm(quant)

##--- list a subtset of fitted quantities
display.names = c(
  "Gamma3", "Gamma5", "Gamma110"
  )

##--- list all fitted quantities
display.names = names(quant.val)

show(rbind(cbind(val=quant.val[display.names],
                 err=quant.err[display.names],
                 "S-factor"=sfact3.types[display.names],
                 "err*S"=quant3.err[display.names])))

