## ////////////////////////////////////////////////////////////////////////////
##
## utility functions for HFAG-tau
##

require(methods, quietly=TRUE)
require(stringr, quietly=TRUE)

## ////////////////////////////////////////////////////////////////////////////
##
## object for doing string translations
##
TrStr = setRefClass("TrStr",
  fields = list(
    .table = "numeric"
    ),
  methods=list(
    initialize = function(set1=character(), set2=character()) {
      ## callSuper(...)
      table = 1:256
      if (length(set1) != 0 || length(set2) != 0) {
        table[as.integer(charToRaw(set1))+1] = as.integer(charToRaw(set2))+1
      }
      .self$.table = table
      .self
    })
  )

##--- translate a string using stored table
rc = TrStr$methods(
  tr = function(str) {
    rawToChar(as.raw(.table[as.integer(charToRaw(str))+1]-1))
  })

##--- translate a string using stored table, delimit numbers by "N"
rc = TrStr$methods(
  trN = function(str) {
    str = gsub("(\\d+)", "N\\1N", str)
    rawToChar(as.raw(.table[as.integer(charToRaw(str))+1]-1))
  })

##--- specialized string translation class for conversion of alphanumeric strings to TeX keywords
TrStr.num2tex = setRefClass("TrStr.num2tex", contains="TrStr",
  methods=list(
    initialize = function() {
      callSuper("0123456789_", "zothfvsneiU")
    })
  )

##
## round data frame to print it in a reproductible way
##
round.data.frame <- function(x, digits = 0) {
  x.names = rownames(x)
  rc = data.frame(lapply(x, function(y) if(is.numeric(y)) round(y, digits) else y))
  rownames(rc) = x.names
  rc
}
  

##
## return numeric id for sorting labels like "Gamma5", "Gamma3by5"
## <n>by<m> are sorted after <n> in ascending order ov <m>
##
alurep.gamma.num.id = function(gamma.name) {
  gamma.name = sub("Unitarity", "Gamma1000", gamma.name, fixed=TRUE)
  gamma.name = sub("GammaAll", "Gamma999", gamma.name, fixed=TRUE)
  num1 = as.numeric(str_match(gamma.name, ("^(\\D+)(\\d+[.]?\\d*)"))[,3])
  tmp = str_match(gamma.name, ("^\\D+\\d+(by|)(\\d*)"))[,3]
  num2 = ifelse(tmp=="", 0, as.numeric(tmp))
  return(1000*num1 + num2)
}

##
## reporting functions
##

##
## return tex label such as \Gamma_1 or \frac{\Gamma_1}{\Gamma_2}
##
alurep.gamma.texlabel = function(str) {
  str = gsub("GammaAll", "\\\\Gamma_{\\\\text{All}}", str, perl=TRUE)
  str = gsub("Gamma(\\d+)by(\\d+)", "\\\\frac{\\\\Gamma_{\\1}}{\\\\Gamma_{\\2}}", str, perl=TRUE)
  str = gsub("Gamma(\\d+)", "\\\\Gamma_{\\1}", str, perl=TRUE)
  return(str)
}

##
## get tex description of a quantity in alucomb2 structure
##
alurep.get.texdescr.nv = function(descr, texdescr) {
  repeat {
    if (!is.null(texdescr) && texdescr != "") {
      ##+++ check for ratio of BR
      texdescr = paste("\\BRF{\\tau^-}{", texdescr, "}", sep="")
      texdescr = gsub("mathrm", "text", texdescr, fixed=TRUE)
      break
    }
    if (is.null(descr) || descr == "") {
      texdescr = ""
      break
    }

    ##--- transform text description into TeX
    descr = sub("\\s+/\\s+G\\(total\\)\\s*", "", descr)
    descr = gsub("ex[.]\\s+K\\(S\\)0 --> pi- pi[+]", "ex. K0", descr, perl=TRUE)
    descr = gsub("\\s*\\(``\\d-prong''\\)", "", descr, perl=TRUE)
    descr = gsub("-->", "\\to", descr, fixed=TRUE)
    descr = gsub("([+-])", "^\\1", descr)
    descr = gsub("\\^-prong", "\\\\text{-prong}", descr)

    ##--- convert <decay mode> into \BRF{\tau}{<decay mode>}
    descr = gsub("G(\\(((?:[^()]++|(?1))+)*+\\))", "\\\\BRF{tau^-}{\\2}", descr, perl=TRUE)

    ##--- transform text description into TeX
    descr = gsub(">=", "\\\\ge{}", descr)
    descr = gsub("nubar\\((mu|tau|e)\\)", "\\\\bar{nu}_\\1", descr)
    descr = gsub("nu\\((mu|tau|e)\\)", "nu_\\1", descr)
    descr = gsub("pi0", "pi^0", descr)
    descr = gsub("Kstar", "K^*", descr)
    descr = gsub("Kbar0", "\\\\bar{K}^0", descr)
    descr = gsub("K0", "K^0", descr)
    descr = gsub("K\\(S\\)0", "K_S^0", descr)
    descr = gsub("K\\(L\\)0", "K_L^0", descr)
    descr = gsub("(nu|tau|mu|pi|omega|eta|gamma)", "\\\\\\1", descr)
    descr = gsub("\\s+\\(ex[.]", "\\\\;(\\\\text{ex.}", descr)
    descr = gsub("(neutrals)", "\\\\text{\\1}", descr)
    descr = gsub("(particles|strange|total)", "\\\\text{\\1}", descr)

    ##--- transform division of BR into a TeX fraction
    descr = gsub("^(.*\\S)\\s*/\\s*(\\S.*)$", "\\\\frac{\\1}{\\2}", descr, perl=TRUE)

    texdescr = descr
    break
  }

  return(texdescr)
}
alurep.get.texdescr = Vectorize(alurep.get.texdescr.nv, USE.NAMES=TRUE)

##
## get tex description of a quantity in alucomb2 structure
##
alurep.tex.quant.descr = function(quant) {
  alurep.get.texdescr.nv(quant$descr, quant$texdescr)
}



##--- return latex command def with specified multi-line body
alurep.tex.cmd = function(cmd, body) {
  paste("\\newcommand{\\", cmd, "}{%\n", body, "}%\n", sep="")
}

##--- return latex command def with specified one-line body
alurep.tex.cmd.short = function(cmd, body) {
  paste("\\newcommand{\\", cmd, "}{", body, "\\xspace}%\n", sep="")
}

##
## return measurement value, stat, syst original values
## all properly formatted for latex printing, separately for
## - quant = val +- stat +- syst
## - val
## - stat
## - syst
##
alurep.tex.meas.val.card.fields = function(meas) {
  if (attr(meas$stat, "input") != "") {
    stat.txt = paste("\\pm", attr(meas$stat, "input"))
  } else {
    stat.txt = paste("{}^{", attr(meas$stat.p, "input"), "}_{", attr(meas$stat.n, "input"), "}", sep="")
  }
  if (attr(meas$syst, "input") != "") {
    syst.txt = paste("\\pm", attr(meas$syst, "input"))
  } else {
    syst.txt = paste("{}^{", attr(meas$syst.p, "input"), "}_{", attr(meas$syst.n, "input"), "}", sep="")
  }
  val.txt = attr(meas$value.orig, "input")

  if (syst.txt == "\\pm 0") {
    quant.tex = c(val.txt, stat.txt)
  } else {
    quant.tex = c(val.txt, stat.txt, syst.txt)
  }

  order = gsub("^[^e]+(|e[+]?([-])?0*(\\d+))$", "\\2\\3", quant.tex, perl=TRUE)
  order = ifelse(order == "", 0, as.numeric(order))
  
  quant.tex = paste(quant.tex, collapse=" ")
  if (max(order) == min(order)) {
    if (max(order) != 0) {
      quant.tex = gsub("(\\S+)e[+-]*\\d+", "\\1", quant.tex, perl=TRUE, ignore.case=TRUE)
      quant.tex = paste("(", quant.tex, ") \\cdot 10^{", max(order), "}")
    }
  }
  
  val.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", val.txt, ignore.case=TRUE)
  stat.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", stat.txt, ignore.case=TRUE)
  syst.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", syst.txt, ignore.case=TRUE)
  quant.tex = gsub("e[+]?([-])?0*(\\d+)", "\\\\cdot 10^{\\1\\2}", quant.tex, ignore.case=TRUE)

  return(list(quant=quant.tex, val=val.tex, stat=stat.tex, syst=syst.tex))
}

##
## return value +- stat +- syst original values
## convert exponential format for latex
##
alurep.tex.meas.val.card = function(meas) {
  rc = alurep.tex.meas.val.card.fields(meas)
  return(paste("\\ensuremath{", rc$quant, "}", sep=""))
}

##
## for the spec vector of measurement, compute the
## appropriate precision and power-of-ten order
##

alurep.precision.order = function(vals, perc=FALSE, signif=4, signif.min=2) {
  vals = vals[vals != 0]
  vals.str = sprintf("%.4g", vals)
  order.str = gsub("^[^eE]+(|[eE]((-[0-9]+)+|[+]([0-9]+)))$", "\\2",
    sprintf("%.4e", vals), perl=TRUE)
  order = ifelse(order.str == "", 0, as.numeric(order.str))
  order.max = max(order)
  order.min = min(order)  
  precision = signif.min - 1 + 1 * max(signif-signif.min, order.max - order.min)
  ## print(cbind(order.max, order.min, precision))
  if (order.max == -1) {
    if (!perc) {
      order.max = 0
      precision = precision+1
    } else {
      order.max = -2
      precision = precision-1
    }
  } else if (order.max == -3) {
    order.max = -2
    precision = precision+1
  } else if (order.max == 1) {
    order.max = 0
    precision = precision-1
  } else if (order.max == 2) {
    order.max = 0
    precision = precision-2
  }
  else {
  }
  return(list(precision=precision, order=order.max))
}

##
## return numeric value formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.val.prec.ord = function(val, precision, order, width=0, perc=FALSE) {
  val = val/10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f", sep=""), val)
  } else if (perc && order == -2) {
    ##--- will act depending on  alurep.precision.order()
    rc = sprintf(paste("%", width, ".", precision, "f\\%%", sep=""), val)
  } else {
    rc = sprintf(paste("%", width, ".", precision, "f\\cdot 10^{%d}", sep=""), val, order)
  }
  return(rc)
}

##
## return numeric value formatted val in a string
## use automatic precision and order of magnitude
##
alurep.tex.val.auto = function(val, width=0, perc=FALSE) {
  rc = alurep.precision.order(val, perc=perc)
  precision = rc$precision
  order = rc$order
  return(alurep.tex.val.prec.ord(val, precision, order, width=width, perc=perc))
}

##
## return quantity formatted val +- stat in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.val.err.prec.ord = function(val, err, precision, order, width=0, perc=FALSE) {
  val = val/10^order
  err = err/10^order
  if (order == 0) {
    rc = sprintf(paste("%", width, ".", precision, "f \\pm %", width, ".", precision, "f", sep=""),
      val, err)
  } else if (perc && order == -2) {
    rc = sprintf(paste("(%", width, ".", precision, "f \\pm %", width, ".", precision, "f)\\%%", sep=""),
      val, err)
  } else {
    rc = sprintf(paste("(%", width, ".", precision, "f \\pm %", width, ".", precision, "f) \\cdot 10^{%d}", sep=""),
      val, err, order)
  }
  return(paste("\\ensuremath{", rc, "}", sep=""))
}

##
## return quantity formatted val +- stat in a string
## according to the self-determined optimal precision and power-of-ten order
##
alurep.tex.val.err.auto = function(val, err, width=0, perc=FALSE) {
  rc = alurep.precision.order(c(val, err), perc=perc)
  precision = rc$precision
  order = rc$order
  return(alurep.tex.val.err.prec.ord(val, err, precision, order, width=width, perc=perc))
}

##
## return measurement formatted val +- stat +- syst in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.meas.val = function(meas, precision, order, width=0, perc=FALSE) {
  norm = 10^order
  str.val = sprintf(paste("%", width, ".", precision, "f", sep=""), meas$value.orig/norm)
  if (meas$stat.p == -meas$stat.n) {
    str.stat = sprintf(paste("\\pm %", width, ".", precision, "f", sep=""), meas$stat/norm)
  } else {
    str.stat = sprintf(paste("{}^{%+", width, ".", precision, "f}_{%+", width, ".", precision, "f}", sep=""),
      meas$stat.p/norm, meas$stat.n/norm)
  }
  if (meas$syst.p == -meas$syst.n) {
    str.syst = sprintf(paste("\\pm %", width, ".", precision, "f", sep=""), meas$syst.orig/norm)
  } else {
    str.syst = sprintf(paste("{}^{%+", width, ".", precision, "f}_{%+", width, ".", precision, "f}", sep=""),
      meas$syst.p/norm, meas$syst.n/norm)
  }
  str.meas = paste(str.val, str.stat, str.syst)
  if (order == 0) {
    rc = paste("", str.meas, "", sep="")
  } else {
    rc = sprintf(paste("(", str.meas, ") \\cdot 10^{%d} ", sep=""), order)
  }
  return(paste("\\ensuremath{", rc, "}", sep=""))
}

##
## for the spec. quantity return the associated measurements
##
alurep.meas.quant = function(quant.name, delta) {
  delta.names = delta[, quant.name]
  meas.names = delta.names[delta.names!=0]
  return(names(meas.names))
}

##
## compute appropriate precision and order of magnitude for printing a measurement
##
alurep.precision.order.meas = function(meas, perc=FALSE) {
  vals = c(meas$value.orig, meas$stat.p, meas$stat.n, meas$syst.p, meas$syst.n)
  return(alurep.precision.order(vals, perc=perc))
}

##
## return measurement formatted val +- stat +- syst in a string
## according to the specified precision and power-of-ten order
##
alurep.tex.meas.val.auto = function(meas, width=0, perc=FALSE) {
  rc = alurep.precision.order.meas(meas, perc=perc)
  precision = rc$precision
  order = rc$order
  return(alurep.tex.meas.val(meas, precision, order, width, perc=perc))
}

##
## get all measurements related to a quantity
## compute appropriate precision and order of magnitude
##
alurep.precision.order.quant = function(quant.name, perc=FALSE, with.meas=FALSE) {
  vals = c(quant.val[quant.name], quant.err[quant.name])
  if (with.meas) {
    meas.names = alurep.meas.quant(quant.name, delta)
    meas.vals = unlist(lapply(measurements[meas.names], function(m) c(m$value.orig, m$stat.p, m$stat.n, m$syst.p, m$syst.n)))
    vals = c(vals, meas.vals)
  }
  return(alurep.precision.order(vals, perc=perc))
}

##
## substitute parameter labes with TeX printable names
##
alurep.subst.params = function(str) {
  str = gsub("BR_eta_2gam", "\\Gamma_{\\eta\\to\\gamma\\gamma}", str, fixed=TRUE)
  str = gsub("BR_eta_neutral", "\\Gamma_{\\eta\\to\\text{neutral}}", str, fixed=TRUE)
  str = gsub("BR_eta_3piz", "\\Gamma_{\\eta\\to3\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_eta_pimpippiz", "\\Gamma_{\\eta\\to\\pi^+\\pi^-\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_eta_charged", "\\Gamma_{\\eta\\to\\text{charged}}", str, fixed=TRUE)
  str = gsub("BR_KS_2piz", "\\Gamma_{K_S\\to\\pi^0\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_KS_pimpip", "\\Gamma_{K_S\\to\\pi^+\\pi^-}", str, fixed=TRUE)
  str = gsub("BR_om_pimpippiz", "\\Gamma_{\\omega\\to\\pi^+\\pi^-\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_om_pimpip", "\\Gamma_{\\omega\\to\\pi^+\\pi^-}", str, fixed=TRUE)
  str = gsub("BR_om_pizgamma", "\\Gamma_{\\omega\\to\\pi^0\\gamma}", str, fixed=TRUE)
  str = gsub("BR_phi_KmKp", "\\Gamma_{\\phi\\to K^+K^-}", str, fixed=TRUE)
  str = gsub("BR_phi_KSKL", "\\Gamma_{\\phi\\to K_S K_L}", str, fixed=TRUE)
  str = gsub("BRA_Kz_KS_KET", "\\Gamma_{<K^0|K_S>}", str, fixed=TRUE)
  str = gsub("BRA_Kz_KL_KET", "\\Gamma_{<K^0|K_L>}", str, fixed=TRUE)
  str = gsub("BRA_Kzbar_KS_KET", "\\Gamma_{<\\bar{K}^0|K_S>}", str, fixed=TRUE)
  str = gsub("BRA_Kzbar_KL_KET", "\\Gamma_{<\\bar{K}^0|K_L>}", str, fixed=TRUE)
  str = gsub("BRA_KzKzbar_KLKL_KET_by_BRA_KzKzbar_KSKS_KET", "R_{0\\bar{0}SS/LL}", str, fixed=TRUE)

  str = gsub("BR_f1_2pizpippim", "\\Gamma_{f_1\\to\\pi^+\\pi^-2\\pi^0}", str, fixed=TRUE)
  str = gsub("BR_f1_2pip2pim", "\\Gamma_{f_1\\to2\\pi^+2\\pi^-}", str, fixed=TRUE)
}

##
## get left and right elements of a constraint equation as TeX expressions
##
alurep.tex.constraint.zeroval = function(constr.str) {
  ##--- replace "0 = -x + y" into "x = y"
  rc = str_match(unlist(constr.str), "-([[:alnum:]]+)\\s+[+]\\s+(.*\\S)\\s*$")
  constr.left = ifelse(!is.na(rc[,2]), rc[,2], "")
  constr.right = ifelse(!is.na(rc[,3]), rc[,3], "")
  ##--- remove outer braces
  constr.right = gsub("^\\s*(\\(((?:[^()]++|(?1))+)*+\\))\\s*$", "\\2", constr.right, perl=TRUE)

  ##--- force GammaXbyY constraint to be GammaX/GammaY
  constr.right = ifelse(regexpr("^Gamma(\\d+)by(\\d+)$", constr.left, perl=TRUE) != -1, constr.left, constr.right)

  ##--- convert Gamma# notation for TeX
  constr.left = alurep.gamma.texlabel(constr.left)
  constr.right = alurep.gamma.texlabel(constr.right)
  return(list(left=constr.left, right=constr.right))
}
alurep.tex.constraint.nonzeroval = function(constr.val, constr.str) {
  constr.left = as.character(constr.val)
  constr.right = constr.str
  ##--- convert Gamma# notation for TeX
  constr.right = alurep.gamma.texlabel(constr.right)
  return(list(left=constr.left, right=constr.right))
}
alurep.tex.constraint = function(constr.val, constr.str) {
  sel.zero = which(constr.val == 0)
  sel.nonzero = which(constr.val != 0)
  rc.zero = alurep.tex.constraint.zeroval(constr.str[sel.zero])
  rc.nonzero = alurep.tex.constraint.nonzeroval(constr.val[sel.nonzero], constr.str[sel.nonzero])

  constr.left = character(0)
  constr.right = character(0)
  constr.left[c(sel.zero, sel.nonzero)] = c(rc.zero$left, rc.nonzero$left)
  constr.right[c(sel.zero, sel.nonzero)] = c(rc.zero$right, rc.nonzero$right)
  names(constr.left) = names(constr.val)
  names(constr.right) = names(constr.val)
  
  ##--- split long eqs
  constr.right.split = gsub("(([^+*]+[+*]){6}[^+]+)[+]", "\\1 \\\\\\\\ \n  {}& +", constr.right, perl=TRUE)
  
  constr.right = gsub("*", "\\cdot{}", constr.right, fixed=TRUE)
  constr.right = alurep.subst.params(constr.right)
  constr.right.split = gsub("*", "\\cdot{}", constr.right.split, fixed=TRUE)
  constr.right.split = alurep.subst.params(constr.right.split)
  
  return(list(left = constr.left, right=constr.right, right.split=constr.right.split))
}
