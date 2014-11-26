#!/usr/bin/env Rscript

##
## mkreport.r
##

require("optparse", quietly=TRUE)
require("stringr", quietly=TRUE)
source("../../../Common/bin/alureport.r")
source("../../../Common/bin/alucomb2-fit.r")
source("../../../Common/bin/aluelab3.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## example function to use vars of an env in the function env
##
mkreport.rc = function(alucomb.env) {
  attach(alucomb.env)
  rc = try({
  })
  detach(alucomb.env)
  return(invisible(NULL))
}

##
## get the inspire bibtex key from the measurement tags
##
get.reference = function(tags) {
  bib.pdg.table = list(
    "BARTELT.96" = "Bartelt:1996iv",
    "DEL-AMO-SANCHEZ.11E" = "delAmoSanchez:2010pc",
    "AUBERT.10B" = "Aubert:2009ag",
    "AUBERT.10F" = "Aubert:2009qj",
    "HAYASAKA.10" = "Hayasaka:2010np",
    "LEE.10" = "Lee:2010tc",
    "LEES.10A" = "not found",
    "MIYAZAKI.10" = "Miyazaki:2009wc",
    "MIYAZAKI.10A" = "Miyazaki:2010qb",
    "AUBERT.09AK" = "Aubert:2009ra",
    "AUBERT.09D" = "Aubert:2009ys",
    "AUBERT.09W" = "Aubert:2009ap",
    "GROZIN.09A" = "Grozin:2008nw",
    "INAMI.09" = "Inami:2008ar",
    "MIYAZAKI.09" = "Miyazaki:2008mw",
    "AUBERT.08" = "Aubert:2007mh",
    "AUBERT.08AE" = "Aubert:2008nj",
    "AUBERT.08K" = "Aubert:2007kx",
    "FUJIKAWA.08" = "Fujikawa:2008ma",
    "HAYASAKA.08" = "Hayasaka:2007vc",
    "MIYAZAKI.08" = "Miyazaki:2007zw",
    "NISHIO.08" = "Nishio:2008zx",
    "ANASHIN.07" = "Anashin:2007zz",
    "AUBERT.07AP" = "Aubert:2007jh",
    "AUBERT.07BK" = "Aubert:2007pw",
    "AUBERT.07I" = "Aubert:2006cz",
    "BELOUS.07" = "Abe:2006vf",
    "EPIFANOV.07" = "Epifanov:2007rf",
    "MIYAZAKI.07" = "Miyazaki:2007jp",
    "ABDALLAH.06A" = "Abdallah:2003cw",
    "AUBERT.06C" = "Aubert:2005wa",
    "AUBERT,B.06" = "Aubert:2006hw",
    "INAMI.06" = "Inami:2006vd",
    "MIYAZAKI.06" = "Miyazaki:2005ng",
    "MIYAZAKI.06A" = "Miyazaki:2006sx",
    "YUSA.06" = "Yusa:2006qq",
    "ARMS.05" = "Arms:2005qg",
    "AUBERT,B.05A" = "Aubert:2005ye",
    "AUBERT,B.05F" = "Aubert:2005nx",
    "AUBERT,B.05W" = "Aubert:2005waa",
    "AUBERT,BE.05D" = "Aubert:2005tp",
    "ENARI.05" = "Enari:2005gc",
    "HAYASAKA.05" = "Hayasaka:2005xw",
    "SCHAEL.05C" = "Schael:2005am",
    "ABBIENDI.04J" = "Abbiendi:2004xa",
    "ABDALLAH.04K" = "Abdallah:2003xd",
    "ABDALLAH.04T" = "Abdallah:2003yq",
    "ABE.04B" = "Abe:2003sx",
    "ACHARD.04G" = "Achard:2004jj",
    "AUBERT.04J" = "Aubert:2003pc",
    "ENARI.04" = "Enari:2004ax",
    "YUSA.04" = "Yusa:2004gm",
    "ABBIENDI.03" = "Abbiendi:2002jw",
    "BRIERE.03" = "Briere:2003fr",
    "HEISTER.03F" = "Heister:2002ik",
    "INAMI.03" = "Inami:2002ah",
    "CHEN.02C" = "not found",
    "ABBIENDI.01J" = "Abbiendi:2000ee",
    "ABREU.01M" = "Abreu:2001wc",
    "ACCIARRI.01F" = "Acciarri:2001sg",
    "ACHARD.01D" = "Achard:2001pk",
    "ANASTASSOV.01" = "Anastassov:2000xu",
    "HEISTER.01E" = "Heister:2001me",
    "ABBIENDI.00A" = "Abbiendi:2000kz",
    "ABBIENDI.00C" = "Abbiendi:1999pm",
    "ABBIENDI.00D" = "Abbiendi:1999cq",
    "ABREU.00L" = "Abreu:2000sg",
    "ACCIARRI.00B" = "Acciarri:2000vq",
    "AHMED.00" = "not found",
    "ALBRECHT.00" = "Albrecht:2000yg",
    "ASNER.00" = "Asner:1999kj",
    "ASNER.00B" = "Asner:2000nx",
    "BERGFELD.00" = "Bergfeld:1999yh",
    "BROWDER.00" = "Browder:1999fr",
    "EDWARDS.00A" = "Edwards:1999fj",
    "GONZALEZ-SPRINBERG.00" = "GonzalezSprinberg:2000mk",
    "ABBIENDI.99H" = "Abbiendi:1998cx",
    "ABREU.99X" = "Abreu:1999rb",
    "ACKERSTAFF.99D" = "Ackerstaff:1998yk",
    "ACKERSTAFF.99E" = "Ackerstaff:1998ia",
    "BARATE.99K" = "Barate:1999hi",
    "BARATE.99R" = "Barate:1999hj",
    "BISHAI.99" = "Bishai:1998gf",
    "GODANG.99" = "Godang:1999ge",
    "RICHICHI.99" = "Richichi:1998bc",
    "ACCIARRI.98C" = "Acciarri:1998zc",
    "ACCIARRI.98E" = "Acciarri:1998iv",
    "ACCIARRI.98R" = "Acciarri:1998as",
    "ACKERSTAFF.98M" = "Ackerstaff:1997tx",
    "ACKERSTAFF.98N" = "Ackerstaff:1998mt",
    "ALBRECHT.98" = "Albrecht:1997gn",
    "BARATE.98" = "Barate:1997ma",
    "BARATE.98E" = "Barate:1997tt",
    "BLISS.98" = "Bliss:1997iq",
    "ABE.97O" = "Abe:1997dy",
    "ACKERSTAFF.97J" = "Ackerstaff:1997is",
    "ACKERSTAFF.97L" = "Ackerstaff:1996gy",
    "ACKERSTAFF.97R" = "Ackerstaff:1997dv",
    "ALEXANDER.97F" = "Alexander:1997bv",
    "AMMAR.97B" = "Ammar:1996xh",
    "ANASTASSOV.97" = "Anastassov:1996tc",
    "ANDERSON.97" = "Anderson:1997kb",
    "AVERY.97" = "not found",
    "BARATE.97I" = "Barate:1997hw",
    "BARATE.97R" = "Barate:1997be",
    "BERGFELD.97" = "Bergfeld:1997zt",
    "BONVICINI.97" = "Bonvicini:1997bw",
    "BUSKULIC.97C" = "Buskulic:1996qs",
    "COAN.97" = "Coan:1997am",
    "EDWARDS.97" = "not found",
    "EDWARDS.97B" = "not found",
    "ESCRIBANO.97" = "Escribano:1996wp",
    "ABREU.96B" = "Abreu:1995xt",
    "ACCIARRI.96H" = "Acciarri:1996ht",
    "ACCIARRI.96K" = "Acciarri:1996kc",
    "ALAM.96" = "Alam:1995mt",
    "ALBRECHT.96E" = "Albrecht:1996gr",
    "ALEXANDER.96D" = "Alexander:1995tw",
    "ALEXANDER.96E" = "Alexander:1996nc",
    "ALEXANDER.96S" = "Alexander:1996um",
    "BAI.96" = "Bai:1995hf",
    "BALEST.96" = "Balest:1996cr",
    "BARTELT.96" = "Bartelt:1996iv",
    "BUSKULIC.96" = "Buskulic:1995ty",
    "BUSKULIC.96C" = "Buskulic:1995rh",
    "COAN.96" = "Coan:1996iu",
    "ABE.95Y" = "Abe:1995yj",
    "ABREU.95T" = "Abreu:1995gr",
    "ABREU.95U" = "Abreu:1995gs",
    "ACCIARRI.95" = "Acciarri:1994vr",
    "ACCIARRI.95F" = "Acciarri:1995kx",
    "AKERS.95F" = "Akers:1994xf",
    "AKERS.95I" = "Akers:1995nh",
    "AKERS.95P" = "Akers:1995vy",
    "AKERS.95Y" = "Akers:1995ry",
    "ALBRECHT.95" = "Albrecht:1994nm",
    "ALBRECHT.95C" = "Albrecht:1995un",
    "ALBRECHT.95G" = "Albrecht:1995ht",
    "ALBRECHT.95H" = "Albrecht:1995wx",
    "BALEST.95C" = "Balest:1995kq",
    "BUSKULIC.95C" = "Buskulic:1994yh",
    "BUSKULIC.95D" = "Buskulic:1994hi",
    "ABREU.94K" = "Abreu:1994fi",
    "AKERS.94E" = "Akers:1994ex",
    "AKERS.94G" = "Akers:1994td",
    "ALBRECHT.94E" = "Albrecht:1994es",
    "ARTUSO.94" = "Artuso:1994ii",
    "BARTELT.94" = "Bartelt:1994hn",
    "BATTLE.94" = "Battle:1994by",
    "BAUER.94" = "Bauer:1993wn",
    "BUSKULIC.94D" = "Buskulic:1993mn",
    "BUSKULIC.94E" = "Buskulic:1994je",
    "BUSKULIC.94F" = "Buskulic:1994jf",
    "GIBAUT.94B" = "Gibaut:1994ik",
    "ADRIANI.93M" = "Adriani:1993gk",
    "ALBRECHT.93C" = "Albrecht:1992ka",
    "ALBRECHT.93G" = "Albrecht:1993fr",
    "BALEST.93" = "not found",
    "BEAN.93" = "Bean:1992hh",
    "BORTOLETTO.93" = "Bortoletto:1993px",
    "ESCRIBANO.93" = "Escribano:1993pq",
    "PROCARIO.93" = "Procario:1992hd",
    "ABREU.92N" = "Abreu:1992gn",
    "ACTON.92F" = "Acton:1992ff",
    "ACTON.92H" = "Acton:1992by",
    "AKERIB.92" = "Akerib:1992hb",
    "ALBRECHT.92D" = "Albrecht:1991rh",
    "ALBRECHT.92K" = "Albrecht:1992uba",
    "ALBRECHT.92M" = "Albrecht:1992td",
    "ALBRECHT.92Q" = "Albrecht:1992pe",
    "AMMAR.92" = "Ammar:1991jc",
    "ARTUSO.92" = "Artuso:1992qu",
    "BAI.92" = "Bai:1992bu",
    "BATTLE.92" = "Battle:1992qv",
    "BUSKULIC.92J" = "Buskulic:1992jj",
    "DECAMP.92C" = "Decamp:1991jp",
    "ADEVA.91F" = "Adeva:1991qq",
    "ALBRECHT.91D" = "Albrecht:1990nc",
    "ALEXANDER.91D" = "Alexander:1991am",
    "ANTREASYAN.91" = "Antreasyan:1990fq",
    "GRIFOLS.91" = "Grifols:1990ha",
    "ABACHI.90" = "Abachi:1989vr",
    "ALBRECHT.90E" = "Albrecht:1990zj",
    "BEHREND.90" = "Behrend:1989qx",
    "BOWCOCK.90" = "Bowcock:1989mq",
    "DELAGUILA.90" = "delAguila:1990jg",
    "GOLDBERG.90" = "Goldberg:1990ep",
    "WU.90" = "Wu:1989hx",
    "ABACHI.89B" = "Abachi:1988gx",
    "BEHREND.89B" = "Behrend:1989wc",
    "JANSSEN.89" = "Janssen:1989wg",
    "KLEINWORT.89" = "Kleinwort:1988sb",
    "ADEVA.88" = "Adeva:1988eq",
    "ALBRECHT.88B" = "Albrecht:1987zf",
    "ALBRECHT.88L" = "Albrecht:1988bt",
    "ALBRECHT.88M" = "Albrecht:1988ik",
    "AMIDEI.88" = "Amidei:1987zj",
    "BEHREND.88" = "Behrend:1987fa",
    "BRAUNSCHWEIG.88C" = "Braunschweig:1988yy",
    "KEH.88" = "Keh:1988gs",
    "TSCHIRHART.88" = "Tschirhart:1987uh",
    "ABACHI.87B" = "Abachi:1987hw",
    "ABACHI.87C" = "Abachi:1987nh",
    "ADLER.87B" = "Adler:1987bf",
    "AIHARA.87B" = "Aihara:1986mw",
    "AIHARA.87C" = "Aihara:1987zm",
    "ALBRECHT.87L" = "not found",
    "ALBRECHT.87P" = "Albrecht:1987fb",
    "BAND.87" = "Band:1987nv",
    "BAND.87B" = "Band:1987bm",
    "BARINGER.87" = "Baringer:1987tr",
    "BEBEK.87C" = "Bebek:1987nh",
    "BURCHAT.87" = "Burchat:1986yp",
    "BYLSMA.87" = "Bylsma:1986zy",
    "COFFMAN.87" = "Coffman:1987da",
    "DERRICK.87" = "Derrick:1987sp",
    "FORD.87" = "Ford:1986zk",
    "FORD.87B" = "Ford:1987ha",
    "GAN.87" = "Gan:1987mr",
    "GAN.87B" = "Gan:1987xs",
    "AIHARA.86E" = "Aihara:1986nj",
    "BARTEL.86D" = "Bartel:1986un",
    "RUCKSTUHL.86" = "Ruckstuhl:1986qg",
    "SCHMIDKE.86" = "Schmidke:1986gp",
    "YELTON.86" = "Yelton:1985sx",
    "ALTHOFF.85" = "Althoff:1984ja",
    "ASH.85B" = "Ash:1985hp",
    "BALTRUSAITIS.85" = "Baltrusaitis:1985fh",
    "BARTEL.85F" = "not found",
    "BEHRENDS.85" = "Behrends:1985pm",
    "BELTRAMI.85" = "Beltrami:1985vb",
    "BERGER.85" = "Berger:1985in",
    "BURCHAT.85" = "Burchat:1985mv",
    "FERNANDEZ.85" = "Fernandez:1984if",
    "MILLS.85" = "Mills:1985mh",
    "AIHARA.84C" = "Aihara:1984hr",
    "BEHREND.84" = "Behrend:1984up",
    "MILLS.84" = "Mills:1984dn",
    "BEHREND.83C" = "not found",
    "SILVERMAN.83" = "Silverman:1982ft",
    "BEHREND.82" = "not found",
    "BLOCKER.82B" = "Blocker:1982rw",
    "BLOCKER.82D" = "not found",
    "FELDMAN.82" = "Feldman:1981md",
    "HAYES.82" = "Hayes:1981bn",
    "BERGER.81B" = "not found",
    "DORFAN.81" = "Dorfan:1980mk",
    "BRANDELIK.80" = "not found",
    "BACINO.79B" = "Bacino:1979fz",
    "BACINO.78B" = "Bacino:1978gb",
    "BRANDELIK.78" = "not found",
    "JAROS.78" = "Jaros:1978um",
    "DAVIER.06" = "Davier:2005xq",
    "RAHAL-CALLOT.98" = "RahalCallot:1998ki",
    "GENTILE.96" = "Gentile:1995ue",
    "WEINSTEIN.93" = "not found",
    "PERL.92" = "Perl:1991gd",
    "PICH.90" = "Pich:1990gh",
    "BARISH.88" = "Barish:1987nj",
    "GAN.88" = "Gan:1987zi",
    "HAYES.88" = "Hayes:1988ut",
    "PERL.80" = "not found",
    "LEES.2012X" = "Lees:2012ks",
    "LEES.2012Y" = "Lees:2012de",
    "Belous:2013dba" = "Belous:2013dba",
    "Ryu:2014vpc" = "Ryu:2014vpc"
    )

  bib.conf.table = list(
    "BaBar.ICHEP08" = "Aubert:2008an",
    "BaBar.DPF09" = "Paramesvaran:2009ec",
    "Belle.PHIPSI11" = "SooRyu:2011aa"
    )

  ref.conf = bib.conf.table[[paste(tags[1], tail(tags, -3), collapse=".", sep=".")]]
  if (!is.null(ref.conf)) return(ref.conf)

  tag = paste(tail(tags, -3), collapse=".")
  ref.insp = bib.pdg.table[[tag]]
  if (is.null(ref.insp) || ref.insp=="not found") {
    ref.insp = paste("not found:", paste(tags, collapse="."))
  }
  return(ref.insp)
}

##
## for the spec. measurement return list with:
## - experiment
## - latex \cite{} reference
## - type: prelim or pub
## - formatted val +- stat +- syst
##
get.tex.meas = function(meas, precision, order) {
  meas.item = list()
  meas.item$exp = meas$tags[1]
  meas.item$exp = sub("BaBar", "\\\\babar", meas.item$exp, ignore.case=TRUE)
  meas.item$ref = paste("\\cite{", get.reference(meas$tags), "}", sep="")
  meas.item$type = meas$tags[3]
  meas.item$value = alurep.tex.meas.val(meas, precision, order)
  return(meas.item)
}

##
## latex defs for all measurements
##
get.tex.meas.defs = function() {
  meas.used = measurements[combination$measurements]
  meas.used.names = names(meas.used)
  meas.used.names = meas.used.names[order(meas.used.names)]

  rc = mapply(function(meas.name, meas) {
    exp = meas$tags[1]
    ref = get.reference(meas$tags)
    type = meas$tags[3]
    quant = meas$quant
    rc2 = alurep.tex.meas.val.card.fields(meas)
    paste("\\htmeasdef{", meas.name, "}{",
          quant, "}{", exp, "}{", ref, "}{",
          rc2$quant, "}{",
          rc2$val, "}{", rc2$stat, "}{", rc2$syst, "}%",
          sep="")
  }, meas.used.names, measurements[meas.used.names])
  return(paste(rc, collapse="\n"))
}

##
## return TeX defs for listing a quantity and its measurements
## with the same precision and order
##
## - create TeX defs for the quantiy and all its measurements
## - create a
## return body of latex tabular environment for the requested quantities
## include
## - quantity description, HFAG average
## - list of experimental measurements, with reference
##
get.tex.table = function(quant.names, perc=FALSE) {
  quant.order = order(alurep.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]

  ##--- get TeX defs for quantities and their measurements
  qm.defs = mapply(function(quant.name, quant) {
    ##--- quantity value according to precision / order found for quantity and its measurements
    rc = alurep.precision.order.quant(quant.name, perc=perc, with.meas=TRUE)
    precision = rc$precision
    order = rc$order

    ##--- TeX def for the quantity
    rc = paste(
      "\\htdef{", quant.name, ".qt}{",
      alurep.tex.val.err.prec.ord(quant.val[quant.name], quant.err[quant.name], precision, order, perc=FALSE),
      "}%", sep=""
      )
    
    ##--- TeX defs for the quantity measurements
    meas.names = alurep.meas.quant(quant.name, delta)
    meas.lines = sapply(meas.names, function(meas.name) {
      meas.item = get.tex.meas(measurements[[meas.name]], precision, order)
      rc = paste(
        "\\htdef{", meas.name, ",qt}{", meas.item$value, "}%",
        sep=""
        )
    })
    
    ##--- join all defs into one string
    if (length(meas.lines) > 0) {
      meas.lines = paste(meas.lines, collapse="\n")
      rc = paste(rc, meas.lines, sep=" \n")
    }

    return(rc)
  },
    quant.names, combination$quantities[quant.names])
  qm.defs = paste(qm.defs, collapse=" \n")

  ##--- get TeX defs for a tabular body listing a quantity and its measurements (using previous defs)
  qm.defs2 = mapply(function(quant.name, quant) {
    ##--- TeX expr for quantity
    quant.descr = paste("\\htuse{", quant.name, ".gn}", " = ", "\\htuse{", quant.name, ".td}", sep="")
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}\n", sep="")

    ##--- get quantity value with proper precision and order
    rc = paste(quant.descr, " & \\htuse{", quant.name, ".qt} & \\hfagFitLabel", sep="")

    ##--- get measurements value, experiment and citation
    meas.names = alurep.meas.quant(quant.name, delta)
    meas.lines = sapply(meas.names, function(meas.name) {
      rc = paste(
        paste("\\htuse{", meas.name, ",qt}", sep=""),
        paste("\\htuse{", meas.name, ",exp}", sep=""),
        paste("\\htuse{", meas.name, ",ref}", sep=""),
        sep=" & ")
      return(rc)
    })
    meas.lines = paste(meas.lines, collapse=" \\\\\n")
    
    ##--- define macro for tabular body
    quant.meas.tex = paste("\\htdef{", quant.name, ".qm}{%\n", rc, sep="")
    ##--- add measurements if there are any
    if (nchar(meas.lines) > 0) {
      quant.meas.tex = paste(quant.meas.tex, "\\\\\n", meas.lines, "\n", sep="")
    }
    quant.meas.tex = paste(quant.meas.tex, "}%", sep="")
    return(quant.meas.tex)
  },
    quant.names, combination$quantities[quant.names])
  qm.defs2 = paste(qm.defs2, collapse=" \n")

  ##--- join the two previous blocks of definitions
  qm.defs = paste(qm.defs, qm.defs2, sep="\n")
  
  ##--- build complete tabular body of all quantities and measurements
  qm.table = mapply(function(quant.name, quant) {
    rc = paste("\\htuse{", quant.name, ".qm}", sep="")
    return(rc)
  },
    quant.names, combination$quantities[quant.names])
  
  with.meas = TRUE
  if (with.meas) {
    qm.table = paste(qm.table, collapse=" \\\\\n\\hline\n")
  } else {
    qm.table = paste(qm.table, collapse=" \\\\\n")
  }

  return(list(defs=qm.defs, table=qm.table))
}

##
## return body of latex tabular environment for the requested quantities
## include
## - quantity description and HFAG average
##
get.tex.table.simple = function(quant.names, precision, order) {
  quant.order = order(alurep.gamma.num.id(quant.names))
  quant.names = quant.names[quant.order]
  rc = mapply(function(quant.name, quant) {
    quant.descr = paste("\\htuse{", quant.name, ".gn}", " = ", "\\htuse{", quant.name, ".td}", sep="")
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}\n", sep="")
    rc = paste(
      quant.descr, "&",
       alurep.tex.val.err.prec.ord(quant.val[quant.name], quant.err[quant.name], precision, order)
      )
  }, quant.names, combination$quantities[quant.names])
  return(paste(rc, collapse=" \\\\\n"))
}

##
## return measurements by collaboration
##
get.tex.meas.by.collab = function() {
  toTex = TrStr.num2tex$new()
  collab.meas = sapply(measurements[combination$measurements], function(meas) meas$tags[1])
  collabs = sort(unique(collab.meas))
  collab.nmeas = sapply(collabs, function(collab) sum(collab.meas == collab))
  return(alurep.tex.cmd.short(paste("NumMeas", collabs, sep=""), collab.nmeas))
}

##
## return table body with measurements by reference
##
get.tex.meas.by.ref = function() {
  meas.paper = list()
  rc = mapply(function(meas.name, meas) {
    quant.name = meas$quant

    pub.flag = (regexpr("^pub", meas$tags[3], ignore.case=TRUE) != -1)
    if (pub.flag) {
      pub = "pub"
      cit = paste(meas$tags[-(1:3)], collapse=" ")
    } else {
      pub = ""
      cit = paste(meas$tags[1], "prelim.", meas$tags[-(1:3)], sep=" ")
    }
    expnt = meas$tags[1]
    ref = get.reference(meas$tags)
    
    return(list(cit=cit, pub.flag=pub, expnt=expnt, ref=ref, quant=quant.name, meas=meas.name))
  }, combination$measurements, measurements[combination$measurements])

  ##--- convert to data.frame
  rc = do.call(rbind.data.frame, apply(rc, 2, function(x) {as.data.frame(t(unlist(x)))}))

  ##--- get citations grouped by experiment
  expnt = unique(rc$expnt)
  ##--- sort experiments alphabetically
  expnt = expnt[order(expnt)]
  cits = lapply(expnt, function(x) {
    df = subset(rc, expnt == x & pub.flag == "pub")
    cits = as.character(unique(df$cit))
    cits = cits[order(cits)]
    df = subset(rc, expnt == x & pub.flag != "pub")
    cits.prelim = as.character(unique((df$cit)))
    cits.prelim = cits.prelim[order(cits.prelim)]
    return(c(cits, cits.prelim))
  })
  cits = unlist(cits, recursive = TRUE)
  
  meas.by.ref = lapply(cits, function(x) {
    df = subset(rc, cit == x)
    ##--- reorder measurements according to quantity
    quant.order = order(alurep.gamma.num.id(df$quant))
    df = df[quant.order,]

    expnt = df$expnt[1]
    cite.tex = paste("\\cite{", df$ref[1], "}", sep="")
    if (df$pub.flag[1] == "pub") {
      expnt = sub("^BaBar$", "\\\\babar", expnt, perl=TRUE, ignore.case=TRUE)
      ref.tex = paste(df$cit[1], " (", expnt, ") ", cite.tex, sep="")
    } else {
      cit = df$cit[1]
      cit = sub("^BaBar", "\\\\babar", cit, perl=TRUE, ignore.case=TRUE)
      ref.tex = paste(cit, cite.tex, sep=" ")
    }

    quant.descr = paste("\\htuse{", df$quant, ".gn}", " = ", "\\htuse{", df$quant, ".td}", sep="")
    quant.descr = paste("\\begin{ensuredisplaymath}\n", quant.descr, "\n\\end{ensuredisplaymath}", sep="")

    val.lines = paste(quant.descr, " & ", "\\htuse{", df$meas, "}", sep="")
    val.tex = paste(val.lines, collapse="\n\\\\\n")

    list(cit=x, ref.tex=ref.tex, val.tex=val.tex, cite.tex=cite.tex, collab.tex=expnt)
  })

  meas.by.ref.defs = lapply(meas.by.ref, function(x) {
    def.cite = paste("\\htdef{", x$cit, ".cite", "}{", x$cite.tex, "}%", sep="")
    def.collab = paste("\\htdef{", x$cit, ".collab", "}{", x$collab.tex, "}%", sep="")
    def.ref = paste("\\htdef{", x$cit, ".ref", "}{", x$ref.tex, "}%", sep="")
    def.meas = paste("\\htdef{", x$cit, ".meas", "}{%\n", x$val.tex, "}%", sep="")
    paste(def.cite, def.collab, def.ref, def.meas, sep="\n")
  })
  meas.by.ref.defs.tex = paste(meas.by.ref.defs, collapse="\n")

  meas.by.ref.tex = lapply(meas.by.ref, function(x) {
    rc = paste(
      "\\multicolumn{2}{l}{\\htuse{", x$cit, ".ref", "}} \\\\\n",
      "\\htuse{", x$cit, ".meas", "}",
      sep="")
  })
  meas.by.ref.tex = paste(meas.by.ref.tex, collapse=" \\\\\\hline\n")

  list(defs=meas.by.ref.defs.tex, table=meas.by.ref.tex)
}

##
## return latex code with the base node correlation coefficients
##
get.tex.base.nodes.corr = function() {
  tex.all.tau.br.corr = NULL
  ##--- names of all quantities
  quant.names = combination$combine
  ##--- get quantity names defined via constraint
  quant.constr.names = sub(".c", "", names(combination$constr.all.expr), fixed=TRUE)
  ##--- take out quantities defined with a constraint
  quant.names = setdiff(quant.names, quant.constr.names)
  ##--- also remove the unitarity complement
  quant.names = setdiff(quant.names, "Gamma998")
  ##--- sort by ascending Gamma number
  quant.names = quant.names[order(alurep.gamma.num.id(quant.names))]

  ##--- tex code preceding the correlation table content
  corr.pre = c(
    "%%",
    "%% base nodes correlation, @@num@@",
    "%%",
    "\\ifhevea\\begin{table}\\fi%% otherwise cannot have normalsize caption",
    "\\begin{center}",
    "\\ifhevea",
    "\\caption{Base nodes correlation coefficients in percent, section @@num@@\\label{tab:tau:br-fit-corr@@num@@}}%",
    "\\else",
    "\\begin{minipage}{\\linewidth}",
    "\\begin{center}",
    "\\captionof{table}{Base nodes correlation coefficients in percent, section @@num@@}\\label{tab:tau:br-fit-corr@@num@@}%",
    "\\fi",
    "\\begin{envsmall}",
    "\\begin{center}",
    "\\renewcommand*{\\arraystretch}{1.1}%",
    "\\begin{tabular}{@@tabcols@@}",
    "\\hline")

  ##--- tex code following the correlation table content
  corr.post = c(
    "\\\\\\hline",
    "\\end{tabular}",
    "\\end{center}",
    "\\end{envsmall}",
    "\\ifhevea\\else",
    "\\end{center}",
    "\\end{minipage}",
    "\\fi",
    "\\end{center}",
    "\\ifhevea\\end{table}\\fi")

  ##--- correlation of base nodes in percent
  quant.corr.base = quant.corr[quant.names, quant.names] * 100

  ##--- write data in tables with at most spec. rows and columns
  coeff.per.row = 14
  coeff.per.col = 14
  inum = 1
  for (j in seq(1, length(quant.names), by=coeff.per.col)) {
    for (i in seq(1, length(quant.names), by=coeff.per.row)) {
      submat.txt = NULL
      for(ii in i:min(i+coeff.per.row-1, length(quant.names))) {
        if (ii<=j) next
        maxcol = min(ii-1,j+coeff.per.col-1)
        row.txt = paste("\\(", alurep.gamma.texlabel(quant.names[ii]), "\\)")
        row.txt = c(row.txt, sprintf("%4.0f", quant.corr.base[ii, j:maxcol]))
        row.txt = c(row.txt, rep("", length.out=min(j+coeff.per.col-1, length(quant.names))-j+1-(maxcol-j+1)))
        submat.txt = c(submat.txt, paste(row.txt, collapse=" & "))
      }
      if (is.null(submat.txt)) next
      label.txt = alurep.gamma.texlabel(quant.names[j:min(j+coeff.per.col-1, length(quant.names))])
      label.num = min(j+coeff.per.col-1, length(quant.names)) - j + 1
      label.txt = paste("\\(", label.txt, "\\)", sep=" ", collapse=" & ")
      label.txt = paste("", label.txt, sep=" & ")
      submat.txt = paste(submat.txt, collapse=" \\\\\n")
      submat.txt = paste(submat.txt, label.txt, sep=" \\\\\n")
      corr.pre.mod = gsub("@@num@@", as.character(inum), corr.pre)
      corr.pre.mod = gsub("@@tabcols@@", paste(rep("r", length.out=label.num+1), collapse=""), corr.pre.mod)
      ## corr.pre.mod = paste(corr.pre.mod, collapse="\n")
      tex.all.tau.br.corr = c(tex.all.tau.br.corr, corr.pre.mod, submat.txt, corr.post)
      inum = inum + 1
    }
  }
  return(paste(tex.all.tau.br.corr, collapse="\n"))
}

##--- get all constraints names that were used in the fit
get.constraints.used.names = function() {
  constr.used = combination$constr.all.lin | combination$constr.all.nl
  constr.used.names = names(constr.used[constr.used])
  constr.order = order(alurep.gamma.num.id(constr.used.names))
  constr.used.names = constr.used.names[constr.order]
  return(constr.used.names)
}

##--- get all used constraint equations
get.tex.constraints.used = function() {
  constr.used.names = get.constraints.used.names()
  comb.str = combination$constr.all.str.expr[constr.used.names]
  comb.val = combination$constr.all.val[constr.used.names]
  comb.nl = intersect(constr.used.names, names(combination$constr.nl.str.expr))
  if (length(comb.nl) > 0) {
    comb.str[comb.nl] = combination$constr.nl.str.expr[comb.nl]
    comb.val[comb.nl] = combination$constr.nl.str.val[comb.nl]
  }
  return(alurep.tex.constraint(unlist(comb.val), unlist(comb.str)))
}

##
## return latex code with constraint equations
## (only constraints not corresponding to ratio of BRs)
##
get.tex.constraint.defs = function() {
  rc = get.tex.constraints.used()
  return(paste("\\htconstrdef{", names(rc$left), "}{", rc$left, "}{", rc$right, "}{", rc$right.split, "}%", sep=""))
}

##
## return latex code with constraint equations
## (only constraints not corresponding to ratio of BRs)
##
get.tex.constraint.equations = function() {
  rc = get.tex.constraints.used()
  constr.used.names = names(rc$left)

  constr.used.names = get.constraints.used.names()
  
  ##
  ## selection of constraints related to BRs that are ratios of 2 BRs
  ## in the report, these constraints are listed in the definition of the BRs not in the list of constraints
  ## also select the dummy constraint to compute the Unitarity discrepancy
  ##
  sel = grepl("Gamma\\d+by\\d+.*|Unitarity", constr.used.names, perl=TRUE)
  constr.names = constr.used.names[!sel]

  return(paste("\\begin{align*}\n",
               "\\htuse{", constr.names, ".left}",
               " ={}& ",
               "\\htuse{", constr.names, ".right.split}",
               "\n\\end{align*}", sep="", collapse="\n"))
}

##
## create .tex files for HFAG report
##
mkreport = function(fname) {
  load(fname, .GlobalEnv)

  ##--- build out file name
  fname = sub("[.][^.]*$", ".tex", fname, perl=TRUE)
  fname = sub("^[^.-]*", "tau-br-fit", fname, perl=TRUE)
  fname = file.path("../report", fname)
  cat("", file=fname)
  cat("file '", fname, "' created\n", sep="")

  ##
  ## HFAG-Tau fit uses some extra quantities not related to measurements:
  ## GammaAll, Gamma998 (unitarity residual), Gamma110
  ## we subtract these three to get the number of variables that
  ## are actually corresponding to a measurement or an estimate
  ##
  dummy.quant.num = 3

  ##
  ## write tex defs of some quantities
  ##
  tex.defs = c(
    ##--- unitarity check
    alurep.tex.cmd.short("UnitarityResid",
                         alurep.tex.val.err.auto(quant.val["Gamma998"], quant.err["Gamma998"], perc=TRUE)),
    ##--- measurements
    alurep.tex.cmd.short("MeasNum", as.character(meas.num)),
    ##--- quantities corresponding to measurements
    alurep.tex.cmd.short("QuantNum", as.character(quant.num-dummy.quant.num)),
    ##--- base quantities
    alurep.tex.cmd.short("BaseQuantNum", as.character(quant.num - constr.num)),
    ##--- constraints used to relate measurements to base quantities
    alurep.tex.cmd.short("ConstrNum", as.character(constr.num - dummy.quant.num)),
    alurep.tex.cmd.short("Chisq", sprintf("%.1f", chisq)),
    alurep.tex.cmd.short("Dof", as.character(dof)),
    alurep.tex.cmd.short("ChisqProb", alurep.tex.val.auto(chisq.prob, perc=TRUE))
    )
  cat(tex.defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', initial defs\n", sep="")

  ##
  ## define measurements
  ##
  rc = get.tex.meas.defs()
  cat(rc, file=fname, append=TRUE)
  cat("\n", file=fname, append=TRUE)
  cat("file '", fname, "', measurement description defs\n", sep="")

  ##
  ## write tex macro containing BR fit data
  ##
  quant.names = combination$combine
  quant.names = setdiff(quant.names, "GammaAll")
  rc = get.tex.table(quant.names)

  cat(rc$defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', quantity - measurements defs\n", sep="")

  tex.all.tau.br.val = alurep.tex.cmd("BrVal", rc$table)
  cat(tex.all.tau.br.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', quantity - measurements table\n", sep="")

  ##
  ## write tex macro containing all measurements by reference
  ##
  rc = get.tex.meas.by.ref()
  cat(rc$defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR meas by ref definitions\n", sep="")
  ##
  tex.meas.paper = alurep.tex.cmd("MeasPaper", rc$table)
  cat(tex.meas.paper, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR meas by ref table content\n", sep="")
  rm(rc)

  ##
  ## write text macro containing all strange BR values and refs
  ##
  gamma110.names = names(combination$constr.all.comb$Gamma110.c)
  gamma110.names = setdiff(gamma110.names, "Gamma110")
  tex.tau.br.strange.val = alurep.tex.cmd("BrStrangeVal", get.tex.table.simple(gamma110.names, 4, -2))
  cat(tex.tau.br.strange.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR strange table content\n", sep="")

  ##
  ## write text macro containing total strange BR
  ##
  gamma110.names = "Gamma110"
  tex.tau.br.strange.tot.val = alurep.tex.cmd("BrStrangeTotVal", get.tex.table.simple(gamma110.names, 4, -2))
  cat(tex.tau.br.strange.tot.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR strange tot table content\n", sep="")

  ##
  ## write text macro containing all strange BR values and refs
  ##
  gammaAll.names = names(combination$constr.all.comb$GammaAll.c)
  gammaAll.names = c(gammaAll.names, "Gamma998")
  gammaAll.names = setdiff(gammaAll.names, "GammaAll")
  tex.tau.unitarity.quants = alurep.tex.cmd("UnitarityQuants", get.tex.table.simple(gammaAll.names, 4, -2))
  cat(tex.tau.unitarity.quants, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', unitarity quantities\n", sep="")

  ##
  ## write text macro containing correlation of base nodes
  ##
  tex.all.tau.br.corr = alurep.tex.cmd("BrCorr", get.tex.base.nodes.corr())
  cat(tex.all.tau.br.corr, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', BR correlations table content\n", sep="")

  ##
  ## write TeX definitions for all constraint equations
  ##
  tex.constr.defs = get.tex.constraint.defs()
  cat(tex.constr.defs, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', constraint definitions\n", sep="")
  
  ##
  ## write text macro containing constraint equations
  ##
  tex.constr.val = alurep.tex.cmd("ConstrEqs", get.tex.constraint.equations())
  cat(tex.constr.val, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', constraint table content\n", sep="")

  ##
  ## measurements by collaboration
  ##
  tex.num.meas.per.collab = get.tex.meas.by.collab()
  cat(tex.num.meas.per.collab, sep="\n", file=fname, append=TRUE)
  cat("file '", fname, "', measurements per collaboration\n", sep="")
}

## ////////////////////////////////////////////////////////////////////////////
## code

args = commandArgs(TRUE)
if (length(args) == 1) {
  rc = mkreport(fname = args[1])
} else {
  mkreport(fname = "average.rdata")
}
