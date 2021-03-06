#!/usr/bin/env Rscript

##
## aluelab-results.r [flags] [<.rdata file>]
## elaborate hfag tau 2011 results for Vus, lepton univ
##
## winter 2012 report:
##   tau/2011/TauFit> ../scripts/aluelab2-results.r -vadirect average2-aleph-hcorr.rdata
##   will produce '../report/tau-br-fit-aleph-hcorr-vadirect-elab.tex'
##
## reproduce 2009/TauFit-Aug2011:
##   replace old f_K and f_K/f_pi
##   tau/2011/TauFit> ../scripts/aluelab2-results.r -lepuniv -u average2-aleph-hcorr.rdata
##

require(stringr, quietly=TRUE)
require(Matrix, quietly=TRUE)

source("../../../Common/bin/aluelab3.r")
source("../../../Common/bin/alureport.r")

## ////////////////////////////////////////////////////////////////////////////
## functions

##
## get list of quantity names in the definition of a quantity
## using the relevant constraint equation
## used since Oct 2010
## since March 2012, use aluelab.get.str.expr, needs alucomb2 format
##
aluelab.get.quant.names = function(quant.name, combination) {
  var.comb = combination$constr.all.comb[[paste(quant.name, ".c", sep="")]]
  if (is.null(var.comb) || is.na(var.comb)) {
    var.comb = combination$constr.all.comb[[paste(quant.name, ".coq", sep="")]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.lin.comb[[paste(quant.name, ".coq", sep="")]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.all.comb[[quant.name]]
  }
  if (is.null(var.comb)) {
    var.comb = combination$constr.lin.comb[[quant.name]]
  }
  var.comb = var.comb[names(var.comb) != quant.name]
  var.names = names(var.comb[var.comb != 0])
  return(var.names)
}

##
## get string expression for a quantity, from its constraint equations
## needs alucomb2 format
##
aluelab.get.str.expr = function(quant.name, combination) {
  str.expr = combination$constr.all.str.expr[[paste(quant.name, ".c", sep="")]]
  str.expr = str_match(str.expr, "-([[:alnum:]]+)\\s+[+]\\s+(.*\\S)\\s*$")[,3]
  ##--- remove outer braces
  str.expr = gsub("(\\(((?:[^()]++|(?1))+)*+\\))", "\\2", str.expr, perl=TRUE)
  return(str.expr)
}

## ////////////////////////////////////////////////////////////////////////////
## initialization

display.names = character(0)

## ////////////////////////////////////////////////////////////////////////////
## code

aluelab.results = function(args) {
  ##--- option -s, use S-factors
  flag.s.factors =  FALSE
  if(any(args == "-h")) {
    ##--- help
    cat("aluelab-results.r [flags] [<.rdata file>]\n")
    cat("  no option: traditional Vus but without universality improved Be for B_had\n")
    cat("  -s use error scale factors\n")
    cat("  -u use unitarity constraint to compute Be_univ, B_had_VA, B_had_s\n")
    cat("  -vadirect determine B_VA directly rather than from 1-Be-Bmu-Bstrange\n")
    cat("  -kmaltman use delta_R from Kim Maltman, private communication, 2014\n")
    cat("  -lepuniv use R_had = (1 - (1+(f_mu/f_e))Be_univ)/Be_univ\n")
    ## args = args[args != "-h"]
    stop()
  }
  if(any(args == "-s")) {
    ##--- use S-factors
    flag.s.factors =  TRUE
    args = args[args != "-s"]
  }
  flag.unitarity =  FALSE
  if(any(args == "-u")) {
    ##--- use unitarity constraint to compute Be_univ, B_had_VA, B_had_s
    flag.unitarity =  TRUE
    args = args[args != "-u"]
  }
  flag.vadirect =  FALSE
  if(any(args == "-vadirect")) {
    ##--- determine B_VA from direct measurements rather than 1-Be-Bmu-Bstrange
    flag.vadirect =  TRUE
    args = args[args != "-vadirect"]
  }
  flag.lepuniv =  FALSE
  if(any(args == "-lepuniv")) {
    ##--- use R_had = (1 - (1+(f_mu/f_e))Be_univ)/Be_univ
    flag.lepuniv =  TRUE
    args = args[args != "-lepuniv"]
  }
  flag.kmaltman =  FALSE
  if(any(args == "-kmaltman")) {
    ##--- use delta_R from Kim Maltman, private communication, 2014
    flag.kmaltman =  TRUE
    args = args[args != "-kmaltman"]
  }

  if (length(args) > 0) {
    file.name = args[1]
  } else {
    file.name = "average.rdata"
  }

  ##--- get alucomb results and data
  load(file.name)

  if (flag.s.factors) {
    quant.cov = quant.sf.cov
  }

  quant.err = sqrt(diag(quant.cov))
  quant.err.nonzero = (quant.err != 0)
  quant.corr = quant.cov
  quant.corr[quant.err.nonzero, quant.err.nonzero] = cov2cor(quant.cov[quant.err.nonzero, quant.err.nonzero])

  quant = StatComb$new(quant.val, quant.cov)

  quant.names = names(quant.val)
  comb.params = lapply(combination$params, function(x) unname(x["value"]))
  quant$param.add(comb.params)
  quant$param.add(c(pi=pi))

  quant$texdescr.add(alurep.get.texdescr.quantities(combination$quantities))

  ##
  ## recover some BRs as function of others
  ##
  if (!any("Gamma89" == quant.names)) {
    rc = quant$quant.expr.add("Gamma89", Gamma803 + BR_om_pimpippiz*Gamma151)
  }

  ##
  ## PDG 2009 definition of Gamma110 = B(tau -> Xs nu)
  ##
  Gamma110_pdg09.comb = c(
    Gamma10=1,  Gamma16=1,
    Gamma23=1,  Gamma28=1,  Gamma35=1,
    Gamma40=1,  Gamma85=1,  Gamma89=1,
    Gamma128=1)
  quant$quant.comb.add("Gamma110_pdg09", Gamma110_pdg09.comb)

  ##
  ## Gamma110 = B(tau -> Xs) -- local definition
  ##
  ## Gamma110 =
  ##   Gamma10+Gamma16+Gamma23+Gamma28+Gamma35+Gamma40+Gamma44+
  ##   Gamma53+Gamma85+Gamma89+Gamma128+Gamma130+Gamma132+
  ##   Gamma96+Gamma96*B(phi -> KS0 KL0)/B(phi -> K+ K-)
  ##
  Gamma110.comb = c(
    Gamma10=1,  Gamma16=1,
    Gamma23=1,  Gamma28=1,  Gamma35=1,
    Gamma40=1,  Gamma44=1,  Gamma53=1,
    Gamma85=1,  Gamma89=1,  Gamma128=1,
    Gamma130=1, Gamma132=1,
    Gamma96=1.699387)

  ##
  ## Oct 2010, Gamma110 updated as follows
  ##
  ## Gamma110 = Gamma10  + Gamma16   + Gamma23   + Gamma28  + Gamma35  + Gamma40 + Gamma85 + Gamma89 + Gamma128
  ##          + Gamma151 + Gamma130  + Gamma132  + Gamma44  + Gamma53  + Gamma801
  ##
  ## - replaced Gamma96 with Gamma801 = Gamma96 * [1 + B(phi -> KS0 KL0)/B(phi -> K+ K-)]
  ##   because of technical issues with readpdg.cc
  ## - added Gamma151 (tau -> K omega nu)
  ##
  Gamma110.comb = c(
    Gamma10=1,  Gamma16=1,  Gamma23=1,
    Gamma28=1,  Gamma35=1,  Gamma40=1,
    Gamma44=1,  Gamma53=1,  Gamma85=1,
    Gamma89=1,  Gamma128=1, Gamma130=1,
    Gamma132=1, Gamma151=1, Gamma801=1)

  ##
  ## end Oct2010 update
  ## define Gamma110 using the alucomb.r Gamma110 constraint
  ## this is more robust against changes in the fitting
  ##
  Gamma110.names = aluelab.get.quant.names("Gamma110", combination)

  ##
  ## March 2012, get Gamma110 from its constraint equation
  ##
  Gamma110.str.expr = aluelab.get.str.expr("Gamma110", combination)
  Gamma110.comb = quant$str.to.comb(Gamma110.str.expr)
  Gamma110.names = names(Gamma110.comb)

  Gamma110.str.expr = paste(Gamma110.comb, Gamma110.names, sep="*", collapse=" + ")
  cat("Gamma110 = ", Gamma110.str.expr, "\n", sep="")

  if (!any("Gamma151" == quant.names)) {
    ##--- remove Gamma151 (tau -> K omega nu) if not present in fit for backward compatibility
    Gamma110.comb = Gamma110.comb[Gamma110.names != "Gamma151"]
    cat("warning, Gamma151 removed from Gamma110 definition\n")
  }

  ##--- list of all tau BRs that are not leptonic and not strange, i.e. not-strange-hadronic
  ## B.tau.VA.names = setdiff(aluelab.get.quant.names("GammaAll", combination), Gamma110.names)
  ## B.tau.VA.names = setdiff(B.tau.VA.names, c("Gamma5", "Gamma3", "Gamma998"))
  ## quant$quant.qexpr.add("B_tau_VA", parse(text=paste(B.tau.VA.names, collapse="+")))

  ##  cannot assume any more that GammaAll is sum of quantities all with coefficient 1
  ##  ##--- compute expression of all tau BRs that are not leptonic and not strange, i.e. not-strange-hadronic
  ##  GammaAll.str.expr = gsub("-GammaAll\\s*\\s*[+]\\s*\\((.*)\\)", "\\1", combination$constr.all.str.expr[["GammaAll.c"]])
  ##  GammaAll.term.names = str_trim(str_extract_all(GammaAll.str.expr, "([^+]+)")[[1]])
  ##  Gamma110.names.not.in.GammaAll = Gamma110.names[!(Gamma110.names %in% GammaAll.term.names)]
  ##  ##--- GammaAll includes Gamma103 instead of its terms, which include Gamma822 and Gamma833 in Gamma110
  ##  B.tau.VA.names = c(
  ##    setdiff(GammaAll.term.names, c("Gamma3", "Gamma5", Gamma110.names)),
  ##    paste("-", Gamma110.names.not.in.GammaAll, sep=""))
  ##  B.tau.VA.str.expr = gsub("[+](\\s*)-", "-\\1", paste(B.tau.VA.names, collapse=" + "))
  ##  quant$quant.qexpr.add("B_tau_VA", parse(text=B.tau.VA.str.expr))
  
  quant$quant.expr.add("B_tau_VA", GammaAll-Gamma5-Gamma3-Gamma110)

  ##
  ## using unitarity, any BR has an additional estimate as BR_uni = 1 - (BR_all - BR)) 
  ##
  quant$quant.expr.add("B_tau_VA_unitarity", 1 - (GammaAll - B_tau_VA))
  quant$quant.expr.add("B_tau_s_unitarity",  1 - (GammaAll - Gamma110))

  ##
  ## add measurements to compute universality improved Be
  ##

  ##--- from PDG 2015 ---upd16
  quant$quant.add.single("m_e",0.510998928, 0.000000011) ## changed14
  quant$quant.add.single("m_mu", 105.6583715, 0.0000035) ## changed14

  ##--- from PDG 2015 ---upd16
  quant$quant.add.single("tau_tau", 290.29e-15, 0.52e-15) ## changed14
  ## quant$quant.add.single("tau_tau", 290.6e-15, 1.0e-15)

  ##--- m_tau PDG 2015 ---upd16
  quant$quant.add.single("m_tau", 1776.86, 0.12)
  ##--- m_tau PDG 2013
  ## quant$quant.add.single("m_tau", 1776.82, 0.16)
  ##--- m_tau HFAG 2009
  ## quant$quant.add.single("m_tau", 1776.7673082, 0.1507259)

  ##--- from PDG 2015 ---upd16
  quant$quant.add.single("m_pi", 139.57018, 0.00035)
  quant$quant.add.single("tau_pi", 2.6033e-8, 0.0005e-8)
  quant$quant.add.single("m_K", 493.677, 0.016)
  quant$quant.add.single("tau_K", 1.2380e-8, 0.0020e-8)

  ##--- from PDG 2015, 2013, 2011 ---upd16
  quant$quant.add.single("m_W", 80.385e3, 0.015*1e3) ## changed14
  quant$quant.add.single("tau_mu", 2.1969811e-6, 0.000022e-6) ## changed14

  ##
  ## Be, from unitarity = 1 - Bmu - B_VA - B_s
  ## Bmu, from unitarity = 1 - Be - B_VA - B_s
  ##
  quant$quant.expr.add("Be_unitarity", 1 - Gamma3 - B_tau_VA - Gamma110)
  quant$quant.expr.add("Bmu_unitarity", 1 - Gamma5 - B_tau_VA - Gamma110)

  ##--- understand if there was no unitarity constraint
  no.unit.constr.flag =
    ("Gamma998" %in% aluelab.get.quant.names("GammaAll", combination) ||
     "Gamma998" %in% aluelab.get.quant.names("Unitarity", combination))

  if (no.unit.constr.flag) {
    ##
    ## unconstrained fit
    ## fit best values for Be, Bmu, B_tau_VA, B_tau_s
    ## using both the direct measurements and the result from the unitarity constraint
    ##
    if (flag.unitarity) {
      ##--- use unitarity constraint to compute Be, Bmu, B_tau_VA/s
      quant$quant.fit.add("Be_fit", c(Gamma5=1, Be_unitarity=1))
      quant$quant.fit.add("Bmu_fit", c(Gamma3=1, Bmu_unitarity=1))
      quant$quant.fit.add("B_tau_VA_fit", c(B_tau_VA=1, B_tau_VA_unitarity=1))
      quant$quant.fit.add("B_tau_s_fit", c(Gamma110=1, B_tau_s_unitarity=1))
    } else {
      ##--- direct measurements for leptonic BRs and strange hadronic BR
      quant$quant.expr.add("Be_fit", Gamma5)
      quant$quant.expr.add("Bmu_fit", Gamma3)
      quant$quant.expr.add("B_tau_s_fit", Gamma110)
      if (flag.vadirect) {
        ##--- direct measurement for non-strange hadronic BR
        quant$quant.expr.add("B_tau_VA_fit", B_tau_VA)
      } else {
        ##
        ## non-strange hadronic BR by subtraction, however note that
        ## here we do not use the improved Be from universality
        ## (a Be fit using Be, Bmu and tau lifetime) as in
        ## M.Davier et al, RevModPhys.78.1043, arXiv:hep-ph/0507078
        ##
        quant$quant.expr.add("B_tau_VA_fit", 1-Gamma5-Gamma3-Gamma110)
      }
    }
  } else {
    ##
    ## unitarity constrained fit
    ## the fitted values have the unitarity constraint already
    ##
    rc = quant$quant.expr.add("Be_fit", Gamma5)
    rc = quant$quant.expr.add("Bmu_fit", Gamma3)
    rc = quant$quant.expr.add("B_tau_VA_fit", B_tau_VA)
    rc = quant$quant.expr.add("B_tau_s_fit", Gamma110)
  }

  rc = quant$quant.expr.add("B_tau_had_fit", B_tau_VA_fit+B_tau_s_fit)

  ##
  ## compute phase space factors for Bmu/Be universality
  ##
  ##--- phase space factor, function of lepton masses
  phspf = quote(1 -8*x + 8*x^3 - x^4 - 12*x^2*log(x))
  ##--- phase space factors for e/tau, mu/tau, e/mu
  rc = quant$quant.qexpr.add("phspf_mebymtau",  esub.expr(phspf, list(x=quote(m_e^2/m_tau^2))))
  rc = quant$quant.qexpr.add("phspf_mmubymtau", esub.expr(phspf, list(x=quote(m_mu^2/m_tau^2))))
  rc = quant$quant.qexpr.add("phspf_mebymmu", esub.expr(phspf, list(x=quote(m_e^2/m_mu^2))))
  rc = quant$quant.expr.add("Bmu_by_Be_th", phspf_mmubymtau/phspf_mebymtau)

  ##--- Be from Bmu
  quant$quant.expr.add("Be_from_Bmu", phspf_mebymtau/phspf_mmubymtau * Bmu_fit)

  ##
  ## rad. corrections from to get Be from tau lifetime
  ## values from 10.1103/RevModPhys.78.1043 p.1047, arXiv:hep-ph/0507078v2 p.7, could be recomputed
  ## - delta^L_gamma = 1 + alpha(mL)/2pi * (25/4 - pi^2)
  ## - delta^L_W = 1 + 3/5* m_L^2/M_W^2 + 9/5 m_l^2/M_W^2
  ## ---upd14
  ## - 2014: no update available
  quant$param.add(c(delta_mu_gamma=(1 - 42.4e-4), delta_tau_gamma=(1 - 43.2e-4)))
  quant$quant.expr.add("delta_mu_W", 1 + 3/5*m_mu^2/m_W^2 + 9/5*m_e^2/m_W^2)
  quant$quant.expr.add("delta_tau_W", 1 + 3/5*m_tau^2/m_W^2 + 9/5*m_mu^2/m_W^2)

  ##
  ## Be from tau lifetime
  ## Be= tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_e/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  ##
  quant$quant.expr.add("Be_from_taulife",
                      tau_tau/tau_mu * (m_tau/m_mu)^5 * phspf_mebymtau/phspf_mebymmu
                      * (delta_tau_gamma*delta_tau_W) / (delta_mu_gamma*delta_mu_W))
  ##
  ## Bmu from tau lifetime
  ## Bmu= tau_tau/tau_mu (m_tau/m_mu)^5 f(m^2_mu/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  ##
  quant$quant.expr.add("Bmu_from_taulife",
                      tau_tau/tau_mu * (m_tau/m_mu)^5 * phspf_mmubymtau/phspf_mebymmu
                      * (delta_tau_gamma*delta_tau_W) / (delta_mu_gamma*delta_mu_W))

  ##
  ## Be_lept = B(tau -> e nu nubar (gamma)) improved using also Bmu
  ## - useful to have Be for universality plot vs. tau_tau
  ## - minimum chisq fit using, Be, Be from Bmu
  ##
  quant$quant.fit.add("Be_lept", c(Gamma5=1, Be_from_Bmu=1))

  ##
  ## Be_univ, universality improved Be = B(tau -> e nu nubar (gamma))
  ## - see arXiv:hep-ph/0507078v2 p.7, doi:10.1103/RevModPhys.78.1043 p.1047
  ## - minimum chisq fit for Be_univ using, Be, Be from Bmu, Be from tau lifetime
  ##
  quant$quant.fit.add("Be_univ", c(Gamma5=1, Be_from_Bmu=1, Be_from_taulife=1))

  ##
  ## Vud ---upd16
  ## Phys. Rev. C91 (2) (2015) 025501
  ## Superallowed 0+→0+0+→0+ nuclear β decays: 2014 critical survey, with precise results for VudVud and CKM unitarity
  ## J.C. Hardy, I.S. Towner (Texas A-M, Cyclotron Inst.)
  ##
  ## Vud ---upd14
  ## PDG13: no change vs PDG11
  ##
  ## arXiv:0710.3181v1 [nucl-th], 10.1103/PhysRevC.77.025501
  ## I.S.Towner, J.C.Hardy, An improved calculation of the isospin-symmetry-breaking corrections to superallowed Fermi beta decay
  ## also PDG 2013 review
  ##
  ## Vud= 0.97425 (8)exp. (10)nucl.dep. (18)RC
  ## only superallowed Beta decays, Hardy & Towner
  ##
  ## Vud = 0.9774 (5)tau-n (16)gA (2)RC
  ## including neutron lifetime

  ##
  ## PDG11: quant$quant.add.single("Vud", 0.97425, 0.00022)
  ## PDG13 with neutron lifetime: quant$quant.add.single("Vud", 0.9774, 0.0001*sqrt(5^2 + 16^2 + 2^2))
  ## PDG15 unchanged
  ## quant$quant.add.single("Vud", 0.97425, 0.00001*sqrt(8^2 + 10^2 + 18^2))
  ## TH2015 0.97417(5)exp.(9)nucl.dep.(18)RC, also 2016 PDG review
  quant$quant.add.single("Vud", 0.97417, 0.00001*sqrt(5^2 + 9^2 + 18^2))
  
  ##--- Moulson CKM2014 procs
  quant$quant.add.single("Vud_moulson_ckm14", 0.97417, 0.00001*21)
  ##--- Moulson CKM2014 procs, using FLAG Nf = 2+1+1
  quant$quant.add.single("Vus_kl3_moulson_ckm14", 0.2232, 0.0001*9)
  ##--- Moulson CKM2014 procs, using his elaboration using FNAL/MILC result Nf = 2+1+1, fK± /fπ± = 1.1960(25)
  quant$quant.add.single("VusbyVud_moulson_ckm14", 0.2308, 0.0001*6)
  quant$quant.expr.add("Vus_kl2_moulson_ckm14", VusbyVud_moulson_ckm14 * Vud_moulson_ckm14)

  ##
  ## Kim Maltman Mainz March 2016 slides
  ## Nf = 2+1+1
  ##
  quant$quant.add.single("Vus_maltman_mainz16", 0.2228, sqrt(0.0023^2 + 0.0005^2))
  quant$quant.add.single("Vus_kl3_maltman_mainz16", 0.2231, sqrt(0.0004^2 + 0.0007^2))
  quant$quant.add.single("Vus_kl2_maltman_mainz16", 0.2250, 0.000985)
  
  ##--- s quark mass, PDG2011
  ## quant$quant.add.single("m_s", 100, sqrt((20.^2 + 30.^2)/2.))
  ##--- s quark mass, PDG2013 ---upd14
  ## quant$quant.add.single("m_s", 93.5, 2.5) ## changed14
  ##--- s quark mass, PDG2015 ---upd16
  quant$quant.add.single("m_s", 95, 5) ## changed16

  if (flag.kmaltman) {
    ## 
    ## K.Maltman, private comm. 2014
    ## - use tau pi Kl2 poles
    ## - FOPT 4 order
    ## the following is just for the record, the value computed by K.Maltman is used for deltaR_su3break
    ## 
    quant$quant.add.single("deltaR_su3break_pheno", 0.1509, 0.0041) 
    quant$quant.add.single("deltaR_su3break_d2pert", 11.36, 1.14)
    quant$quant.add.single("deltaR_su3break_remain", 0.0084, 0.0185)
    
    ## 
    ## K.Maltman, private comm. 2014
    ## - use FOPT 3rd loop order as central value (favoured by Lattice)
    ## - use FOPT uncertainty
    ## - add FOPT-CIPT 3rd order discrepancy as additional theory uncertainty
    ##
    ## quant$quant.expr.add("deltaR_su3break", deltaR_su3break_pheno + deltaR_su3break_d2pert*(m_s/1000)^2 + deltaR_su3break_remain)
    quant$quant.add.single("deltaR_su3break", 0.254, 0.038)
  } else {
    ##
    ## SU3 breaking correction
    ##
    
    ##
    ## POS(KAON)08, A.Pich, Theoretical progress on the Vus determination from tau decays
    ##
    ## deltaR.su3break.val = 0.216
    ## deltaR.su3break.err = 0.016
    
    ##
    ## use deltaR_su3break from
    ## E. Gamiz et al., Nucl.Phys.Proc.Suppl.169:85-89,2007, arXiv:hep-ph/0612154v1
    ##
    ## deltaR.su3break.val = 0.240
    ## deltaR.su3break.err = 0.032
    ## quant$quant.add.single("deltaR_su3break", deltaR.su3break.val, deltaR.su3break.err)
    
    ##
    ## recompute deltaR_su3break using updated m_s
    ## E. Gamiz et al., Nucl.Phys.Proc.Suppl.169:85-89,2007, arXiv:hep-ph/0612154v1
    ##
    ## deltaR_su3break = deltaR_su3break_pheno + deltaR_su3break_d2pert * m_s(GeV)^2 + deltaR_su3break_remain
    ## - m_s = 93.5 +- 2.5 (PDG13) strange quark mass in GeV in the MSbar scheme at a renormalisation scale of mu = 2 GeV
    ## - deltaR_su3break_pheno = 0.1544 +- 0.0037 phenomenological scalar and pseudo-scalarcontributions
    ## - deltaR_su3break_d2pert = 9.3 +- 3.4 rest of the perturbative D=2 contribution
    ## - deltaR_su3break_remain = 0.0034 +- 0.0028 remaining contributions
    ##
    
    ##--- E.Gamiz, M.Jamin, A.Pich, J.Prades, F.Schwab, |V_us| and m_s from hadronic tau decays
    quant$quant.add.single("deltaR_su3break_pheno", 0.1544, 0.0037)
    quant$quant.add.single("deltaR_su3break_d2pert", 9.3, 3.4)
    quant$quant.add.single("deltaR_su3break_remain", 0.0034, 0.0028)

    quant$quant.expr.add("deltaR_su3break", deltaR_su3break_pheno + deltaR_su3break_d2pert*(m_s/1000)^2 + deltaR_su3break_remain)
  }
  
  if (no.unit.constr.flag) {
    ##
    ## if using constrained fit with dummy mode (Gamma998), i.e. unconstrained fit
    ##
    if (flag.lepuniv) {
      ##--- determine R_tau using leptonic BRs and universality
      quant$quant.expr.add("R_tau", 1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau)
      quant$quant.expr.add("R_tau_s", B_tau_s_fit/Be_univ)
      quant$quant.expr.add("R_tau_VA", R_tau - R_tau_s)
    } else {
      ##
      ## use "fit" values computed in this script (possibly incorporating unitarity constraint)
      ##
      ##--- add R_tau_VA = R - R_tau_s
      quant$quant.expr.add("R_tau_VA", B_tau_VA_fit/Be_univ)
      ##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
      quant$quant.expr.add("R_tau_s", B_tau_s_fit/Be_univ)
      ##--- add R_tau as function of quantities
      quant$quant.expr.add("R_tau", R_tau_VA+R_tau_s)
    }
  } else {
    ##
    ## if using constrained fit without dummy mode (Gamma998), i.e. constrained fit
    ##
    ##--- add R_tau as function of quantities
    if (flag.lepuniv) {
      ##--- determine R_tau using leptonic BRs and universality
      quant$quant.expr.add("R_tau", 1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau)
    } else {
      ##--- use values computed in the alucomb.r fit, incorporating unitarity constraints
      quant$quant.expr.add("R_tau", (B_tau_VA+Gamma110)/Be_univ)
    }
    ##--- add R_tau_s = B(tau -> Xs nu) / Be_univ
    quant$quant.expr.add("R_tau_s", Gamma110/Be_univ)
    ##--- add R_tau_VA = R_tau - R_tau_s
    quant$quant.expr.add("R_tau_VA", R_tau - R_tau_s)
  }

  ##--- R_tau using only leptonic BRs
  quant$quant.expr.add("R_tau_leptuniv", 1/Be_univ -1 -phspf_mmubymtau/phspf_mebymtau)
  quant$quant.expr.add("R_tau_leptonly", (1-Gamma3-Gamma5)/Be_univ)

  ## ////////////////////////////////////////
  ##
  ## Vus from tau -> s inclusive
  ##

  ##--- add Vus
  quant$quant.expr.add("Vus", sqrt(R_tau_s/(R_tau_VA/Vud^2 - deltaR_su3break)))
  quant$param.add.single("Vus_err_perc", quant$err.contrib.perc("Vus"))

  ##--- theory, exp errors, with percent
  quant$param.add.single("Vus_err_th", quant$err.contrib("Vus", "deltaR_su3break"))
  quant$param.add.single("Vus_err_th_perc", quant$err.contrib.perc("Vus", "deltaR_su3break"))

  Vus.err.exp = sqrt(quant$err("Vus")^2 - quant$param("Vus_err_th")^2)
  quant$param.add.single("Vus_err_exp", Vus.err.exp)
  quant$param.add.single("Vus_err_exp_perc", Vus.err.exp/quant$val("Vus")*100)

  ##--- Vus from Vud using CKM unitarity
  quant$quant.expr.add("Vus_uni", sqrt(1-Vud^2))

  ##--- Vus vs <Vus from CKM unitarity>
  quant$quant.expr.add("Vus_mism", Vus - Vus_uni)
  Vus_mism_sigma = quant$val("Vus_mism") / quant$err("Vus_mism")
  quant$param.add.single("Vus_mism_sigma", Vus_mism_sigma)
  quant$param.add.single("Vus_mism_sigma_abs", abs(Vus_mism_sigma))

  ##--- add quantities to print
  display.names = c(
    display.names,
    ## Gamma110.names,
    "Gamma5", "Be_unitarity", "Be_fit",
    "Gamma3", "Bmu_unitarity", "Bmu_fit",
    "Bmu_by_Be_th", "Be_from_Bmu", "Be_from_taulife",
    "Be_lept", "Be_univ",
    "Bmu_from_taulife", "B_tau_had_fit",
    "B_tau_VA", "B_tau_VA_unitarity", "B_tau_VA_fit",
    "Gamma110", "B_tau_s_unitarity", "B_tau_s_fit",
    "Gamma110_pdg09",
    "R_tau", "R_tau_s", "R_tau_VA",
    "deltaR_su3break",
    "Vus",
    "Vus_err_exp", "Vus_err_th",
    "Vus_err_perc", "Vus_err_exp_perc", "Vus_err_th_perc",
    "Vus_uni", "Vus_mism_sigma"
    )
  ## display.names = c(Gamma110.names, display.names)

  ## ////////////////////////////////////////
  ##
  ## radiative corrections parameters
  ##

  ##
  ## ratio of radiative corrections for G(K -> mu nu) / G(pi -> mu nu)
  ##
  ## Marciano:2004uf
  ## W. J. Marciano, "Precise determination of |V(us)| from lattice calculations of pseudoscalar decay constants",
  ## Phys. Rev. Lett. 93:231803, 2004, doi:10.1103/PhysRevLett.93.231803, arXiv:hep-ph/0402299.
  ##
  ## actually relying on calculations for rad corr of G(h -> mu nu) by
  ## M. Finkemeier, Phys. Lett. B387, 391 (1996)
  ##
  ## ---upd16
  ## this is now superseeded by dRrad_kmunu_by_pimunu calculated by A. Pich
  quant$quant.add.single("Rrad_kmunu_by_pimunu", 0.9930, 0.0035)

  ##
  ## from A. Pich, “Precision Tau Physics,” Prog. Part. Nucl. Phys. 75 (2014) 41–85, arXiv:1310.7922 [hep-ph].
  ## includes previous calculation by Marciano and Cirigliano delta_SU2 correction
  ##
  quant$quant.add.single("dRrad_kmunu_by_pimunu", -0.0113, 0.0023)

  ##
  ## ratios of radiative corrections R_tau/pi = Gamma(tau -> h nu) / Gamma(h -> mu nu)
  ## for h = pi, K
  ## R_tau/h = 1 + deltaR_tau/h
  ##
  ## Decker:1994dd
  ## R. Decker, M. Finkemeier
  ## Radiative corrections to the decay tau --> pi (K) tau-neutrino.
  ## Phys.Lett. B334 (1994) 199-202
  ## http://inspirehep.net/record/383287?ln=en
  ##
  ## ---upd16
  ## use average of asymmetric error from Decker 1994 paper
  ## quant$quant.add.single("dRrad_taupi_by_pimu", 0.16e-2, sqrt( (0.09e-2^2 + 0.14e-2^2)/2 ))
  ## quant$quant.add.single("dRrad_tauK_by_Kmu", 0.90e-2, sqrt( (0.17e-2^2 + 0.26e-2^2)/2 ))
  ##
  ## --- Summer16
  ## use conservative values used also in A.Pich, Precision Tau Physics
  quant$quant.add.single("dRrad_taupi_by_pimu", 0.16e-2, 0.14e-2)
  quant$quant.add.single("dRrad_tauK_by_Kmu", 0.90e-2, 0.22e-2)

  ##--- compute Rrad_tauK_by_taupi = (1+dRrad_tauK_by_Kmu)/(1+dRrad_taupi_by_pimu) * (1+dRrad_kmunu_by_pimunu
  ## quant$quant.expr.add("Rrad_tauK_by_taupi", (1+dRrad_tauK_by_Kmu)/(1+dRrad_taupi_by_pimu) * Rrad_kmunu_by_pimunu)
  quant$quant.expr.add("Rrad_tauK_by_taupi", (1+dRrad_tauK_by_Kmu)/(1+dRrad_taupi_by_pimu) * (1+dRrad_kmunu_by_pimunu))

  ## ////////////////////////////////////////
  ##
  ## QCD lattice
  ##

  ##
  ## ---upd16
  ## Lattice averages FLAG 2016 Aoki et al.
  ## http://arxiv.org/abs/1607.00299
  ## https://inspirehep.net/record/1473344, Aoki:2016frl
  ## Nf = 2 + 1 + 1 : fK± /fπ± = 1.193(3) Refs. [66, 67, 69], (64) QCD with broken isospin
  ## Nf = 2 + 1 + 1 : fK± = 155.6 (0.4) MeV Refs. [66, 67, 69],
  ## Nf = 2 + 1 :     f+(0) = 0.9677(27) Refs. [37, 39], (59)
  ##
  quant$quant.add.single("f_K_by_f_pi", 1.193, 0.003) ## changed16
  quant$quant.add.single("f_K", 155.6, 0.4) ## changed16
  quant$quant.add.single("fp0_Kpi", 0.9677, 0.0027) ## changed16

  ##
  ## ---upd14
  ## Lattice averages FLAG 2013 Aoki et al.
  ## http://arxiv.org/abs/1310.8555
  ## http://inspirehep.net/record/1262813, Aoki:2013ldr
  ## use values for Nf = 2+1 (compatible with Nf = 2+1+1
  ##
  ## quant$quant.add.single("f_K_by_f_pi", 1.192, 0.005) ## unchanged12
  ## quant$quant.add.single("f_K", 156.3, 0.9) ## changed14
  ## quant$quant.add.single("fp0_Kpi", 0.9661, 0.0032) ## new since 2014

  ##
  ## Lattice averages from http://arxiv.org/abs/0910.2928 and
  ## http://krone.physik.unizh.ch/~lunghi/webpage/LatAves/page7/page7.html
  ##
  ## quant$quant.add.single("f_K_by_f_pi", 1.192, 0.005)
  ## quant$quant.add.single("f_K", 156.1, 1.1)

  ##--- http://arxiv.org/abs/1101.5138
  ## quant$quant.add.single("f_K_by_f_pi", 1.189, 0.007)
  ## quant$quant.add.single("f_K", 157, 2)

  ##+++ check effect of lattice correlations
  ## quant$corr.add.single("f_K_by_f_pi", "f_K", 100/100)

  ## ////////////////////////////////////////
  ##
  ## lepton universality tests
  ##

  ##
  ## gtau/gmu using tau -> hnu / h -> mu nu ---upd16
  ## PDG 2015
  ##
  quant$quant.add.single("pitoENu", 1.230e-4, 0.004e-4)
  quant$quant.add.single("pitoMuNu", 99.98770e-2, 0.00004e-2)
  quant$quant.add.single("KtoENu", 1.582e-5, 0.007e-5) ## changed16
  quant$quant.add.single("KtoMuNu", 63.56e-2, 0.11e-2) ## changed16

  ##--- gtau/gmu using tau -> pi nu / pi -> mu nu
  quant$quant.expr.add("gtaubygmu_pi",
                      sqrt(Gamma9/pitoMuNu *(2*m_pi*m_mu^2*tau_pi) /((1+dRrad_taupi_by_pimu)*m_tau^3*tau_tau)*
                           ((1-m_mu^2/m_pi^2)/(1-m_pi^2/m_tau^2))^2))

  ##--- gtau/gmu using tau -> K nu / K -> mu nu
  rc = quant$quant.expr.add(
    "gtaubygmu_K",
    sqrt(Gamma10/KtoMuNu *(2*m_K*m_mu^2*tau_K) / ((1+dRrad_tauK_by_Kmu)*m_tau^3*tau_tau) * ((1-m_mu^2/m_K^2)/(1-m_K^2/m_tau^2))^2)
  )

  ##--- gtau/gmu using tau lifetime
  quant$quant.expr.add("gtaubygmu_tau", sqrt(Gamma5/Be_from_taulife))

  ##--- gtau/gmu average
  quant$quant.fit.add("gtaubygmu_fit", c(gtaubygmu_tau=1, gtaubygmu_pi=1, gtaubygmu_K=1))

  ##--- gtau/ge using tau lifetime
  quant$quant.expr.add("gtaubyge_tau", sqrt(Gamma3/Bmu_from_taulife))

  ##--- gmu / ge from tau -> mu / tau -> e
  quant$quant.expr.add("gmubyge_tau", sqrt(Gamma3/Gamma5 * phspf_mebymtau/phspf_mmubymtau))

  ##--- add quantities to print
  display.names = c(
    display.names,
    "gtaubygmu_tau", "gtaubygmu_pi", "gtaubygmu_K", "gtaubygmu_fit",
    "gtaubyge_tau", "gmubyge_tau"
    )

  ## ////////////////////////////////////////
  ##
  ## Vus from tau -> Knu
  ##

  if (is.na(quant$val("Gamma10by9"))) {
    ##--- compute intermediate B_tau_K / B_tau_pi
    quant$quant.expr.add("Gamma10by9", Gamma10/Gamma9)
    quant$texdescr.add(c("Gamma10by9"="\\frac{\\BRF{\\tau^-}{K^- \\nu_\\tau}}{\\BRF{\\tau^-}{\\pi^- \\nu_\\tau}}"))
  }

  ##
  ## Vus^2 = Vud^2 * B(tau->Knu)/B(tau->pinu) * f_pi^2/f_K^2 * (1-m_pi^2/m_tau*2)^2/(1-m_K^2/m_tau*2)^2 / (1+dRrad_tauK_by_taupi)
  ##
  rc = quant$quant.expr.add(
    "Vus_tauKpi",
    Vud * sqrt(Gamma10/Gamma9) / f_K_by_f_pi * (1-m_pi^2/m_tau^2)/(1-m_K^2/m_tau^2) / sqrt(Rrad_tauK_by_taupi)
  )
  
  ##--- Vus_tauKpi vs Vus-from-CKM-unitarity
  quant$quant.expr.add("Vus_tauKpi_mism", Vus_tauKpi - Vus_uni)
  Vus_tauKpi_mism_sigma = quant$val("Vus_tauKpi_mism") / quant$err("Vus_tauKpi_mism")
  quant$param.add.single("Vus_tauKpi_mism_sigma", Vus_tauKpi_mism_sigma)
  quant$param.add.single("Vus_tauKpi_mism_sigma_abs", abs(Vus_tauKpi_mism_sigma))

  ##--- theory error contribution
  rc = quant$err.contrib.perc("Vus_tauKpi", "f_K_by_f_pi", "dRrad_taupi_by_pimu", "dRrad_tauK_by_Kmu", "dRrad_kmunu_by_pimunu")
  quant$param.add.single("Vus_tauKpi_err_th_perc", rc)

  lapply(c("f_K_by_f_pi", "dRrad_taupi_by_pimu", "dRrad_tauK_by_Kmu", "dRrad_kmunu_by_pimunu"), function(val) {
    quant$param.add.single(paste("Vus_tauKpi_err_th_perc", val, sep="_"),
                           quant$err.contrib.perc("Vus_tauKpi", val))
  })

  ## ////////////////////////////////////////
  ##
  ## Vus from tau -> K nu
  ##

  ##
  ## J. Erler, “Electroweak radiative corrections to semileptonic tau decays”,
  ## Rev. Mex. Fis. 50:200–202, 2004, arXiv:hep-ph/0211345.
  ## ---upd14
  ##
  quant$quant.add.single("Rrad_SEW_tau_Knu", 1.0201, 0.0003)

  ##
  ## R.Decker, M. Finkemeier, Nucl.Phys.Proc.Suppl. 40 (1995) 453-461
  ##
  quant$quant.add.single("dRrad_k_munu", 1.3e-2, 0.2e-2)  

  ##
  ## PDG 2015, citing CODATA 2014 Mohr:2015ccw ---upd16
  ## http://pdg.lbl.gov/2015/reviews/rpp2015-rev-phys-constants.pdf
  ##
  ## Mohr:2008fa, codata 2006, http://inspirehep.net/record/791091
  ## CODATA Recommended Values of the Fundamental Physical Constants: 2006.
  ##
  ## CODATA 2010 https://inspirehep.net/record/1095471/
  ## Mohr:2008fa CODATA Recommended Values of the Fundamental Physical Constants: 2010
  ##
  ## Plack constant h/2pi in MeV s
  ##
  quant$quant.add.single(
    "hcut",
    6.582119514e-22,
    0.000000040e-22
  )
  ## 2014: quant$quant.add.single("hcut", 6.58211928e-22, 0.00000015e-22)
  ## 2011: quant$quant.add.single("hcut", 6.58211899e-22, 0.00000016e-22)

  ##
  ## from PDG 2015, PDG 2013 ---upd16
  ## http://pdg.lbl.gov/2015/reviews/rpp2015-rev-phys-constants.pdf
  quant$quant.add.single("G_F_by_hcut3_c3", 1.1663787e-5*1e-6, 0.0000006e-5*1e-6)
  ## 2011: quant$quant.add.single("G_F_by_hcut3_c3", 1.16637e-5*1e-6, 1.16637e-5*1e-6 *9e3/1e9)

  ##
  ## Vus from tau -> K nu
  ##
  rc = quant$quant.expr.add("Vus_tauKnu",
    sqrt(Gamma10 * 16*pi * hcut / (m_tau^3*tau_tau*(1+dRrad_tauK_by_Kmu)*(1+dRrad_k_munu))) /
    (G_F_by_hcut3_c3 * f_K * (1 - m_K^2/m_tau^2)))

  ##--- Vus_tauKnu vs Vus-from-CKM-unitarity
  quant$quant.expr.add("Vus_tauKnu_mism", Vus_tauKnu - Vus_uni)
  Vus_tauKnu_mism_sigma = quant$val("Vus_tauKnu_mism") / quant$err("Vus_tauKnu_mism")
  quant$param.add.single("Vus_tauKnu_mism_sigma", Vus_tauKnu_mism_sigma)
  quant$param.add.single("Vus_tauKnu_mism_sigma_abs", abs(Vus_tauKnu_mism_sigma))

  ##--- theory error contribution
  quant$param.add.single(
          "Vus_tauKnu_err_th_perc",
          quant$err.contrib.perc("Vus_tauKnu", "f_K",
                                 ## "Rrad_SEW_tau_Knu",
                                 "dRrad_tauK_by_Kmu",
                                 "dRrad_k_munu"))

  ## ////////////////////////////////////////
  ##
  ## Vus from tau fit
  ##
  ## quant$quant.fit.add("Vus_tau", c(Vus=1, Vus_tauKpi=1, Vus_tauKnu=1))
  ##--- remove tau -> Knu for 2016 HFAG publication
  quant$quant.fit.add("Vus_tau", c(Vus=1, Vus_tauKpi=1))

  ##--- Vus_tau vs Vus-from-CKM-unitarity
  quant$quant.expr.add("Vus_tau_mism", Vus_tau - Vus_uni)
  Vus_tau_mism_sigma = quant$val("Vus_tau_mism") / quant$err("Vus_tau_mism")
  quant$param.add.single("Vus_tau_mism_sigma", Vus_tau_mism_sigma)
  quant$param.add.single("Vus_tau_mism_sigma_abs", abs(Vus_tau_mism_sigma))

  ##
  ## summary
  ##

  ##--- add quantities to print
  display.names = c(
    display.names,
    "Rrad_tauK_by_taupi",
    "Gamma10by9",
    "Vus_tauKpi",
    "Vus_tauKpi_mism_sigma",
    "Vus_tauKpi_err_th_perc",
    paste("Vus_tauKpi_err_th_perc", c("f_K_by_f_pi", "dRrad_taupi_by_pimu", "dRrad_tauK_by_Kmu", "dRrad_kmunu_by_pimunu"), sep="_"),
    "Rrad_SEW_tau_Knu",
    "Vus_tauKnu",
    "Vus_tauKnu_mism_sigma",
    "Vus_tauKnu_err_th_perc",
    "Vus_tau",
    "Vus_tau_mism_sigma",
    "Vud_moulson_ckm14",
    "Vus_kl3_moulson_ckm14",
    "VusbyVud_moulson_ckm14",
    "Vus_kl2_moulson_ckm14",
    "Vus_maltman_mainz16",
    "Vus_kl3_maltman_mainz16",
    "Vus_kl2_maltman_mainz16"
    )

  ##--- show selection of quantities
  print(round.data.frame(data.frame(val=quant$vals(display.names), err=quant$errs(display.names)), digits=6))

  ##--- translate quantity names to get a valid LaTeX command
  toTex = TrStr.num2tex$new()
  ##--- non-default formatting for selected quantities
  specialFormat = c(
    tau_tau="%.1f",
    Vud="%.5f",
    deltaR_su3break="%.3f",
    Be_lept="%.3f",
    Be_univ="%.3f",
    B_tau_had_fit="%.2f",
    B_tau_VA_fit="%.2f",
    B_tau_s_fit="%.3f",
    Vus_mism_sigma="%.1f",
    Vus_mism_sigma_abs="%.1f",
    Vus_tauKpi_mism_sigma="%.1f",
    Vus_tauKpi_mism_sigma_abs="%.1f",
    Vus_tauKnu_mism_sigma="%.1f",
    Vus_tauKnu_mism_sigma_abs="%.1f",
    Vus_tau_mism_sigma="%.1f",
    Vus_tau_mism_sigma_abs="%.1f",
    f_K="%.1f",
    dRrad_kmunu_by_pimunu="%.2f",
    dRrad_taupi_by_pimu="%.2f",
    dRrad_tauK_by_Kmu="%.2f",
    dRrad_k_munu="%.2f",
    Rrad_tauK_by_taupi="%.2f",
    Vus_err_th_perc="%.2f"
    )

  ##--- non-default multiplicative factor for selected quantities
  specialFactor = c(
    tau_tau=1e15,
    Be_lept=100,
    Be_univ=100,
    B_tau_had_fit=100,
    B_tau_VA_fit=100,
    B_tau_s_fit=100,
    dRrad_kmunu_by_pimunu=100,
    dRrad_taupi_by_pimu=100,
    dRrad_tauK_by_Kmu=100,
    dRrad_k_munu=100
    )

  ##--- provide LaTeX commands printing quantity +- error
  quant.all.sorted = quant$vnames(order(quant$vnames()))

  rc = mapply(function(name, val, err, texdescr) {                
    ##--- special rounding for GammAll and Gamma998, hide numeric approximations
    if (name == "GammaAll") {
      val = round(val, digits=12)
      err = round(err, digits=12)
    }
    if (name == "Gamma998") {
      val = round(val, digits=12)
      err = round(err, digits=12)
    }

    if (!is.na(specialFactor[name])) {
      val = val*specialFactor[name]
      err = err*specialFactor[name]
    }
    if (is.na(specialFormat[name])) {
      rc = alurep.precision.order(c(val, err))
      precision = rc$precision
      order = rc$order
      val.err = alurep.tex.val.err.prec.ord(val, err, precision, order)
      val.str = alurep.tex.val.prec.ord(val, precision, order)
      err.str = alurep.tex.val.prec.ord(err, precision, order)
    } else {
      fmt = specialFormat[name]
      val.str = sprintf(fmt, val)
      err.str = sprintf(fmt, err)
      val.err = paste(val.str, "\\pm", err.str)
    }
    if (err == 0) {
      err.str = "0"
      val.err = val.str
    }
    gammaname = alurep.gamma.texlabel(name)
    if (is.na(texdescr)) {texdescr = ""}
    rc = paste("\\htquantdef{", name, "}{",
      gammaname, "}{", texdescr, "}{", val.err, "}{", val.str, "}{", err.str, "}%", sep="")
  },
    quant.all.sorted,
    quant$vals(quant.all.sorted),
    quant$errs(quant.all.sorted),
    quant$texdescr(quant.all.sorted))

  ##--- print correlation of universality results
  quant.list = c(
    "gtaubygmu_tau", "gtaubyge_tau", "gmubyge_tau",
    "gtaubygmu_pi", "gtaubygmu_K"
    )
  univ.corr = quant$corr(quant.list, quant.list)

  ##
  ## Correlation matrix has one eigenvalue = 0 because mu/e is 100%
  ## correlated with tau/mu and tau/e. To avoid getting a non positive
  ## semi-definite matrix we regularize if the minimum eigenvalues is
  ## negative
  min.eigenvalue.univ.corr = min(eigen(univ.corr)$values)
  if (min.eigenvalue.univ.corr < 0) {
    cat(file=stderr(), "minimum eigenvalue of couplings ratios correlation =",
        min.eigenvalue.univ.corr)
    cat("\n")
    cat(file=stderr(), "correlation will be regularized to be positive semi-definite\n")
    univ.corr = nearPD(univ.corr, corr=TRUE, posd.tol=0)$mat
  }
  
  univ.corr = 100*univ.corr
  ## print(univ.corr)

  tex.corr = NULL

  tex.corr = c(tex.corr,
    paste("$\\left( \\frac{g_\\tau}{g_e} \\right)$",
          paste(sprintf("%4.0f", univ.corr[2, 1]), collapse=" & "),
          sep=" & "))

  tex.corr = c(tex.corr,
    paste("$\\left( \\frac{g_\\mu}{g_e} \\right)$",
          paste(sprintf("%4.0f", univ.corr[3, 1:2]), collapse=" & "),
          sep=" & "))

  tex.corr = c(tex.corr,
    paste("$\\left( \\frac{g_\\tau}{g_\\mu} \\right)_\\pi$",
          paste(sprintf("%4.0f", univ.corr[4, 1:3]), collapse=" & "),
          sep=" & "))

  tex.corr = c(tex.corr,
    paste("$\\left( \\frac{g_\\tau}{g_\\mu} \\right)_K$",
          paste(sprintf("%4.0f", univ.corr[5, 1:4]), collapse=" & "),
          sep=" & "))

  tex.corr = c(tex.corr, paste(
    "",
    "$\\left( \\frac{g_\\tau}{g_\\mu} \\right)$",
    "$\\left( \\frac{g_\\tau}{g_e} \\right)$",
    "$\\left( \\frac{g_\\mu}{g_e} \\right)$",
    "$\\left( \\frac{g_\\tau}{g_\\mu} \\right)_\\pi$",
    ## "$\\left( \\frac{g_\\tau}{g_\\mu} \\right)_K$",
     sep=" & ", collapse=" & "))

  rc = c(rc, alurep.tex.cmd("couplingsCorr", paste(tex.corr, collapse="\\\\\n")))

  ##--- assemble file name
  fname = file.name
  fname = basename(fname)
  fname = sub("[.][^.]*$", "", fname, perl=TRUE)
  fname = sub("average[^-]*-*", "", fname)
  if (fname != "") fname = paste0("-", fname)
  fname = file.path("../report", paste0("tau-br-fit", fname))

  fname = paste0(fname, "-elab")
  if (flag.unitarity) fname = paste0(fname, "-uniconstr")
  if (flag.vadirect) fname = paste0(fname, "-vadirect")
  if (flag.lepuniv) fname = paste0(fname, "-lepuniv")
  if (flag.kmaltman) fname = paste0(fname, "-kmaltman")

  fname = paste0(fname, ".tex")
  cat(unlist(rc), sep="\n", file=fname)
  cat("produced file '", fname, "'\n", sep="")

  fname = gsub("[.]tex$", ".rdata", basename(fname))
  rc = save(file=fname, list = c("quant"))
  cat("produced file '", fname, "'\n", sep="")
}

args = commandArgs(TRUE)
aluelab.results(args)
