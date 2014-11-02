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
    cat("  -lepuniv use R_had = (1 - (1+(f_mu/f_e))Be_univ)/Be_univ\n")
    args = args[args != "-s"]
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

  if (length(args) > 0) {
    file.name = args[1]
  } else {
    file.name = "average.rdata"
  }

  ##--- get alucomb results and data
  load(file.name)

  if (flag.s.factors) {
    quant.err = quant.sf.err
    quant.cov = quant.sf.cov
  }

  quant = StatComb$new(quant.val, quant.cov)

  quant.names = names(quant.val)
  comb.params = lapply(combination$params, function(x) unname(x["value"]))
  quant$param.add(comb.params)
  quant$param.add(c(pi=pi))

  quant.texdescr = sapply(combination$quantities, function(x) {if (is.null(x$texdescr)) {""} else {x$texdescr}})
  quant.descr = sapply(combination$quantities, function(x) {if (is.null(x$descr)) {""} else {x$descr}})

  quant$texdescr.add(alurep.get.texdescr(quant.descr, quant.texdescr))

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
  B.tau.VA.names = setdiff(aluelab.get.quant.names("GammaAll", combination), Gamma110.names)
  B.tau.VA.names = setdiff(B.tau.VA.names, c("Gamma5", "Gamma3", "Gamma998"))
  quant$quant.qexpr.add("B_tau_VA", parse(text=paste(B.tau.VA.names, collapse="+")))
  quant$quant.expr.add("B_tau_VA_unitarity", 1-Gamma5-Gamma3-Gamma110)
  quant$quant.expr.add("B_tau_s_unitarity", 1-Gamma5-Gamma3-B_tau_VA)

  ##
  ## add measurements to compute universality improved Be
  ##

  ##--- from PDG 2013 ---upd14
  quant$quant.add.single("m_e",0.510998928, 0.000000011) ## changed14
  quant$quant.add.single("m_mu", 105.6583715, 0.0000035) ## changed14

  ## quant$quant.add.single("tau_tau", 290.6e-15, 1.0e-15)
  ##--- use HFAG 2014 average for tau lifetime (= PDG14)
  quant$quant.add.single("tau_tau", 290.29e-15, 0.52e-15) ## changed14

  ##--- m_tau HFAG 2009
  ## quant$quant.add.single("m_tau", 1776.7673082, 0.1507259)
  ##--- m_tau PDG 2013 ---upd14
  quant$quant.add.single("m_tau", 1776.82, 0.16)

  ##---upd14
  quant$quant.add.single("m_pi", 139.57018, 0.00035)
  quant$quant.add.single("tau_pi", 2.6033e-8, 0.0005e-8)
  quant$quant.add.single("m_K", 493.677, 0.016)
  quant$quant.add.single("tau_K", 1.2380e-8, 0.0021e-8)

  ##--- from PDG 2013, 2011 ---upd14
  quant$quant.add.single("m_W", 80.399385e3, 0.015*1e3) ## changed14
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
  quant$quant.expr.add("Be_from_Bmu", phspf_mebymtau/phspf_mmubymtau *Bmu_fit)

  ##
  ## rad. corrections from to get Be from tau lifetime
  ## values from 10.1103/RevModPhys.78.1043 p.1047, arXiv:hep-ph/0507078v2 p.7, could be recomputed
  ## - delta^L_gamma = 1 + alpha(mL)/2pi * (25/4 - pi^2)
  ## - delta^L_W = 1 + 3/5* m_L^2/M_W^2
  ## ---upd14
  ## - 2014: no update available
  quant$param.add(c(delta_mu_gamma=(1 - 42.4e-4), delta_tau_gamma=(1 - 43.2e-4)))
  quant$quant.expr.add("delta_mu_W", 1 + 3/5*m_mu^2/m_W^2)
  quant$quant.expr.add("delta_tau_W", 1 + 3/5*m_tau^2/m_W^2)

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
  ## universality improved Be = B(tau -> e nu nubar (gamma))
  ## see arXiv:hep-ph/0507078v2 p.7, doi:10.1103/RevModPhys.78.1043 p.1047
  ##
  ## minimum chisq fit for Be_univ using, Be, Be from Bmu, Be from tau lifetime
  ##
  ## Bmu/Be = f(m_mu^2/m_tau^2) / f(m_e^2/m_tau^2)
  ## Be= tau_tau / tau_mu (m_tau/m_mu)^5 f(m^2_e/m^2_tau)/f(m^2_e/m^2_mu) (delta^tau_gamma delta^tau_W)/(delta^mu_gamma delta^mu_W)
  ##
  quant$quant.fit.add("Be_univ", c(Gamma5=1, Be_from_Bmu=1, Be_from_taulife=1))

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
  quant$quant.add.single("Vud", 0.97425, 0.00001*sqrt(8^2 + 10^2 + 18^2))

  ## values of Moulson CKM 2014
  quant$quant.add.single("Vud_moulson_ckm14", 0.97417, 0.00001*21)
  quant$quant.add.single("Vus_kl3_moulson_ckm14", 0.2232, 0.0001*9)
  quant$quant.add.single("VusbyVud_moulson_ckm14", 0.2308, 0.0001*6)
  quant$quant.expr.add("Vus_kl2_moulson_ckm14", VusbyVud_moulson_ckm14 * Vud_moulson_ckm14)

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

  ##--- s quark mass, PDG2011
  ## quant$quant.add.single("m_s", 100, sqrt((20.^2 + 30.^2)/2.))
  ##--- s quark mass, PDG2013 ---upd14
  quant$quant.add.single("m_s", 93.5, 2.5) ## changed14

  ##--- E.Gamiz, M.Jamin, A.Pich, J.Prades, F.Schwab, |V_us| and m_s from hadronic tau decays
  ## quant$quant.add.single("deltaR_su3break_pheno", 0.1544, 0.0037)
  quant$quant.add.single("deltaR_su3break_d2pert", 9.3, 3.4)
  quant$quant.add.single("deltaR_su3break_remain", 0.0034, 0.0028)

  ## 
  ## K.Maltman, private comm. 2014
  ## - use tau pi Kl2 poles
  ## - FOPT 4 order
  ##
  quant$quant.add.single("deltaR_su3break_pheno", 0.1509, 0.0041) 
  quant$quant.add.single("deltaR_su3break_d2pert", 11.36, 1.14)
  quant$quant.add.single("deltaR_su3break_remain", 0.0084, 0.0185)
  
  ## quant$quant.expr.add("deltaR_su3break", deltaR_su3break_pheno + deltaR_su3break_d2pert*(m_s/1000)^2 + deltaR_su3break_remain)

  ## 
  ## K.Maltman, private comm. 2014
  ## - use FOPT 3rd loop order as central value (favoured by Lattice)
  ## - use FOPT uncertainty
  ## - add FOPT-CIPT 3rd order discrepancy as additional theory uncertainty
  ##
  quant$quant.add.single("deltaR_su3break", 0.254, 0.038)
  
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

  ##--- add Vus
  quant$quant.expr.add("Vus", sqrt(R_tau_s/(R_tau_VA/Vud^2 - deltaR_su3break)))
  quant$quant.add.single("Vus_err_perc", quant$err.contrib.perc("Vus"), 0)

  ##--- theory, exp errors, with percent
  quant$quant.add.single("Vus_err_th", quant$err.contrib("Vus", "deltaR_su3break"), 0)
  quant$quant.add.single("Vus_err_th_perc", quant$err.contrib.perc("Vus", "deltaR_su3break"), 0)

  Vus.err.exp = sqrt(quant$err("Vus")^2 - quant$param("Vus_err_th")^2)
  quant$quant.add.single("Vus_err_exp", Vus.err.exp, 0)
  quant$quant.add.single("Vus_err_exp_perc", Vus.err.exp/quant$val("Vus")*100, 0)

  ##--- Vus from Vud using CKM unitarity
  quant$quant.expr.add("Vus_uni", sqrt(1-Vud^2))

  ##--- Vus vs <Vus from CKM unitarity>
  quant$quant.expr.add("Vus_mism", Vus - Vus_uni)
  Vus_mism_sigma = quant$val("Vus_mism") / quant$err("Vus_mism")
  quant$quant.add.single("Vus_mism_sigma", Vus_mism_sigma, 0)
  quant$quant.add.single("Vus_mism_sigma_abs", abs(Vus_mism_sigma), 0)

  ##--- select quantities to print
  display.names = c(
    ## Gamma110.names,
    "Gamma5", "Be_unitarity", "Be_fit",
    "Gamma3", "Bmu_unitarity", "Bmu_fit",
    "Bmu_by_Be_th", "Be_from_Bmu", "Be_from_taulife", "Be_univ",
    "Bmu_from_taulife",
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

  ##
  ## gtau/gmu using tau -> hnu / h -> mu nu ---upd14
  ## PDG 2013
  ##
  quant$quant.add.single("pitoENu", 1.230e-4, 0.004e-4)
  quant$quant.add.single("pitoMuNu", 99.98770e-2, 0.00004e-2)
  quant$quant.add.single("KtoENu", 1.581e-5, 0.008e-5) ## changed14
  quant$quant.add.single("KtoMuNu", 63.55e-2, 0.11e-2)

  ##--- from Marciano:1993sh,Decker:1994ea,Decker:1994dd ---upd14
  quant$quant.add.single("delta_pi", 0.16e-2, 0.14e-2)
  quant$quant.add.single("delta_K", 0.90e-2, 0.22e-2)

  ##--- gtau/gmu using tau -> pi nu / pi -> mu nu
  quant$quant.expr.add("gtaubygmu_pi",
                      sqrt(Gamma9/pitoMuNu *(2*m_pi*m_mu^2*tau_pi) /((1+delta_pi)*m_tau^3*tau_tau)*
                           ((1-m_mu^2/m_pi^2)/(1-m_pi^2/m_tau^2))^2))

  ##--- gtau/gmu using tau -> K nu / K -> mu nu
  quant$quant.expr.add("gtaubygmu_K",
                      sqrt(Gamma10/KtoMuNu *(2*m_K*m_mu^2*tau_K) /((1+delta_K)*m_tau^3*tau_tau)*
                           ((1-m_mu^2/m_K^2)/(1-m_K^2/m_tau^2))^2))

  ##--- gtau/gmu using tau lifetime
  quant$quant.expr.add("gtaubygmu_tau", sqrt(Gamma5/Be_from_taulife))

  ##--- gtau/gmu average
  quant$quant.fit.add("gtaubygmu_fit", c(gtaubygmu_tau=1, gtaubygmu_pi=1, gtaubygmu_K=1))

  ##--- gtau/ge using tau lifetime
  quant$quant.expr.add("gtaubyge_tau", sqrt(Gamma3/Bmu_from_taulife))

  ##--- gmu / ge from tau -> mu / tau -> e
  quant$quant.expr.add("gmubyge_tau", sqrt(Gamma3/Gamma5 * phspf_mebymtau/phspf_mmubymtau))

  ##--- select quantities to print
  display.names = c(
    display.names,
    "gtaubygmu_tau", "gtaubygmu_pi", "gtaubygmu_K", "gtaubygmu_fit",
    "gtaubyge_tau", "gmubyge_tau"
    )

  ## ////////////////////////////////////////
  ##
  ## Vus from tau -> Knu
  ##

  ##
  ## QCD lattice
  ##
  ## ---upd14
  ## Lattice averages FLAG 2013 Aoki et al.
  ## http://arxiv.org/abs/1310.8555
  ## http://inspirehep.net/record/1262813, Aoki:2013ldr
  ## use values for Nf = 2+1 (compatible with Nf = 2+1+1
  ##
  quant$quant.add.single("f_K_by_f_pi", 1.192, 0.005) ## unchanged12
  quant$quant.add.single("f_K", 156.3, 0.9) ## changed14
  quant$quant.add.single("fp0_Kpi", 0.9661, 0.0032) ## new since 2014

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

  ##
  ## Marciano:2004uf
  ## W. J. Marciano, "Precise determination of |V(us)| from lattice calculations of pseudoscalar decay constants",
  ## Phys. Rev. Lett. 93:231803, 2004, doi:10.1103/PhysRevLett.93.231803, arXiv:hep-ph/0402299.
  ## ---upd14
  ##
  quant$quant.add.single("rrad_LD_kmu_pimu", 0.9930, 0.0035)

  ##
  ## Decker:1994ea
  ## R. Decker and M. Finkemeier, "Short and long distance effects in the decay tau -> pi nu_tau (gamma)",
  ## Nucl. Phys. B438:17-53, 1995, doi:10.1016/0550-3213(95)00597-L, arXiv:hep-ph/9403385.
  ## delta_LD(tau -> h nu / h -> mu nu)
  ## ---upd14
  ##
  quant$quant.add.single("delta_LD_taupi_pimu", 0.16/100, 0.14/100)
  quant$quant.add.single("delta_LD_tauK_Kmu", 0.90/100, 0.22/100)

  ##--- compute delta_LD_tauK_taupi = delta_LD_tauK_Kmu/delta_LD_taupi_pimu * delta_LD_kmu_pimu
  quant$quant.expr.add("rrad_LD_tauK_taupi", (1+delta_LD_tauK_Kmu)/(1+delta_LD_taupi_pimu) * rrad_LD_kmu_pimu)

  ##--- compute intermediate B_tau_K / B_tau_pi
  quant$quant.expr.add("Gamma10by9", Gamma10/Gamma9);

  ##
  ## Vus^2 = Vud^2 * B(tau -> Knu)/B(tau -> pinu) * f_pi^2/f_K^2 * (1-m_pi^2/m_tau*2)/(1-m_K^2/m_tau*2) /(1+delta_LD)
  ##
  quant$quant.expr.add("Vus_tauKpi", Vud*sqrt(Gamma10/Gamma9) /f_K_by_f_pi
                      * (1-m_pi^2/m_tau^2)/(1-m_K^2/m_tau^2) / sqrt(rrad_LD_tauK_taupi))

  ##--- Vus_tauKpi vs Vus-from-CKM-unitarity
  quant$quant.expr.add("Vus_tauKpi_mism", Vus_tauKpi - Vus_uni)
  Vus_tauKpi_mism_sigma = quant$val("Vus_tauKpi_mism") / quant$err("Vus_tauKpi_mism")
  quant$quant.add.single("Vus_tauKpi_mism_sigma", Vus_tauKpi_mism_sigma, 0)
  quant$quant.add.single("Vus_tauKpi_mism_sigma_abs", abs(Vus_tauKpi_mism_sigma), 0)

  ##--- theory error contribution
  rc = quant$err.contrib.perc("Vus_tauKpi", "f_K_by_f_pi", "delta_LD_taupi_pimu", "delta_LD_tauK_Kmu", "rrad_LD_kmu_pimu")
  quant$quant.add.single("Vus_tauKpi_err_th_perc", rc, 0)

  lapply(c("f_K_by_f_pi", "delta_LD_taupi_pimu", "delta_LD_tauK_Kmu", "rrad_LD_kmu_pimu"), function(val) {
    quant$quant.add.single(paste("Vus_tauKpi_err_th_perc", val, sep="_"),
                          quant$err.contrib.perc("Vus_tauKpi", val), 0)
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
  quant$quant.add.single("rrad_tau_Knu", 1.0201, 0.0003)

  ##
  ## Mohr:2008fa, codata 2006, http://inspirehep.net/record/791091
  ## CODATA Recommended Values of the Fundamental Physical Constants: 2006.
  ##
  ## CODATA 2010 https://inspirehep.net/record/1095471/
  ## Mohr:2008fa CODATA Recommended Values of the Fundamental Physical Constants: 2010
  ##
  ## --- Plack h/ in MeV s
  ## h/ = 6.58211928(15)e-16 eV s
  quant$quant.add.single("hcut", 6.58211928e-22, 0.00000015e-22)
  ## 2011:quant$quant.add.single("hcut", 6.58211899e-22, 0.00000016e-22)

  ## --- G_F / (hcut c)^3 from PGD11 [1.1663787(6)] in GeV^-2, converted to MeV^-2 ---PDG13
  quant$quant.add.single("G_F_by_hcut3_c3", 1.1663787e-5*1e-6, 0.0000006e-5*1e-6)
  ## 2011: quant$quant.add.single("G_F_by_hcut3_c3", 1.16637e-5*1e-6, 1.16637e-5*1e-6 *9e3/1e9)

  ##
  ## Vus from tau -> K nu
  ##
  rc = quant$quant.expr.add("Vus_tauKnu",
    sqrt(Gamma10 * 16*pi * hcut / (m_tau^3*tau_tau*rrad_tau_Knu)) / (G_F_by_hcut3_c3 * f_K * (1 - m_K^2/m_tau^2)))

  ##--- Vus_tauKnu vs Vus-from-CKM-unitarity
  quant$quant.expr.add("Vus_tauKnu_mism", Vus_tauKnu - Vus_uni)
  Vus_tauKnu_mism_sigma = quant$val("Vus_tauKnu_mism") / quant$err("Vus_tauKnu_mism")
  quant$quant.add.single("Vus_tauKnu_mism_sigma", Vus_tauKnu_mism_sigma, 0)
  quant$quant.add.single("Vus_tauKnu_mism_sigma_abs", abs(Vus_tauKnu_mism_sigma), 0)

  ##--- theory error contribution
  quant$quant.add.single("Vus_tauKnu_err_th_perc", quant$err.contrib.perc("Vus_tauKnu", "f_K", "rrad_tau_Knu"), 0)

  ## ////////////////////////////////////////
  ##
  ## Vus from tau fit
  ##
  quant$quant.fit.add("Vus_tau", c(Vus=1, Vus_tauKpi=1, Vus_tauKnu=1))

  ##--- Vus_tau vs Vus-from-CKM-unitarity
  quant$quant.expr.add("Vus_tau_mism", Vus_tau - Vus_uni)
  Vus_tau_mism_sigma = quant$val("Vus_tau_mism") / quant$err("Vus_tau_mism")
  quant$quant.add.single("Vus_tau_mism_sigma", Vus_tau_mism_sigma, 0)
  quant$quant.add.single("Vus_tau_mism_sigma_abs", abs(Vus_tau_mism_sigma), 0)

  ##
  ## summary
  ##
  display.names = c(
    display.names,
    "rrad_LD_tauK_taupi",
    "Gamma10by9",
    "Vus_tauKpi",
    "Vus_tauKpi_mism_sigma",
    "Vus_tauKpi_err_th_perc",
    paste("Vus_tauKpi_err_th_perc", c("f_K_by_f_pi", "delta_LD_taupi_pimu", "delta_LD_tauK_Kmu", "rrad_LD_kmu_pimu"), sep="_"),
    "rrad_tau_Knu",
    "Vus_tauKnu",
    "Vus_tauKnu_mism_sigma",
    "Vus_tauKnu_err_th_perc",
    "Vus_tau",
    "Vus_tau_mism_sigma",
    "Vud_moulson_ckm14",
    "Vus_kl3_moulson_ckm14",
    "VusbyVud_moulson_ckm14",
    "Vus_kl2_moulson_ckm14"
    )

  ##--- show selection of quantities
  print(rbind(cbind(val=quant$vals(display.names), err=quant$errs(display.names))))

  ##--- translate quantity names to get a valid LaTeX command
  toTex = TrStr.num2tex$new()
  ##--- non-default formatting for selected quantities
  specialFormat = c(
    tau_tau="%.1f",
    Vud="%.5f",
    deltaR_su3break="%.3f",
    Be_univ="%.3f",
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
    delta_LD_taupi_pimu="%.2f",
    delta_LD_tauK_Kmu="%.2f",
    delta_LD_tauK_taupi="%.2f"
    )
  ##--- non-default multiplicative factor for selected quantities
  specialFactor = c(
    tau_tau=1e15,
    Be_univ=100,
    B_tau_VA_fit=100,
    B_tau_s_fit=100,
    delta_LD_taupi_pimu=100,
    delta_LD_tauK_Kmu=100,
    delta_LD_tauK_taupi=100
    )

  ##--- provide LaTeX commands printing quantity +- error
  quant.all.sorted = quant$vnames(order(quant$vnames()))

  rc = mapply(function(name, val, err, texdescr) {
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

  rc = c(rc, paste("\\newcommand{\\couplingsCorr}{", paste(tex.corr, collapse="\\\\\n"), "}", sep=""))

  ##--- make file name for report .tex output
  fname.short = gsub("average[^-]*-*", "", file.name)
  fname.short = gsub("[.]rdata", "", fname.short)
  fname = "../report/tau-br-fit"
  if (fname.short != "") fname = paste(fname, "-", fname.short, sep="")
  if (flag.unitarity) fname = paste(fname, "-uniconstr", sep="")
  if (flag.vadirect) fname = paste(fname, "-vadirect", sep="")
  if (flag.lepuniv) fname = paste(fname, "-lepuniv", sep="")
  fname = paste(fname, "-elab.tex", sep="")
  cat(unlist(rc), sep="\n", file=fname)
  cat("produced file '", fname, "'\n", sep="")
}

args = commandArgs(TRUE)
aluelab.results(args)
