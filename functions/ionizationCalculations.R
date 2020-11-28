#####log D
#reference: Lipophilicity in Drug Action and Toxicology Applications of a Solvation Equation to Drug Transport Properties M. H. Abraham H. S. Chadha
f.logD <- function(ionization, pH, pKa1, pKa2, logP){
  if (ionization=="neutral") {
    logD<-logP
  } else if (ionization =="monoproticAcid") {
    logD<-logP-log10(1+10^(pH-pKa1))
  } else if (ionization =="monoproticBase") {
    logD<-logP-log10(1+10^(pKa1-pH))
  } else if (ionization == "diproticAcid") {
    logD<-logP-log10(1+10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2))
  } else if (ionization == "diproticBase") {
    logD<-logP-log10(1+10^(pKa2-pH)+10^(pKa2-pH+pKa1-pH))
  } else {
    logD<-logP-log10(1+10^(pH-pKa1)+10^(pKa2-pH))
  }
  return(logD) 
}

#####fractions Ionized (acid/base) vs unionized

#reference: Lipophilicity in Drug Action and Toxicology Applications of a Solvation Equation to Drug Transport Properties M. H. Abraham H. S. Chadha
f.fracUnionized <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="neutral") {
    fracUnionized<-1
  } else if (ionization =="monoproticAcid") {
    fracUnionized<-1/(1+10^(pH-pKa1))
  } else if (ionization =="monoproticBase") {
    fracUnionized<-1/(1+10^(pKa1-pH))
  } else if (ionization == "diproticAcid") {
    fracUnionized<-1/(1+10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2))
  } else if (ionization == "diproticBase") {
    fracUnionized<-1/(1+10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH))
  } else {
    fracUnionized<-1/(1+10^(pH-pKa1)+10^(pKa2-pH))
  }
  return(fracUnionized) 
}


f.fracIonizedAcid <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="zwitterionic"||ionization=="monoproticAcid") {
    fracIonizedAcid<-1-1/(1+10^(pH-pKa1))
  } else if (ionization=="diproticAcid"){
    fracIonizedAcid<-1-1/(1+10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2))
  } else {
    fracIonizedAcid<-0
  }
  return(fracIonizedAcid) 
}

f.fracIonizedBase <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="zwitterionic") {
    fracIonizedBase<-1-1/(1+10^(pKa2-pH))
  } else if (ionization=="monoproticBase") {
    fracIonizedBase<-1-1/(1+10^(pKa1-pH))
  } else if (ionization=="diproticBase") {
    fracIonizedBase<-1-1/(1+10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH))
  } else {
    fracIonizedBase<-0
  }
  return(fracIonizedBase) 
}

