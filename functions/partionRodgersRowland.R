##### partition coefficients according to the method of Rodgers and Rowland 2006. doi:10.1002/jps.20502
f.XRodgersRowland <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="neutral") {
    XRodgersRowland<-0
  } else if (ionization =="monoproticAcid") {
    XRodgersRowland<-10^(pH-pKa1)
  } else if (ionization =="monoproticBase") {
    XRodgersRowland<-10^(pKa1-pH)
  } else if (ionization == "diproticAcid") {
    XRodgersRowland<-10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2)
  } else if (ionization == "diproticBase") {
    XRodgersRowland<-10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH)
  } else {
    XRodgersRowland<-10^(pH-pKa1)+10^(pKa2-pH)
  }
  return(XRodgersRowland) 
}


f.YRodgersRowland <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="neutral") {
    YRodgersRowland<-0
  } else if (ionization =="monoproticAcid") {
    YRodgersRowland<-10^(pH-pKa1)
  } else if (ionization =="monoproticBase") {
    YRodgersRowland<-10^(pKa1-pH)
  } else if (ionization == "diproticAcid") {
    YRodgersRowland<-10^(pH-pKa1)+10^(pH-pKa1+pH-pKa2)
  } else if (ionization == "diproticBase") {
    YRodgersRowland<-10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH)
  } else {
    YRodgersRowland<-10^(pH-pKa1)+10^(pKa2-pH)
  }
  return(YRodgersRowland) 
}

f.ZRodgersRowland <- function(ionization, pH, pKa1, pKa2){
  if (ionization=="monoproticBase" & pKa1>=7) {
    ZRodgersRowland<-10^(pKa1-pH)
  } else if (ionization == "diproticBase" & pKa1>=7) {
    ZRodgersRowland<-10^(pKa1-pH)+10^(pKa2-pH+pKa1-pH)
  } else if (ionization == "zwitterionic"&& pKa2>=7) {
    ZRodgersRowland<-10^(pH-pKa1)+10^(pKa2-pH)
  } else {
    ZRodgersRowland<-NA
  }
  return(ZRodgersRowland) 
}

f.KaPRRodgersRowland <- function(ionization, pKa1, pKa2, logP, BP, Fup, 
                                 YRodgersRowland, ZRodgersRowland, 
                                 Bfnl, Bfnp, Bfiw, APbc, H){
  KaPRRodgersRowlandcc<- ifelse(ionization == "monoproticBase"&& pKa1>=7 || ionization == "diproticBase"&& pKa1>=7|| ionization == "zwitterionic"&& pKa2>=7, 
                              ((BP+H-1)/(H*Fup) -
                                 (1+ZRodgersRowland)/(1+YRodgersRowland)*Bfiw-
                                 (10^logP*Bfnl+(0.3*10^logP+0.7)*Bfnp)/(1+YRodgersRowland))*
                                (1+YRodgersRowland)/(APbc*ZRodgersRowland)
                              , NA)
  return(ifelse(KaPRRodgersRowlandcc<0, 0, KaPRRodgersRowlandcc)) 
}


f.adipose<-function(ionization, XRodgersRowland,YRodgersRowland, ZRodgersRowland,KaPRRodgersRowland, 
                    logP, pH ,pKa1, pKa2, Fup, BP, fnl,
                    fnp,	fiw, few, albuminRatio, lipoproteinRatio,
                    APT, Pfnl,Pfnp,	Pfiw,
                    APbc,	H){
  if (ionization == "monoproticAcid" ) {
    logDvow<- 1.115*logP-1.35-log10(1+10^(pH-pKa1))
  } else if (ionization == "monoproticBase") {
    logDvow<- 1.115*logP-1.35-log10(1+10^(pKa1-pH))  
  } else if (ionization == "diproticAcid") {
    logDvow<- 1.115*logP-1.35-log10(1+10^(pH-pKa1+pH-pKa2))
  } else if (ionization == "diproticBase") {
    logDvow<- 1.115*logP-1.35-log10(1+10^(pKa1-pH+pKa2-pH))
  } else if (ionization == "zwitterionic") {
    logDvow<- 1.115*logP-1.35-log10(1+10^(-pKa2+pH+pKa1-pH))
  } else
    logDvow<- 1.115*logP-1.35
  
  
  if (ionization == "monoproticBase"&& pKa1>=7 || ionization == "diproticBase"&& pKa1>=7|| ionization == "zwitterionic"&& pKa2>=7) {
    Kp<- (1+XRodgersRowland)/(1+YRodgersRowland)*fiw+few+
      (10^logDvow*fnl +(0.3*10^logDvow+0.7)*fnp)/(1+YRodgersRowland)+
      (KaPRRodgersRowland*APT*XRodgersRowland)/(1+YRodgersRowland)
  } else if (ionization == "neutral") {
    Kp<-(1+XRodgersRowland)/(1+YRodgersRowland)*fiw+few+
      (10^logDvow*fnl +(0.3*10^logDvow+0.7)*fnp)/(1+YRodgersRowland)+
      (1/Fup-1-(10^logDvow*Pfnl+(0.3*10^logDvow+0.7)*Pfnp)/(1+YRodgersRowland))*lipoproteinRatio
  } else {
    Kp<-(1+XRodgersRowland)/(1+YRodgersRowland)*fiw+few+
      (10^logDvow*fnl +(0.3*10^logDvow+0.7)*fnp)/(1+YRodgersRowland)+
      (1/Fup-1-(10^logDvow*Pfnl+(0.3*10^logDvow+0.7)*Pfnp)/(1+YRodgersRowland))*albuminRatio
  }
  return(Kp*Fup)
}

f.rest<-function(ionization, XRodgersRowland,YRodgersRowland, ZRodgersRowland,KaPRRodgersRowland, 
                 logP, pH ,pKa1, pKa2, Fup, BP, fnl,
                 fnp,	fiw, few, albuminRatio, lipoproteinRatio,
                 APT, Pfnl,Pfnp,	Pfiw,
                 APbc,	H){
  if (ionization == "monoproticBase"&& pKa1>=7 ||ionization ==  "diproticBase"&& pKa1>=7||ionization == "zwitterionic"&& pKa2>=7) {
    Kp<- (1+XRodgersRowland)/(1+YRodgersRowland)*fiw+few+
      (10^logP*fnl +(0.3*10^logP+0.7)*fnp)/(1+YRodgersRowland)+
      (KaPRRodgersRowland*APT*XRodgersRowland)/(1+YRodgersRowland)
   } else if (ionization == "neutral") {
    Kp<-(1+XRodgersRowland)/(1+YRodgersRowland)*fiw+few+
      (10^logP*fnl +(0.3*10^logP+0.7)*fnp)/(1+YRodgersRowland)+
      (1/Fup-1-(10^logP*Pfnl+(0.3*10^logP+0.7)*Pfnp)/(1+YRodgersRowland))*lipoproteinRatio
  } else {
    Kp<-(1+XRodgersRowland)/(1+YRodgersRowland)*fiw+
      few+
      (10^logP*fnl +(0.3*10^logP+0.7)*fnp)/(1+YRodgersRowland)+
      (1/Fup-1-(10^logP*Pfnl+(0.3*10^logP+0.7)*Pfnp)/(1+YRodgersRowland))*albuminRatio
  }
  return(Kp*Fup)
}











  
  
  




