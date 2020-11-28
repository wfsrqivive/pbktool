##### partition coefficients according to the method of Berezhkovskiy 2004. doi: 10.1002/jps.20073.  

f.adiposeBerezhkovskiy<-function(ionization, pH, pKa1, pKa2, logP, Fup, Fut, Vnl, Vph, Vw,PVnl, PVph, PVw){
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
  
  Kp<-(10^logDvow*(Vnl+0.3*Vph)+0.7*Vph+Vw/Fut)/(10^logDvow*(PVnl+0.3*PVph)+0.7*PVph+PVw/Fup)
  return(Kp)
}



f.restBerezhkovskiy<-function(logP, Fup, Fut, Vnl, Vph, Vw,PVnl, PVph, PVw){
  Kp<-(10^logP*(Vnl+0.3*Vph)+0.7*Vph+Vw/Fut)/(10^logP*(PVnl+0.3*PVph)+0.7*PVph+PVw/Fup)
  return(Kp)
}











  
  
  




