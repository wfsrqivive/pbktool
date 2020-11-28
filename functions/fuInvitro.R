##### fraction unbound in the in vitro hepatic clearance incubation according to Kilford et al. 2008

f.Fuinc<-function(enzymeSource, fracUnionized, fracIonizedAcid, fracIonizedBase, logD, logP, concHep){
    if (enzymeSource == "hepatocytes1"||enzymeSource == "hepatocytes2" & fracUnionized >= 0.95 || enzymeSource == 1 & fracIonizedAcid >= 0.95){
      Fuinc<-1/(1+125*(concHep/1*0.005)*10^(0.072*logD^2+0.067*logD-1.126))
  } else if(enzymeSource == "hepatocytes1" ||enzymeSource == "hepatocytes2" & fracIonizedBase >= 0.95 ){
    Fuinc<-1/(1+125*(concHep/1*0.005)*10^(0.072*logP^2+0.067*logP-1.126))
  } else if(enzymeSource == "hepatocytes1" || enzymeSource == "hepatocytes2" & fracUnionized < 0.95 || enzymeSource == 1 &fracIonizedAcid < 0.95 || enzymeSource == 1 & fracIonizedBase < 0.95){
    Fuinc<-(fracIonizedAcid+fracUnionized)*(1/(1+125*(concHep/1*0.005)*10^(0.072*logD^2+0.067*logD-1.126))) + fracIonizedBase*(1/(1+125*(concHep/1*0.005)*10^(0.072*logP^2+0.067*logP-1.126))) 
  } else if (enzymeSource == "microsomes1" || enzymeSource == "microsomes2"& fracUnionized >= 0.95 || enzymeSource == 1 & fracIonizedAcid >= 0.95){
    Fuinc<-1/(1+(concHep*10^(0.072*logD^2+0.067*logD-1.126)))
  } else if(enzymeSource == "microsomes1" ||enzymeSource == "microsomes2" & fracIonizedBase >= 0.95 ){
    Fuinc<-1/(1+(concHep*10^(0.072*logP^2+0.067*logP-1.126)))
  } else if(enzymeSource == "microsomes1" || enzymeSource == "microsomes2" & fracUnionized < 0.95 || enzymeSource == 1 &fracIonizedAcid < 0.95 || enzymeSource == 1 & fracIonizedBase < 0.95){
    Fuinc<-(fracIonizedAcid+fracUnionized)*(1/(1+(concHep*10^(0.072*logD^2+0.067*logD-1.126)))) + fracIonizedBase*(1/(1+(concHep*10^(0.072*logP^2+0.067*logP-1.126)))) 
  } else if (enzymeSource == "S91" || enzymeSource == "S92"& fracUnionized >= 0.95 || enzymeSource == 1 & fracIonizedAcid >= 0.95){
    Fuinc<-1/(1+(concHep*10^(0.072*logD^2+0.067*logD-1.126)))
  } else if(enzymeSource == "S91" ||enzymeSource == "S92" & fracIonizedBase >= 0.95 ){
    Fuinc<-1/(1+(concHep*10^(0.072*logP^2+0.067*logP-1.126)))
  } else if(enzymeSource == "S91" || enzymeSource == "S92" & fracUnionized < 0.95 || enzymeSource == 1 &fracIonizedAcid < 0.95 || enzymeSource == 1 & fracIonizedBase < 0.95){
    Fuinc<-(fracIonizedAcid+fracUnionized)*(1/(1+(concHep*10^(0.072*logD^2+0.067*logD-1.126)))) + fracIonizedBase*(1/(1+(concHep*10^(0.072*logP^2+0.067*logP-1.126)))) 
    
  } else {
    Fuinc<-200
  }
  return(Fuinc)
}

