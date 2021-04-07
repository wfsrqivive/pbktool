model <- function (pars, tout, state, dosing){#input required for the solver
  
  derivs <- function(t, state, parms) { #returns rate of change
    with(as.list(c(state, parms)), {
      #differential equations intestinal tract
      Cst.diss = Ast.diss/Vst #concentration dissolved stomach
      Ci1.diss = Ai1.diss/Vi1
      Ci2.diss = Ai2.diss/Vi2
      Ci3.diss = Ai3.diss/Vi3
      Ci4.diss = Ai4.diss/Vi4
      Ci5.diss = Ai5.diss/Vi5
      Ci6.diss = Ai6.diss/Vi6
      Ci7.diss = Ai7.diss/Vi7
      CCo.diss = Aco.diss/Vco
      
      dAst.diss = - GER*Ast.diss #amount dissolved stomach
      dAi1.diss = GER*Ast.diss - Kt*Ai1.diss - Ka*Ai1.diss
      dAi2.diss = Kt*Ai1.diss - Kt*Ai2.diss - Ka*Ai2.diss
      dAi3.diss = Kt*Ai2.diss - Kt*Ai3.diss - Ka*Ai3.diss
      dAi4.diss = Kt*Ai3.diss - Kt*Ai4.diss - Ka*Ai4.diss
      dAi5.diss = Kt*Ai4.diss - Kt*Ai5.diss - Ka*Ai5.diss
      dAi6.diss = Kt*Ai5.diss - Kt*Ai6.diss - Ka*Ai6.diss
      dAi7.diss = Kt*Ai6.diss - Kt*Ai7.diss - Ka*Ai7.diss
      dAco.diss = Kt*Ai7.diss
      
      sumUptake = Ka*Ai1.diss + Ka*Ai2.diss + 
        Ka*Ai3.diss+ Ka*Ai4.diss+ Ka*Ai5.diss+ 
        Ka*Ai6.diss + Ka*Ai7.diss 
      
      #differential equations body
      Cadipose = Aad/Vad		
      Cbone = Abo/Vbo		 
      Cbrain = Abr/Vbr	
      Cgut = Agu/Vgu		
      Cheart = Ahe/Vhe	 
      Ckidney = Aki/Vki		
      Cliver = Ali/Vli	 
      Clung = Alu/Vlu			 
      Cmuscle = Amu/Vmu		
      Cskin = Ask/Vsk		 
      Cspleen = Asp/Vsp		 
      Ctestes = Ate/Vte		 
      Cvenous = Ave/Vve		
      Carterial = Aar/Var	
      Crest = Are/Vre 		 
      Cplasmavenous = Cvenous/BP
      Cliverfree = Cliver*fup		
      Cplasmavenousfree = Cplasmavenous*fup
      CLmet = (CLint/fuhep)*SF*Vli*60/1000
      
      dAad = Qad*(Carterial - Cadipose/Kpad*BP)										
      dAbo = Qbo*(Carterial - Cbone/Kpbo*BP)										
      dAbr = Qbr*(Carterial - Cbrain/Kpbr*BP)										
      dAgu = Qgu*(Carterial - Cgut/Kpgu*BP) 	
      dAhe = Qhe*(Carterial - Cheart/Kphe*BP)										
      dAki = Qki*(Carterial - Ckidney/Kpki*BP) - GFR*Cplasmavenousfree							
      dAli= sumUptake + Qha*Carterial + Qgu*(Cgut/Kpgu*BP) + Qsp*(Cspleen/Kpsp*BP) - Qh*(Cliver/Kpli*BP)  - Cliverfree*CLmet	
      dAlu = Qlu*Cvenous - Qlu*(Clung/Kplu*BP)
      dAmu = Qmu*(Carterial - Cmuscle/Kpmu*BP)
      dAsk = Qsk*(Carterial - Cskin/Kpsk*BP)	
      dAsp = Qsp*(Carterial - Cspleen/Kpsp*BP)
      dAte = Qte*(Carterial - Ctestes/Kpte*BP)
      dAve = Qad*(Cadipose/Kpad*BP) + 
             Qbo*(Cbone/Kpbo*BP) + 
             Qbr*(Cbrain/Kpbr*BP) +  
             Qhe*(Cheart/Kphe*BP) + 
             Qki*(Ckidney/Kpki*BP) + 
             Qh*(Cliver/Kpli*BP)  + 
             Qmu*(Cmuscle/Kpmu*BP) + 
             Qsk*(Cskin/Kpsk*BP) + 
             Qte*(Ctestes/Kpte*BP) + 
             Qre*(Crest/Kpre*BP) - Qlu*Cvenous 						
      dAar = Qlu*(Clung/Kplu*BP) - Qlu*Carterial
      dAre = Qre*(Carterial - Crest/Kpre*BP)		
      
      #{Defining amount metabolitezed and cleared by the kidney for the mass balance equation}
      dAliClearance =   Cliverfree*CLmet
      dAkiClearance =   GFR*Cplasmavenousfree
      
      #{Mass Balance}
      Total = dose*BW
      Calculated = 
        Ast.diss+Ai1.diss+
        Ai2.diss+Ai3.diss+
        Ai4.diss+Ai5.diss+
        Ai6.diss+Ai7.diss+
        Aco.diss+
        Aad+Abo+Abr+Agu+Ahe+Aki+Ali+
        Alu+Amu+Ask+Asp+Ate+Ave+Aar+Are+
        AliClearance+AkiClearance
      
      ERROR= ((Total-Calculated)/Total+1E-30)*100
      MASSBAL=Total-Calculated + 1 
      
      res<-c(dAst.diss, 
             dAi1.diss,
             dAi2.diss,
             dAi3.diss,
             dAi4.diss,
             dAi5.diss,
             dAi6.diss,
             dAi7.diss,
             dAco.diss,
             dAad, dAbo, dAbr, dAgu, 
             dAhe, dAki, dAli, dAlu, 
             dAmu, dAsk, dAsp, dAte, 
             dAve, dAar, dAre, 
             dAliClearance, dAkiClearance) 
      
      return(list(res, ERROR = ERROR, 
                  MASSBAL = MASSBAL,
                  #concentrations in mg/L
                  Calculated = Calculated,
                  Cplasmavenous = Cplasmavenous,
                  Cadipose = Cadipose, 
                  Cbone = Cbone,	 
                  Cbrain = Cbrain,	
                  Cgut = Cgut,		
                  Cheart = Cheart,	 
                  Ckidney = Ckidney,		
                  Cliver = Cliver,	 
                  Clung = Clung,			 
                  Cmuscle = Cmuscle,		
                  Cskin = Cskin,	 
                  Cspleen = Cspleen,		 
                  Cgonads = Ctestes,		 
                  Cvenous = Cvenous,		
                  Carterial = Carterial, 
                  Cst.diss = Cst.diss, 
                  Ci1.diss = Ci1.diss, 
                  Ci2.diss = Ci2.diss, 
                  Ci3.diss = Ci3.diss, 
                  Ci4.diss = Ci4.diss, 
                  Ci5.diss = Ci5.diss,
                  Ci6.diss = Ci6.diss, 
                  Ci7.diss = Ci7.diss, 
                  CCo.diss = CCo.diss))
      
    }) #end with
  } #end derivs
  
  as.data.frame(ode(y = state, times = tout, parms = pars,
                    events = list(data = dosing), 
                    func = derivs))
  
}

