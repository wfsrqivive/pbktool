solveModel <- function (pars, tout, state, dosing){#input required for the solver
  
  derivs <- function(t, state, parms) { # returns rate of change
    with(as.list(c(state, parms)), {
      
      Vad = BW*FVad	 
      Vbo = BW*FVbo        	 
      Vbr = BW*FVbr		
      Vgu = BW*FVgu          	 
      Vhe = BW*FVhe        	
      Vki = BW*FVki
      Vli = BW*FVli 
      Vlu = BW*FVlu
      Vmu = BW*FVmu
      Vsk = BW*FVsk 
      Vsp = BW*FVsp
      Vte = BW*FVte
      Vve = BW*FVve
      Var = BW*FVar
      Vpl = BW*FVpl 
      Vrb = BW*FVrb 
      Vre = BW*FVre
      Vst = BW*FVst
      Vi1 = BW*FVi1
      Vi2 = BW*FVi2
      Vi3 = BW*FVi3
      Vi4 = BW*FVi4
      Vi5 = BW*FVi5
      Vi6 = BW*FVi6
      Vi7 = BW*FVi7
      Vco = BW*FVco
      
      Qad = QC*FQad            
      Qbo = QC*FQbo             
      Qbr = QC*FQbr             
      Qgu = QC*FQgu            
      Qhe = QC*FQhe             
      Qki = QC*FQki             
      Qh = QC*FQh              
      Qha = QC*FQh - QC*FQgu - QC*FQsp 
      Qlu = QC*FQlu		
      Qmu = QC*FQmu    
      Qsk = QC*FQsk    
      Qsp = QC*FQsp    
      Qte = QC*FQte    
      Qre = QC*FQre		
      
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
      #Cplasmavenous_nM = Cplasmavenous/MW*1e6	 	
      CLmet = (CLint/fuhep)*SF*Vli*60/1000
      
      
      Cliverfree = Cliver*fup		
      Ckidneyfree = Ckidney*fup
      Cplasmavenousfree = Cplasmavenous*fup
      
      
      
      Cst.diss = Ast.diss/Vst
      Ci1.diss = Ai1.diss/Vi1
      Ci2.diss = Ai2.diss/Vi2
      Ci3.diss = Ai3.diss/Vi3
      Ci4.diss = Ai4.diss/Vi4
      Ci5.diss = Ai5.diss/Vi5
      Ci6.diss = Ai6.diss/Vi6
      Ci7.diss = Ai7.diss/Vi7
      CCo.diss = Aco.diss/Vco
      Kdst = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Cst.diss) #Si = intrinsict solubility
      Kdi1 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci1.diss) #units cm, mg/l
      Kdi2 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci2.diss)
      Kdi3 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci3.diss)
      Kdi4 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci4.diss)
      Kdi5 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci5.diss)
      Kdi6 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci6.diss)
      Kdi7 = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-Ci7.diss)
      KdCo = 3*1E-4*60/(0.0005*1E6*0.003)*(Si*1000-CCo.diss)
      
      
      #differential equations intestinal tract
      dAst.undiss = -Kdst*Ast.undiss- GER*Ast.undiss  #GER = gastric emtying rate
      dAst.diss =  Kdst*Ast.undiss - GER*Ast.diss #- Ka*Ast.diss
      
      dAi1.undiss =  - Kdi1*Ai1.undiss + GER*Ast.undiss - Kt*Ai1.undiss 
      dAi1.diss = Kdi1*Ai1.undiss + GER*Ast.diss - Kt*Ai1.diss - Ka*Ai1.diss
      
      dAi2.undiss = - Kdi2*Ai2.undiss + Kt*Ai1.undiss - Kt*Ai2.undiss 
      dAi2.diss = Kdi2*Ai2.undiss + Kt*Ai1.diss - Kt*Ai2.diss - Ka*Ai2.diss
      
      dAi3.undiss =  - Kdi3*Ai3.undiss + Kt*Ai2.undiss - Kt*Ai3.undiss 
      dAi3.diss = Kdi3*Ai3.undiss + Kt*Ai2.diss - Kt*Ai3.diss - Ka*Ai3.diss
      
      dAi4.undiss =  - Kdi4*Ai4.undiss + Kt*Ai3.undiss - Kt*Ai4.undiss 
      dAi4.diss = Kdi4*Ai4.undiss + Kt*Ai3.diss - Kt*Ai4.diss - Ka*Ai4.diss
      
      dAi5.undiss =  - Kdi5*Ai5.undiss + Kt*Ai4.undiss - Kt*Ai5.undiss 
      dAi5.diss = Kdi5*Ai5.undiss + Kt*Ai4.diss - Kt*Ai5.diss - Ka*Ai5.diss
      
      dAi6.undiss =  - Kdi6*Ai6.undiss + Kt*Ai5.undiss - Kt*Ai6.undiss 
      dAi6.diss = Kdi6*Ai6.undiss + Kt*Ai5.diss - Kt*Ai6.diss - Ka*Ai6.diss
      
      dAi7.undiss =  - Kdi7*Ai7.undiss + Kt*Ai6.undiss - Kt*Ai7.undiss 
      dAi7.diss = Kdi7*Ai7.undiss + Kt*Ai6.diss - Kt*Ai7.diss - Ka*Ai7.diss
      
      dAco.undiss =  - KdCo*Aco.undiss+ Kt*Ai7.undiss- Ka*Aco.undiss 
      dAco.diss = KdCo*Aco.undiss+Kt*Ai7.diss- Ka*Aco.diss
      
      sumUptake = Ka*Ai1.diss + Ka*Ai2.diss + 
        Ka*Ai3.diss+ Ka*Ai4.diss+ Ka*Ai5.diss+ 
        Ka*Ai6.diss + Ka*Ai7.diss #+Ka*Ast.diss
      
      
      
      
      #differential equations body
      dAad = Qad*(Carterial - Cadipose/Kpad*BP)										
      dAbo = Qbo*(Carterial - Cbone/Kpbo*BP)										
      dAbr = Qbr*(Carterial - Cbrain/Kpbr*BP)										
      dAgu = Ka*D + Qgu*(Carterial - Cgut/Kpgu*BP) #			
      dAhe = Qhe*(Carterial - Cheart/Kphe*BP)										
      dAki = Qki*(Carterial - Ckidney/Kpki*BP) - CLrenal*Ckidneyfree							
      dAli= sumUptake+ Qha*Carterial + Qgu*(Cgut/Kpgu*BP) + Qsp*(Cspleen/Kpsp*BP) - Qh*(Cliver/Kpli*BP) - Cliverfree*CLmet	 
      dAlu = Qlu*Cvenous - Qlu*(Clung/Kplu*BP)
      dAmu = Qmu*(Carterial - Cmuscle/Kpmu*BP)
      dAsk = Qsk*(Carterial - Cskin/Kpsk*BP)	
      dAsp = Qsp*(Carterial - Cspleen/Kpsp*BP)
      dAte = Qte*(Carterial - Ctestes/Kpte*BP)
      dAve = Qad*(Cadipose/Kpad*BP) + Qbo*(Cbone/Kpbo*BP) + Qbr*(Cbrain/Kpbr*BP) +  Qhe*(Cheart/Kphe*BP) + Qki*(Ckidney/Kpki*BP) + Qh*(Cliver/Kpli*BP)  + Qmu*(Cmuscle/Kpmu*BP) + Qsk*(Cskin/Kpsk*BP) + Qte*(Ctestes/Kpte*BP) + Qre*(Crest/Kpre*BP) - Qlu*Cvenous 						
      dAar = Qlu*(Clung/Kplu*BP) - Qlu*Carterial
      dAre = Qre*(Carterial - Crest/Kpre*BP)		
      dD = - Ka*D
      
      #{Defining amount metabolitezed and cleared by the kidney for the mass balance equation}
      dAliClearance =   Cliverfree*CLmet
      dAkiClearance =   CLrenal*Cplasmavenousfree
      
      #{Mass Balance}
      Total = dose*BW
      Calculated = Aad+Abo+Abr+Agu+Ahe+Aki+Ali+
        Alu+Amu+Ask+Asp+Ate+Ave+Aar+Are+D+
        AliClearance+AkiClearance+
        Ast.undiss+Ast.diss+Ai1.undiss+Ai1.diss+
        Ai2.undiss+Ai2.diss+Ai3.undiss+Ai3.diss+
        Ai4.undiss+Ai4.diss+Ai5.undiss+Ai5.diss+
        Ai6.undiss+Ai6.diss+Ai7.undiss+Ai7.diss+
        Aco.undiss+Aco.diss
      
      ERROR= ((Total-Calculated)/Total+1E-30)*100
      MASSBBAL=Total-Calculated + 1 
      
      
      
      res<-c(Ast.undiss, dAst.diss, 
             dAi1.undiss, dAi1.diss,
             dAi2.undiss, dAi2.diss,
             dAi3.undiss, dAi3.diss,
             dAi4.undiss, dAi4.diss,
             dAi5.undiss, dAi5.diss,
             dAi6.undiss, dAi6.diss,
             dAi7.undiss, dAi7.diss,
             dAco.undiss, dAco.diss,
             dAad, dAbo, dAbr, dAgu, 
             dAhe, dAki, dAli, dAlu, 
             dAmu, dAsk, dAsp, dAte, 
             dAve, dAar, dAre, dD, 
             dAliClearance, dAkiClearance) 
      
      return(list(res, ERROR = ERROR, 
                  MASSBBAL = MASSBBAL,
                  Calculated = Calculated,
                  Cplasmavenous = Cplasmavenous,#model of Jones and Rowland is in mg/L
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
                  sumUptake = sumUptake, 
                  Kdst = Kdst, 
                  Si = Si))
      
    }) #end with
  } #end derivs
  
  as.data.frame(ode(y = state, times = tout, parms = pars,
                    atol=1e-50,
                    events = list(data = dosing), 
                    func = derivs))#
  
}

