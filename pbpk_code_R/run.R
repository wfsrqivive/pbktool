library(tidyverse)
library(deSolve)
source("model/model.R")


chemical_specific_parameters<-
  readxl::read_excel("input/chemical_specific_parameters.xlsx") %>%
  #the following  lines convert the data so that each row represents a scenario
  #this is needed to perform the PBK-model simulations per row/per scenario
  #Scenarios can be added in the chemical_specific_parameters.csv file, but do not change the format of first 4 columns of the .csv file. 
  gather(Scenario, value, 5:last_col(), -parameter) %>%
  select (parameter, Scenario, value) %>%
  spread(parameter, value)%>% 
  # conversion to numeric where needed
  mutate_at(vars(-one_of("Scenario", "compound", "species")), as.numeric)

physiological_parameters<-
  readxl::read_excel("input/physiological_parameters.xlsx", sheet ="physiological_parameters") %>%
  #the following three lines convert the data so that each row represents 
  #a species. This to combine the data with the chemical specific data as imported above 
  gather(species, value, rat:human, -parameter) %>%
  select (parameter, species, value) %>%
  spread(parameter, value)

input<-full_join(chemical_specific_parameters, physiological_parameters, by = "species")

runPBPKmodel <- input %>%
  group_by(Scenario, compound, MW, dose,species) %>%
  summarise(model(pars = c(dose = dose,                 
                  #physiological data as defined for rat and human in the physiological input data
                  QC = QC, 
                  Qad = QC*FQad, Qbo = QC*FQbo, Qbr = QC*FQbr, Qgu = QC*FQgu,  
                  Qhe = QC*FQhe, Qki = QC*FQki, Qh = QC*FQh, Qha = QC*FQh-QC*FQgu-QC*FQsp, Qlu = QC*FQlu,	         	
                  Qmu = QC*FQmu, Qsk = QC*FQsk, Qsp = QC*FQsp, Qte = QC*FQte, Qre = QC*FQre,
                  BW = BW,
                  Vad = BW*FVad,  Vbo = BW*FVbo, Vbr = BW*FVbr,	Vgu = BW*FVgu,	
                  Vhe = BW*FVhe, Vki = BW*FVki,	Vli = BW*FVli,	Vlu = BW*FVlu,		
                  Vmu = BW*FVmu, Vsk = BW*FVsk, Vsp = BW*FVsp, Vte = BW*FVte,                	 
                  Vve = BW*FVve, Var = BW*FVar, Vre = BW*FVre,          	 
                  GER = GER, Kt = Kt, Kco = Kco,
                  Vst = BW*FVst, Vi1 = BW*FVi1, Vi2 = BW*FVi2,
                  Vi3 = BW*FVi3, Vi4 = BW*FVi4, Vi5 = BW*FVi5, 
                  Vi6 = BW*FVi6, Vi7 = BW*FVi7,
                  Vco = BW*FVco,
                  #partition coefficients 
                  Kpad = Kpad, Kpbo = Kpbo, Kpbr = Kpbr, Kpgu = Kpgu, 
                  Kphe = Kphe,  Kpki = Kpki,  Kpli = Kpli,	Kplu = Kplu,	
                  Kpmu = Kpmu, Kpsk = Kpsk, Kpsp = Kpsp,  Kpte = Kpte, Kpre = Kpre,  	
                  #remaining chemical specific input parameters
                  BP = BP, 
                  fup= fup, 
                  fuhep = fuhep, 
                  CLint = CLint,
                  GFR = GFR,
                  SF= SF,
                  MW = MW, 
                  Ka = ka), 
                #defining the output time-frame, in this case every 0.1 hour, a result is predicted 
                tout  = seq(0, 24, by = 0.1),
                #what is the amount of chemical present in a compartment at t = 0
                state = c(Ast.diss = 0,
                          Ai1.diss = 0,
                          Ai2.diss = 0,
                          Ai3.diss = 0,
                          Ai4.diss = 0,
                          Ai5.diss = 0,
                          Ai6.diss = 0,
                          Ai7.diss = 0,
                          Aco.diss = 0,
                          Aad=0, Abo=0, Abr=0, Agu=0, 
                          Ahe=0, Aki=0, Ali=0, Alu=0, 
                          Amu=0, Ask=0, Asp=0, Ate=0, 
                          Ave=0, Aar=0, Are=0, 
                          AliClearance = 0, AkiClearance = 0), # starting amounts 
                #doseing frequency    
                dosing =data.frame(var = c("Ast.diss"),
                                   time = c(seq(0, 24, by = 24)), 
                                   value = c(dose*BW), 
                                   method = c("add"))
  )) %>%
  #postrun modifications
  mutate(Time.hr = time) %>%
  mutate(Time.min = time*60) %>%
  mutate(Cplasmavenous.ng.p.ml = Cplasmavenous*1000) %>%
  mutate(Cplasmavenous.nM = Cplasmavenous.ng.p.ml/MW*1000) %>%
  mutate(Scenario = paste0(compound, ";",species,";",dose, "mg/kg bw")) %>%
  #because there is a delay between the dosing and the calculating of the massbalance, the following line sets the ERROR and MASSBAL to the expected defaults
  mutate(ERROR = ifelse (time ==0, 0, ERROR), 
         MASSBAL = ifelse (time ==0, 1, MASSBAL))

#outputplot 
ggplot(data = runPBPKmodel)+
  geom_line(aes(x = Time.hr, y = Cplasmavenous.ng.p.ml, color = Scenario))+
  theme_light()

# Saving the plot and data
# ggsave (saving the plot): path corresponds the folder in which the figure is saved
# filename can be changed
ggsave(path = "results", filename = "pbpkPlot.png",width=8, height=4,dpi=300)
# write_xlsx
openxlsx::write.xlsx(runPBPKmodel, file = "results/pbpkOutput.xlsx")