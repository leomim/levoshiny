#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

#All libraries needed for the application
library(shiny)
library(ggplot2) 
library(saemix)
library(shinythemes)
library(shinydashboard)
library(bslib)
library(DT)

############FONCTION SAEMIXPREDICTNEWDATA######################
#SaemixPredictNewdata is a function needed to make bayesian prediction and will be implemented in the new version of the saemix package
saemixPredictNewdata<-function(saemixObject, newdata, type=c("ipred", "ypred", "ppred", "icpred"),nsamp=1) {
  set.seed(12345)
  saemixObject<-replaceData.saemixObject(saemixObject,newdata)
  if(sum(is.na(saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]))>0) {
    if(saemixObject["model"]["modeltype"]=="likelihood") {
      if(saemixObject["options"]$warnings) cat("Please provide values of the response to obtain predictions for a model defined by loglikelihood\n")
      return(NULL)
    }
    type<-type[type!="ipred" & type!="icpred"]
  }
  if(length(type)==0) type<-"ppred"
  
  # Estimate population parameters (Ci*mu) for the new subjects
  saemixObject<-estimateMeanParametersNewdata(saemixObject) # updates mean.phi (normally...)
  
  # Predictions using the mean parameters ppred=f(E(theta,x))
  newdata<-saemixObject["data"]
  chdat<-saemixObject["rep.data"]
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,c(newdata["name.predictors"],newdata["name.cens"],newdata["name.mdv"],newdata["name.ytype"]),drop=FALSE]
  mean.phi<-saemixObject["results"]["mean.phi"]
  psiM<-transphi(mean.phi,saemixObject["model"]["transform.par"])
  fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
  colnames(psiM)<-saemixObject["model"]["name.modpar"]
  predictions<-data.frame(IdM,XM,ppred=fpred)
  colnames(predictions)[1]<-newdata["name.group"]
  parameters<-list(id=unique(newdata["data"][,newdata["name.group"]]), population=psiM)
  saemixObject["results"]["ppred"]<-fpred
  
  # Mean predictions over the population ypred=E(f(theta,x))
  # Technically... we can obtain ypred as E(f()) by simulating etas in their distribution
  if(length(grep(c("ypred"),type))>0) {
    ind.eta<-saemixObject["model"]["indx.omega"]
    nb.etas<-length(ind.eta)
    omega<-saemixObject["results"]["omega"]
    chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
    
    etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)
    ypred<-matrix(data=0,nrow=dim(XM)[1],ncol=saemixObject["options"]$nb.sim)
    mean.phiM<-do.call(rbind,rep(list(mean.phi),1))
    phiMc<-mean.phiM
    if(length(grep(c("ypred"),type))==1) {
      for(isim in 1:saemixObject["options"]$nb.sim) {
        etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
        phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
        psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
        ypred[,isim]<-fpred
      }
      ypred<-rowMeans(ypred)
      predictions$ypred<-ypred
      saemixObject["results"]["ypred"]<-ypred
    }
  }
  if(sum(!is.na(match(c("icpred","ipred"),type)))==0) { # only population predictions
    return(list(param=parameters,predictions=predictions, object=saemixObject))
  }
  # Estimate individual parameters, if type contains ipred and/or icpred
  
  ctype<-c()
  if(length(grep("ipred",type))==1) ctype<-c(ctype,"mode")
  if(length(grep("icpred",type))==1) ctype<-c(ctype,"mean")
  saemixObject<-estimateIndividualParametersNewdata(saemixObject,type=ctype,nsamp=nsamp) # updates cond.mean.psi, map.psi (normally...)
  
  if(length(grep("icpred",type))==1) {
    psiM<-parameters$cond.mean.psi<-saemixObject["results"]["cond.mean.psi"]
    parameters$cond.var.phi<-saemixObject["results"]["cond.var.phi"]
    parameters$cond.mean.phi<-saemixObject["results"]["cond.mean.phi"]
    fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
    predictions<-cbind(predictions,icpred=fpred)
    saemixObject["results"]["icpred"]<-fpred
    if(nsamp>1) {
      samp.pred<-array(dim=c(length(fpred),nsamp))
      samp.par<-array(dim=c(dim(psiM),nsamp))
      for(isamp in 1:nsamp) {
        phiM<-saemixObject["results"]["phi.samp"][,,isamp]
        if(is.null(dim(phiM))) phiM<-matrix(phiM, nrow=1)
        psiM<-samp.par[,,isamp]<-transphi(phiM,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
        samp.pred[,isamp]<-fpred
      }
    }
  }
  if(length(grep("ipred",type))==1) {
    psiM<-parameters$map.psi<-saemixObject["results"]["map.psi"]
    fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
    predictions<-cbind(predictions,ipred=fpred)
    saemixObject["results"]["ipred"]<-fpred
  }
  saemixObject["results"]["predictions"]<-predictions
  rlist<-list(param=parameters,predictions=predictions, object=saemixObject)
  if(length(grep("icpred",type))==1 & nsamp>1) 
    rlist<-list(param=parameters,predictions=predictions, object=saemixObject, parSample=samp.par, predSample=samp.pred)
  
  return(rlist)
}

###################END FUNCTION############

#User interface
ui <- dashboardPage(
  dashboardHeader(title="Levoshiny"),
  dashboardSidebar(
    sidebarMenu(id="tabs",
      menuItem("Description",tabName="Description",icon=icon("home")),
      menuItem("Inputs", tabName = "Inputs",icon=icon("info-circle")),
      menuItem("Plot", tabName = "Plot",icon=icon("chart-line")),
      menuItem("Predictions",tabName = "Predictions",icon=icon("award") ),
      menuItem('Conditional distributions',tabName="Conditional_distributions",icon=icon("warning")),
      menuItem('Credits',tabName='Credits',icon=icon("book"))
      )
  ),
  dashboardBody(
    
    #Sets the size of the ValueBoxOutput AUC
    tags$style(HTML("
        #auc {
          width: 100% !important; 
          height: 200px; 
        }

        #auc h3 {
          text-align: center; 
          line-height: 100px;
          font-size: 40px; 
        }

        #auc p {
          text-align: center; 
          font-size: 18px;
          
        #auc i {
          font-size: 18px !important;
        }
      ")),
    tabItems(
      #Sidebar Description (description of the app)
      tabItem(tabName="Description",
              title="Description", 
                  h3(tags$strong("Bayesian estimation of oral ofloxacin and levofloxacin pharmacokinetic parameters in osteoarticular infections")),
                  h4("Fill in the following fields in the 'Inputs' panel : daily dose, number of concentrations available, time of sampling and value of the concentrations obtained, gender, age and creatinine levels. You will obtain the predicted individual pharmacokinetic parameters and the AUC24h, using a modified version of the model of Lemaitre et al (1)."),
                  h4("The 'Data' table in the 'Inputs' panel will automatically upload the informations entered in the adjacent 'Dosing' box. The 'Plot', 'Predictions' and 'Conditional distributions' panels will be displayed/modified only after clicking on the 'Run model' button in the 'Inputs' panel."),
                  h6("1. Lemaitre F, Fily F, Foulquier JB, Revest M, Jullien V, Petitcollin A, et al. Development of a dosing-adjustment tool for fluoroquinolones in osteoarticular infections: The Fluo-pop study. Biomed Pharmacother. oct 2021;142:112053. "),
              box(title="WARNING",
                  h5(tags$strong('TO USE THIS APP, THE PATIENT MUST BE AT STEADY-STATE (minimum 72 hours after the introduction of levofloxacine).'),style = "color: red;"),
                  h5(tags$strong("GFR was calculated using the CKD-EPI formula as implemented in the Lemaitre et al model. The model can only be used for GFR values between 38.9 and 129.5ml/min/1.73m2.",style = "color: red;")),width="99%",status="danger",solidHeader = TRUE)),
   
      #Sidebar Inputs (main inputs with different types of inputs (numeric, select, slider) )
      tabItem(tabName = "Inputs",
              fluidRow(
                #box Dosage (inputs)
                box(title = "Dosage",
                    numericInput('Posology','Total daily dose (mg/24h)',1000,min = 0,max=4000,step=125,width=NULL),
                    selectInput('Admin',"Number of administrations per day", choices=c(1,2,3,4)),
                    selectInput('nb_dose', 'Number of concentrations available', choices = c(1,2,3,4)),
                    numericInput('Time1',"Time since last administration for the 1st concentration (hours)",12,min = 0,max=24,step=1,width=NULL),
                    numericInput('Conc1',
                   'First concentration value (mg/L)',
                   min= 0,
                   max= 20,
                   step= 0.1,
                   value=10),
                    conditionalPanel(
                    condition ="input.nb_dose>=2",
                    numericInput('Time2',"Time since last administration for the 2nd concentration (hours)",12,min = 0,max=24,step=1,width=NULL),
                    numericInput('Conc2',
                     'Second concentration value (mg/L)',
                     min= 0,
                     max= 20,
                     step= 0.1,
                     value=10)),
                      conditionalPanel(
                      condition ="input.nb_dose>=3",
                      numericInput('Time3',"Time since last administration for the 3rd concentration (hours)",12,min = 0,max=24,step=1,width=NULL),
                      numericInput('Conc3',
                     'Third concentration value (mg/L)',
                     min= 0,
                     max= 20,
                     step= 0.1,
                     value=10)),
                      conditionalPanel(
                      condition ="input.nb_dose==4",
                      numericInput('Time4',"Time since last administration for the 4th concentration (hours)",12,min = 0,max=24,step=1,width=NULL),
                      numericInput('Conc4',
                     'Fourth concentration value (mg/L)',
                     min= 0,
                     max= 20,
                     step= 0.1,
                     value=10)),
                   selectInput("Sex",
                               "Gender",choices=c("Male","Female")),
                   sliderInput("Age",
                               "Age (years)",
                               min = 18,
                               max = 120,
                               value = 50),
                   sliderInput("Creatinine",
                               "Creatinine (µmol/L)",
                               min = 0,
                               max = 600,
                               value = 80),
                    actionButton('run','Run model',style="color: #fff; background-color: #e95420; border-color: #c34113;border-radius: 10px; border-width: 2px"),solidHeader = TRUE,status="primary"),
              #Box Data (real-time data)
              box(title="Data", 
                  DTOutput("Tab"),solidHeader = TRUE,status="success")
    
  )),
          #Sidebar Plot
      tabItem(tabName = "Plot",
              conditionalPanel(
               condition= "input.run>0",
              box(title="Predicted profile",
              plotOutput('Predind'),uiOutput('textplot'),width = "100%", height = "100%",solidHeader = TRUE,status="info")
              )),
          #Sidebar Predictions
      tabItem(tabName="Predictions",
              conditionalPanel(
                condition="input.run>0",
              fluidRow(
              #box AUC
              valueBoxOutput("auc")),
              #box Parameters
              box(title="Parameters",
                  textOutput('title'),tableOutput('ip'),uiOutput('textparameters'),solidHeader = TRUE,status="info",width="100%"),
              #box Prediction errors
              box(title = "Predictions errors",
                  tableOutput('tabconc'),uiOutput('textpe'),solidHeader = TRUE,status="info",width="100%"
                  )
              )
              ),
  
          #Sidebar Conditional_distributions
      tabItem(tabName = "Conditional_distributions",
              conditionalPanel(
                condition="input.run>0",
              box(title="Conditional distributions of simulated subjects and of the concerned patients",
              plotOutput('conddist'), uiOutput('textconddist'),width="100%",height="100%",solidHeader = TRUE,status="info")
              )),
  
          #Sidebar Credits
      tabItem(tabName="Credits",
              h5("Léo MIMRAM, IAME Inserm UMR1137, Bichat Hospital, France, leo.mimram@aphp.fr"),
              h5("Emmanuelle COMETS, IAME Inserm UMR1137, Bichat Hospital, France, emmanuelle.comets@inserm.fr"),
              h5("Julie BERTRAND, IAME Inserm UMR1137, Bichat Hospital, France, julie.bertrand@inserm.fr"),
              h5("Sophie MAGREAULT, Pharmacology Department, Jean Verdier Hospital, Bondy, France and IAME Inserm UMR1137, Bichat Hospital, France, sophie.magreault@aphp.fr"),
              h5("Vincent JULLIEN, Pharmacology Department, Jean Verdier hospital, Bondy, France and IAME Inserm UMR1137, Bichat Hospital, France, vincent.jullien@aphp.fr"))
  )
)
)


# Server function
server <- function(input, output,session) {
  data <- reactive({
    if (input$Sex == "Male") {
      if (input$Creatinine <= 80) {
        GFR <- 141 * ((input$Creatinine/88.4 / 0.9) ^ -0.411) * (0.993 ^ input$Age)
      } else {
        GFR <- 141 * ((input$Creatinine/88.4 / 0.9) ^ -1.209) * (0.993 ^ input$Age)
      }
    } else {
      if (input$Creatinine <= 62) {
        GFR <- 144 * ((input$Creatinine/88.4 / 0.7) ^ -0.329) * (0.993 ^ input$Age)
      } else {
        GFR <- 144 * ((input$Creatinine/88.4 / 0.7) ^ -1.209) * (0.993 ^ input$Age)
      }
    } # calculates CKD-EPI with Creatinine (µmol/L) input
    
    posology <- as.numeric(input$Posology)
    admin <- as.numeric(input$Admin)
    time1 <- as.numeric(input$Time1)
    conc1 <- as.numeric(input$Conc1)
    Dfg<-GFR
    Gender<-ifelse(input$Sex=="Male",1,0)
    times <- c(time1)
    concentrations <- c(conc1)
    
    if (input$nb_dose >= 2) {
      time2 <- as.numeric(input$Time2)
      conc2 <- as.numeric(input$Conc2)
      times <- c(times, time2)
      concentrations <- c(concentrations, conc2)
    }
    if (input$nb_dose >= 3) {
      time3 <- as.numeric(input$Time3)
      conc3 <- as.numeric(input$Conc3)
      times <- c(times, time3)
      concentrations <- c(concentrations, conc3)
    }
    if (input$nb_dose == 4) {
      time4 <- as.numeric(input$Time4)
      conc4 <- as.numeric(input$Conc4)
      times <- c(times, time4)
      concentrations <- c(concentrations, conc4)
    }
    
    tau_value <- 24 / admin
    
    data.frame(
      Id = rep(1, length(times)),
      Time = times,
      Dose = rep((posology/admin), length(times)),
      tau = rep(tau_value, length(times)),
      Concentration = concentrations,
      GFR=round(Dfg,2),
      Sex=Gender,
      stringsAsFactors = FALSE  
    )
    
  })
  
  ############### MODELE SAEMIX###############
  
  #3 random patients to generate the saemix.model
  #the objective here is simply to create a saemix object with random patients and data
  #the only way to create a saemix object is first to invent patients in a saemix.data (dataframe), used with an saemix.model to create an object (here : lemaitre.fit)
  #The values of the parameters will then be replaced by the values of the chosen model
  
  ID<-c(1,1,2,2,3,3,4,4,5,5)
  DOSE<-c(500,500,200,200,750,750,750,750,200,200)
  TIME<-c(3,12,3,12,3,12,3,24,3,8)
  CONCENTRATION<-c(9,2.5,3,1,12,4,10,3,4,2)
  Creatinine<-c(80,80,90,90,50,50,60,60,100,100)
  Sex<-c(0,0,0,0,1,1,1,1,0,0)
  Age<-c(22,22,66,66,88,88,50,50,40,40)
  tau<-c(12,12,12,12,12,12,24,24,8,8)
  
  #Calcul DFG CKD-EPI
  DFG<-ifelse(Sex == 0,
              ifelse(Creatinine <= 80,
                     141 * ((Creatinine/88.4/ 0.9) ^ -0.411) * (0.993 ^ Age),
                     141 * ((Creatinine/88.4 / 0.9) ^ -1.209) * (0.993 ^ Age)),
              ifelse (Creatinine <= 62,
                      144 * ((Creatinine[1]/88.4/ 0.7) ^ -0.329) * (0.993 ^ Age),
                      144 * ((Creatinine[1]/88.4/ 0.7) ^ -1.209) * (0.993 ^ Age)))
  
  
  essai<-data.frame(Id=ID,Dose=DOSE,Time=TIME,Concentration=CONCENTRATION,DFG=DFG,Sex=Sex,tau=tau)
  essai2<-cbind(essai,lDFG=log(essai$DFG/100))
  
  #saemix data
  saemix.data<-saemixData(name.data=essai2,header=TRUE,sep=" ",na=NA,
                          name.group=c("Id"),name.predictors=c("Time","Dose","tau"),
                          name.response=c("Concentration"),name.covariates=c("Sex","lDFG"),
                          units=list(x="hr",y="mg/L",covariates=c("-","ml/min/1.73m2")), name.X="Time")
  
  plot(saemix.data)
  
  #model used for simulations (here a 1-compartment model with a T-lag) 
  #analytical expressions are needed for saemix to work
  
  modelSS<-function(psi,id,xidep) {
    tim<-xidep[,1]   # Time AFTER dose
    dose<-xidep[,2]
    tau<-xidep[,3]
    ka<-psi[id,1]
    tlag<-psi[id,2]
    V<-psi[id,3]
    CL<-psi[id,4]
    V2<-psi[id,5]
    Q<-psi[id,6]
    k<-CL/V
    k12<-Q/V
    k21<-Q/V2
    beta = 1/2*(k12+k21+k-sqrt((k12+k21+k)**2-4*k21*k))
    alpha = (k21*k)/beta
    A=(ka/V)*((k21-alpha)/((ka-alpha)*(beta-alpha)))
    B=(ka/V)*((k21-beta)/((ka-beta)*(alpha-beta)))
    ypred<-dose*(((A*exp(-alpha*(tim-tlag+tau)))/(1-exp(-alpha*tau)))+((B*exp(-beta*(tim-tlag+tau)))/(1-exp(-beta*tau)))-(((A+B)*exp(-ka*(tim-tlag+tau)))/(1-exp(-ka*tau))))
    ypred2<-dose*(((A*exp(-alpha*(tim-tlag)))/(1-exp(-alpha*tau)))+((B*exp(-beta*(tim-tlag)))/(1-exp(-beta*tau)))-(((A+B)*exp(-ka*(tim-tlag)))/(1-exp(-ka*tau))))
    ypred[tim>=tlag]<-ypred2[tim>=tlag]
    return(ypred)
  }
  saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE,nbiter.saemix=c(10,10), fim=FALSE, map=FALSE)
  
  #Fixed parameters from the chosen model (Lemaitre modified here)
  psi<-data.frame(ka=1.92,tlag=0.42,V=30.4,CL=6.67,V2=42,Q=74)
  psi.lemaitre<-unlist(psi,use.names=FALSE)
  omega.lemaitre<-diag(c(0.6**2,0.94**2,0.75**2,0.25**2,0,0))
  sig.lemaitre<-c(0.2,0.02)
  covmat<-diag(6)
  covmat[5,5]<-covmat[6,6]<-0
  covmodel<-psi.cov<-matrix(data=0, nrow=2, ncol=6)
  covmodel[1,4]<-1
  covmodel[2,4]<-1
  beta.cov<-c(0.4,1.23)
  psi.cov[covmodel==1]<-beta.cov
  psi0<-rbind(psi.lemaitre,psi.cov)
  dimnames(psi0)<-list(NULL, c("ka","tlag","V","CL","V2","Q"))
  
  saemix.model<-saemixModel(model=modelSS,modeltype="structural",
                            description="Two-compartment model with first-order absorption",
                            psi0=psi0,
                            covariate.model = covmodel,
                            transform.par=c(1,1,1,1,1,1), covariance.model =covmat, 
                            omega.init=omega.lemaitre,
                            error.model="combined", error.init = sig.lemaitre) # 
  
  
  lemaitre.fit<-saemix(saemix.model,saemix.data,saemix.options)
  
  saemix.fit<-lemaitre.fit
  psi0.pop<-c(psi.lemaitre[1:4],beta.cov[1:2],psi.lemaitre[5:6])
  saemix.fit@results@fixed.psi<-c(psi.lemaitre[1:6])
  saemix.fit@results@fixed.effects<-psi0.pop
  saemix.fit@results@omega<-omega.lemaitre
  saemix.fit@results@respar<-sig.lemaitre
  saemix.fit@results@betas<-t(c(log(psi.lemaitre[1:4]),beta.cov[1:2],log(psi.lemaitre[5:6])))
  xsim1 <- estimateMeanParametersNewdata(saemix.fit)
  saemix.fit@results@mean.phi <- xsim1@results@mean.phi
  saemix.fit@options$warnings<-TRUE
  
  output$Tab <- renderDT({
    datatable(data(),
              options = list(
                scrollX = TRUE,  
                paging = FALSE, 
                searching = FALSE 
              ))
  })
  
  observeEvent(input$run, { #when the "runbutton" is activated
    
  updateTabItems(session, "tabs", "Plot")
    
  #implementation of the patient's data entered in the application  
  newdata<-data()
  newdata<-cbind(newdata,lDFG=log(newdata$GFR/100))  
  xpred<-saemixPredictNewdata(saemix.fit, newdata, type=c("ipred", "ypred", "ppred", "icpred"),nsamp=200)
  
  #Computation of the predicted AUC
  xpred$parSample[,,c(1:200)]->pr
  AUC200<-(as.numeric(input$Posology)/pr[4,])
  
  #prediction_interval AUC
  lower_bound <- round(quantile(AUC200, 0.025),2)
  upper_bound <- round(quantile(AUC200, 0.975),2)
  
  #prediction_interval KA
  pr[1,]->KA200
  lKA<-quantile(KA200, 0.025)
  uKA<-quantile(KA200,0.975)
  
  #prediction_interval Tlag
  pr[2,]->Tlag200
  lTlag<-quantile(Tlag200, 0.025)
  uTlag<-quantile(Tlag200,0.975)
  
  #prediction_interval V
  pr[3,]->V200
  lV<-quantile(V200, 0.025)
  uV<-quantile(V200,0.975)
  
  #prediction_interval CL
  pr[4,]->CL200
  lCL<-quantile(CL200, 0.025)
  uCL<-quantile(CL200,0.975)
  
  #prediction_interval V2
  pr[5,]->V2200
  lV2<-quantile(V2200, 0.025)
  uV2<-quantile(V2200,0.975)
  
  #prediction_interval Q
  pr[6,]->Q200
  lQ<-quantile(Q200, 0.025)
  uQ<-quantile(Q200,0.975)
  
  AUC<-(as.numeric(input$Posology)/mean(CL200))
  tab<-xpred$predictions[order(as.numeric(rownames(xpred$predictions))), ]
  tabfin<-cbind(tab,newdata$Concentration)
  RPE<-(tabfin$`newdata$Concentration`-tabfin$ipred)/tabfin$ipred
  TAB<-data.frame(Time=tabfin$Time,"Observed_concentration"=tabfin$`newdata$Concentration`,"Predicted_concentration"=tabfin$ipred,RPE)
  colnames(TAB)[colnames(TAB) == "Observed_concentration"] <- "Observed concentration"
  colnames(TAB)[colnames(TAB) == "Predicted_concentration"] <- "Predicted concentration"
  colnames(TAB)[colnames(TAB) == "RPE"] <- "Relative Precision Error"
  
  #Dataframe of the indivual pharmacokinetic parameters
  tab1<-data.frame(KAi=paste(round(mean(KA200),2),"[",round(lKA,2),";",round(uKA,2),"]"),
                   Tlagi=paste(round(mean(Tlag200),2),"[",round(lTlag,2),";",round(uTlag,2),"]"),
                   CLi=paste(round(mean(CL200),2),"[",round(lCL,2),";",round(uCL,2),"]"),
                   Vi=paste(round(mean(V200),2),"[",round(lV,2),";",round(uV,2),"]"),
                   V2i=paste(round(mean(V2200),2),"[",round(lV2,2),";",round(uV2,2),"]"),
                   Qi=paste(round(mean(Q200),2),"[",round(lQ,2),";",round(uQ,2),"]"))
  colnames(tab1)[colnames(tab1) == "KAi"] <- "KAi (/h)"
  colnames(tab1)[colnames(tab1) == "Tlagi"] <- "Tlagi (h)"
  colnames(tab1)[colnames(tab1) == "CLi"] <- "CLi (L/h)"
  colnames(tab1)[colnames(tab1) == "Vi"] <- "Vi (L)"
  colnames(tab1)[colnames(tab1) == "V2i"] <- "V2i (L)"
  colnames(tab1)[colnames(tab1) == "Qi"] <- "Qi (L/h)"
  
  tpar<-xpred$param
  xtim<-seq(0,24,0.5)
  xidep1<-data.frame(Id=rep(1,length(xtim)),Time=xtim,AMT=newdata$Dose[1],II=newdata$tau[1])
  psi1<-NULL
  id1<-xidep1[,1]
  meanparam<-c(mean(KA200),mean(Tlag200),mean(V200),mean(CL200),mean(V2200),mean(Q200))
  meanparam<-matrix(meanparam,nrow=1,ncol=6,byrow=TRUE)
  
  sm.icpred<-modelSS(tpar$cond.mean.psi,id1,xidep1[,-c(1)])
  sm.ipred<-modelSS(tpar$map.psi,id1,xidep1[,-c(1)])
  sm.ppred<-modelSS(tpar$population,id1,xidep1[,-c(1)])
  ypl<-data.frame(Id=xidep1$Id,Time=xidep1$Time,icpred=sm.icpred,ipred=sm.ipred, ppred=sm.ppred)
  
  # Predicting individual profile with variability using the envelope of the predictions
  parsamp<-xpred$parSample
  xp<-NULL
  for(i in 1:dim(parsamp)[3]) {
    psi1 <- parsamp[,,i]
    xp<-cbind(xp,modelSS(matrix(psi1, nrow=1),id1,xidep1[,-c(1)]))
  }
  predint<-apply(xp,1,quantile,c(0.025,0.5,0.975))
  
  ypl$y025<-predint[1,]
  ypl$y50<-predint[2,]
  ypl$y975<-predint[3,]
  
  #Plot
  if (input$Admin == 1) {
    obj<-ggplot(data=ypl,aes(x=Time, y=y50)) + scale_x_continuous("Time (hr)",limits = c(0, 24)) + scale_y_continuous("Predicted levofloxacine concentration (mg/L)") + geom_line(aes(x=Time, y=y50), linetype=1,colour="blue") + 
      geom_line(aes(x=Time, y=ppred),linetype=3, colour="blue")  + geom_ribbon(data=ypl,aes(x=Time, ymin=y025, ymax=y975),fill="lightblue",alpha=0.3)  + geom_point(data=newdata,aes(x=Time, y=Concentration), colour="blue")  +  scale_linetype_manual(name = "Predictions", values = c(1,3), labels=c("Individual","Population")) + theme(plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept=c(3,24),color="red",linetype="dashed")+ 
      theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
            axis.title=element_text(size=11),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9),
            title =element_text(size=12, face='bold'),
            strip.background = element_blank(),
            axis.line=element_line(colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "top")
  }
  
  if (input$Admin == 2){
    obj<-ggplot(data=ypl,aes(x=Time, y=y50)) + scale_x_continuous("Time (hr)",limits=c(0,12)) + scale_y_continuous("Predicted levofloxacine concentration (mg/L)") + geom_line(aes(x=Time, y=y50), linetype=1,colour="blue") + 
      geom_line(aes(x=Time, y=ppred),linetype=3, colour="blue")  + geom_ribbon(data=ypl,aes(x=Time, ymin=y025, ymax=y975),fill="lightblue",alpha=0.3)  + geom_point(data=newdata,aes(x=Time, y=Concentration), colour="blue")  +  scale_linetype_manual(name = "Predictions", values = c(1,3), labels=c("Individual","Population")) + theme(plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept=c(3,12),color="red",linetype="dashed")+ 
      theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
            axis.title=element_text(size=11),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9),
            title =element_text(size=12, face='bold'),
            strip.background = element_blank(),
            axis.line=element_line(colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "top")
  }
  
  if (input$Admin == 3){
    obj<-ggplot(data=ypl,aes(x=Time, y=y50)) + scale_x_continuous("Time (hr)",limits=c(0,8)) + scale_y_continuous("Predicted levofloxacine concentration (mg/L)") + geom_line(aes(x=Time, y=y50), linetype=1,colour="blue") + 
      geom_line(aes(x=Time, y=ppred),linetype=3, colour="blue")  + geom_ribbon(data=ypl,aes(x=Time, ymin=y025, ymax=y975),fill="lightblue",alpha=0.3)  + geom_point(data=newdata,aes(x=Time, y=Concentration), colour="blue")  +  scale_linetype_manual(name = "Predictions", values = c(1,3), labels=c("Individual","Population")) + theme(plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept=c(3,8),color="red",linetype="dashed")+
      theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
            axis.title=element_text(size=11),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9),
            title =element_text(size=12, face='bold'),
            strip.background = element_blank(),
            axis.line=element_line(colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "top")
  }
  
  if (input$Admin == 4){
    obj<-ggplot(data=ypl,aes(x=Time, y=y50)) + scale_x_continuous("Time (hr)",limits=c(0,6)) + scale_y_continuous("Predicted levofloxacine concentration (mg/L)") + geom_line(aes(x=Time, y=y50), linetype=1,colour="blue") + 
      geom_line(aes(x=Time, y=ppred),linetype=3, colour="blue")  + geom_ribbon(data=ypl,aes(x=Time, ymin=y025, ymax=y975),fill="lightblue",alpha=0.3)  + geom_point(data=newdata,aes(x=Time, y=Concentration), colour="blue") +  scale_linetype_manual(name = "Predictions", values = c(1,3), labels=c("Individual","Population")) + theme(plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept=c(3,6),color="red",linetype="dashed")+ 
      theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
            axis.title=element_text(size=11),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9),
            title =element_text(size=12, face='bold'),
            strip.background = element_blank(),
            axis.line=element_line(colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "top")
  }
  
  #simulation of 1000 x 5 (5 initial subjects)
  saemix.fit <- simulate(saemix.fit, nsim=1000)
  saemix.fit@options$displayProgress<-TRUE
  saemix.fit@options$warnings<-TRUE
  
  #conditional distributions of the simulated patients
  dcond2 <- conddist.saemix(saemix.fit, nsamp=200)
  par1 <- xpred$parSample
  par2 <- dcond2@results@psi.samp
  
  ypl<-NULL
  ypl5<-NULL
  for(ipar in 1:4) {
    ypl0<-data.frame(id=1:dcond2@data@N,ysamp=c(par2[,ipar,]),param=saemix.fit@model@name.modpar[ipar])
    ypl5<-rbind(ypl5, ypl0)
  }
  ypl5<-cbind(ypl5,method="Model PK")
  ypl<-rbind(ypl, ypl5)
  ypl5<-NULL
  for(ipar in 1:4) {
    ypl0<-data.frame(id=(1+dcond2@data@N),ysamp=c(par1[,ipar,]),param=saemix.fit@model@name.modpar[ipar])
    ypl5<-rbind(ypl5, ypl0)
  }
  ypl5<-cbind(ypl5,method="Patient PK")
  ypl<-rbind(ypl, ypl5)
  
  plot_distribcond <- ggplot(ypl,aes(ysamp,fill=as.factor(method),colour=as.factor(method)))+geom_density(alpha=0.4) +  facet_wrap(vars(param), ncol=2, scales="free")  + theme(legend.position = "top") + 
    scale_fill_discrete(name="Conditional distributions") + guides(colour="none")+ 
    labs(x = "values of the concerned parameters")+
    theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
          axis.title=element_text(size=11),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9),
          title =element_text(size=12, face='bold'),
          strip.background = element_blank(),
          axis.line=element_line(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")
  
  plot_distribcond2 <- ggplot(ypl,aes(ysamp,fill=as.factor(id),colour=as.factor(id)))+
    geom_density(alpha=0.4) + facet_wrap(vars(param), ncol=2, scales="free")  + 
    scale_fill_discrete(name="",labels=c("Model subject 1","Model subject 2","Model subject 3","Model subject 4","Modele subject 5","Patient PK")) + 
    guides(colour="none")+labs(x = "values of the concerned parameters")
    theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
          axis.title=element_text(size=11),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9),
          title =element_text(size=12, face='bold'),
          strip.background = element_blank(),
          axis.line=element_line(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")  
  
  ################DISPLAYS THE OTHER MAIN PANELS#########
  
  ######PANEL Parameters######
  output$title<-renderText({
    paste("The individual predicted pharmacokinetic parameters")
    
  })
  
  
  output$ip<-renderTable({
    tab1
  })
  
  output$textparameters<-renderUI({
    tags$i(" Values are mean of the conditional distribution [95% confidence interval]",style = "font-size: 12px;")
  })
  
  ######PANEL AUC########
  output$auc<-renderValueBox({
    valueBox(
    HTML(paste(tags$b(round(AUC, 2)),"[",lower_bound,";",upper_bound,"]")),HTML("<b>The predicted AUC24h of levofloxacine in mg*h/L</b><br>Values are mean of the conditional distribution [95% confidence interval]"),color="orange")
  })

  ######PANEL Prediction errors########
  output$tabconc<-renderTable({
    TAB
  })
  
  output$textpe<-renderUI({
    div(tags$i("The individual predicted concentration is the Maximum A Posteriori (MAP) estimate.",style = "font-size: 12px;"),
        tags$br(),
        tags$i("Relative prediction error= (Observed concentration-Predicted concentration)/Predicted concentration",style = "font-size: 12px;"))
    
  })
  
  ######PANEL Plot######
  output$Predind<-renderPlot({
    obj
    
  })
  
  
  output$textplot<-renderUI({
    tagList(
      div(paste("The continuous curve corresponds to the individual predicted concentrations using conditional mean estimates.",
                "The dashed curve corresponds to the population predicted concentrations.",
                "The blue ribbon corresponds to the 95% prediction interval.",
                "The red vertical dashed lines represent the optimal sampling times.")),
      
      tags$br(),
      
      div(tags$i("PS : The individual predicted concentrations represented are the median of all concentrations predicted per time for each sample of the conditional distribution.",style = "font-size: 12px;"),
          tags$i("The lower part of the 95% prediction interval corresponds to the 2.5 percentile of all concentrations predicted per time for each sample of the conditional distribution.",style = "font-size: 12px;"),
          tags$i("The upper part of the 95% prediction interval corresponds to the 97.5 percentile of all concentrations predicted per time for each sample of the conditional distribution.",style = "font-size: 12px;")))
  })
  
  #####PANEL conditional distributions plots######
  
  output$conddist<-renderPlot({
    plot_distribcond2
    
  })
  
  output$textconddist<-renderUI({
    tagList(
      div(paste("The pink distributions corresponds to the conditional distributions of the parameters for your patient.",
                "The other distributions correpond to the conditional distributions of the parameters for 5 simulated subjects under the Lemaitre modified model."),
          
          tags$br(),
          div(tags$p("WARNING : If your patient's conditional distributions of parameters differ significantly from the simulated patients, it may indicate an abnormal pharmacokinetic profile and therefore a sampling error.",style="color:red;"))
      ))
  })
})
}
  

# Run the application 
shinyApp(ui = ui, server = server)
