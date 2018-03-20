## ----setup1, message=FALSE-----------------------------------------------
if (!require(gClinBiomarker)) { 
  install_github("RPackages/gClinBiomarker",  host="https://github.roche.com/api/v3")
  library(gClinBiomarker)
}

library(knitr)
library(devtools)
library(ggplot2)
data(input)
sample.data <- input

## ------------------------------------------------------------------------
head(sample.data)

str(sample.data)

## ------------------------------------------------------------------------
kable(
SummaryVars(data=sample.data,trt='Arm', subgroup='BEP', var=c('Age','Sex'), 
       var.class=c("numeric","categorical"))
)

## ------------------------------------------------------------------------
kable(
SummaryVars(data=sample.data,trt='Arm', subgroup='BEP', var=c('Age','Sex'), 
       var.class=c("numeric","categorical"), compare.subgroup=TRUE)
)

## ------------------------------------------------------------------------
kable(
SummaryVars(data=sample.data,trt='Arm', subgroup='BEP', var=c('Age','Sex'), 
       var.class=c("numeric","categorical"), test.subgroup=TRUE)
)

## ------------------------------------------------------------------------
kable(
SummaryVars(data=sample.data,trt='Arm', subgroup='BEP', var=c('Age','Sex'), 
       var.class=c("numeric","categorical"), trt.order=c("TRT","CTRL"))
)

## ------------------------------------------------------------------------
kable(
SummaryVars(data=sample.data, subgroup='BEP', var=c('Age','Sex'), 
       var.class=c("numeric","categorical"))
)

## ------------------------------------------------------------------------
sample.data$BEP2 <- ifelse(sample.data$BEP==1,"Yes","No")
kable(
SummaryVars(data=sample.data,trt='Arm', subgroup='BEP2', var=c('Age','Sex'), 
       var.class=c("numeric","categorical"), subgroup.indicator="Yes")
)

## ----include=TRUE,fig.width=10, fig.height=5-----------------------------
PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric", log2=TRUE)

## ----include=TRUE,fig.width=10, fig.height=5-----------------------------
PlotProperty(data=input, biomarker.var=NULL, biomarker.class=NULL,
             var=c("Weight","Age"), var.class=c("numeric", "numeric"),
             log2=c(TRUE, FALSE), par.param = list(mfrow=c(1,2), cex.axis=1.2))

## ----include=TRUE, fig.width=12, fig.height=12---------------------------
PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric",
             var=c("ECOG", "Country"), var.class=c("categorical", "categorical"),
             log2=TRUE, par.param = list(mfrow=c(3,2), cex.axis=1.2),
             show.clinical.uni=TRUE)

## ----include=TRUE, fig.width=10, fig.height=12---------------------------
PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric",
             var=c("Sex", "Age"), var.class=c("categorical", "numeric"),
             log2=c(TRUE, FALSE, FALSE), par.param = list(mfrow=c(3,2), cex.axis=1.4),
             show.clinical.uni=TRUE)

## ----include=TRUE, fig.width=14, fig.height=10---------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var="Country", var.class="categorical", 
             par.param = list(mfrow=c(2,2), cex.axis=1.2),
             show.clinical.uni=TRUE)

## ----include=TRUE, fig.width=14, fig.height=12---------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(mfrow=c(3,2), cex.axis=1.2),
             show.clinical.uni=TRUE)

## ----include=TRUE, fig.width=14, fig.height=10---------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(mfrow=c(2,2), cex.axis=1.2), show.biomarker.uni = FALSE,
             show.clinical.uni=TRUE)

## ----include=TRUE, fig.width=14, fig.height=10---------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(mfrow=c(2,2), cex.axis=1.2))

## ----include=TRUE, fig.width=14, fig.height=10---------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(mfrow=c(2,2), cex.axis=1.2), show.association = FALSE,
             show.clinical.uni=TRUE)

## ----include=TRUE, fig.width=14, fig.height=6----------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(mfrow=c(1,2), cex.axis=1.2), 
             show.biomarker.uni = FALSE)

## ----include=TRUE, fig.width=14, fig.height=6----------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(mfrow=c(1,2), cex.axis=1.1),
             show.biomarker.uni = FALSE, show.association = FALSE, 
             show.clinical.uni = TRUE)

## ----include=TRUE, fig.width=14, fig.height=6----------------------------
PlotProperty(data=input, biomarker.var="KRAS.mutant", biomarker.class="categorical",
             var=c("Country", "Age"), var.class=c("categorical", "numeric"),
             par.param = list(cex.axis=1.2),
             show.association = FALSE)

## ---- fig.width=4, fig.height=7------------------------------------------
PlotRspBar (input, outcome.var="Response",
            rsp.levels=c("CR", "PR","SD","NON CR/PD", "PD","NE"))

## ---- fig.width=4, fig.height=7------------------------------------------
PlotRspBar (input, outcome.var="Response", 
            binary=TRUE,
            rsp.response = c("CR","PR"),
            rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA))

## ---- fig.width=4, fig.height=7------------------------------------------
PlotRspBar (input, outcome.var="Response", 
            binary=TRUE,
            rsp.response = c("CR","PR"),
            rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
            col=c("green","orange"))

## ---- fig.width=7, fig.height=4------------------------------------------
PlotRspBar (input, outcome.var="Response", 
            binary=TRUE,
            rsp.response = c("CR","PR"),
            rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
            horiz=TRUE)

## ----fig.width=7, fig.height=7-------------------------------------------
PlotRspBar (input, outcome.var="Response", 
            binary=FALSE,
            rsp.levels=c("CR", "PR","SD","NON CR/PD", "PD","NE"),
            trt="Arm")

PlotRspBar (input, outcome.var="Response", 
            binary=TRUE,
            rsp.response = c("CR","PR"),
            rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
            trt="Arm")

## ----fig.width=7, fig.height=7-------------------------------------------

PlotRspBar (input, outcome.var="Response", 
            binary=TRUE,
            rsp.response = c("CR","PR"),
            rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
            trt="Arm", compare.bep=TRUE,bep="BEP")



## ----fig.width=7, fig.height=7-------------------------------------------
PlotRspBar (input, outcome.var="Response", 
            binary=TRUE,
            rsp.response = c("CR","PR"),
            rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
            trt="Arm", compare.var=TRUE,var="KRAS.mutant")

## ----fig.width=7, fig.height=7-------------------------------------------
PlotRspBar (input, outcome.var="Response", 
            binary=FALSE,
            rsp.levels=c("CR", "PR","SD","NON CR/PD", "PD","NE"),
            trt="Arm", 
            compare.var=TRUE,var="Sex", plot.count = TRUE)


## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS ITT"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, bep="BEP",
             tte="PFS",cen="PFS.event", main="PFS BEP"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS ITT by treatment", trt="Arm"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS ITT by treatment", trt="Arm",
             col=c("orange","brown"),lty=c(1,2)))


## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS ITT by treatment", 
       trt="Arm", 
       plot.grid = FALSE, 
       plot.median=T))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="OS BEP by treatment, by KRAS mutation", 
              var="KRAS.mutant"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="OS BEP by treatment, by KRAS mutation", 
              var="KRAS.mutant", col=c("skyblue","darkgray")))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="OS BEP by treatment, by KRAS mutation", 
             trt="Arm", var="KRAS.mutant"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS BEP by treatment, by KRAS mutation", 
             trt="Arm", var="KRAS.mutant",
             col=c("orange","orange","brown","brown"),
             lty=c(3,1,3,1)))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS BEP by treatment, by KRAS expression", 
             trt="Arm", var="KRAS.exprs", var.class="numeric",
             percentile.cutoff=0.5,xlim=c(0,18))
             )

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS BEP by treatment, by KRAS expression", 
             trt="Arm", var="KRAS.exprs", var.class="numeric",
             numerical.cutoff=100)
             )

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS BEP by treatment, by KRAS expression", 
             trt="Arm", var="KRAS.exprs", var.class="numeric",
             numerical.cutoff=100, equal.in.high = F)
             )

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event",  
             main="PFS BEP by treatment, by KRAS expression", 
             trt="Arm", var="KRAS.exprs", var.class="numeric",
             numerical.cutoff=c(100,500), xlim=c(0, 20))
             )

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             main="PFS BEP by treatment, by KRAS expression", 
             trt="Arm", var="KRAS.exprs", var.class="numeric",
             numerical.cutoff=c(100,500),
             col=c("green","green","green","brown","brown","brown"),
             xlim=c(0,20))
             )

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="OS",cen="OS.event", 
             main="OS BEP by treatment, by KRAS mutation", 
             varlist=c("Arm","KRAS.mutant"), 
             varlist.levels=list(c("TRT","CTRL"),c("Wild Type","Mutant")),
             legend.loc="left"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="OS",cen="OS.event", bep="BEP",
             main="OS BEP by treatment, by KRAS mutation", varlist=c("Arm","KRAS.mutant"), 
             varlist.levels=list(c("TRT","CTRL"),c("Wild Type","Mutant")),
             varlist.labels=list(c("trt","ctrl"),c("wt","mut")),
               plot.median=T,legend.loc="left"))

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
print(PlotKM(data=sample.data, tte="PFS",cen="PFS.event", 
             var=c("Arm","KRAS.mutant"), 
             bep="BEP", legend.loc=NULL, legend.x=5, legend.y=.8))

## ---- fig.height=5,fig.width=7-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical")

## ---- fig.height=10,fig.width=9------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical",
                                  tabforest = TRUE)

## ---- fig.height=5,fig.width=7-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical",
                                  bep = 'BEP', 
                                  bep.indicator=1)

## ---- fig.height=4,fig.width=6-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical", 
                                  var.name="KRAS mut",
                                  show.itt=FALSE, 
                                  show.bep=FALSE)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                              #  cols=c("black","black","darkgreen","darkgreen","darkgreen"),
                                  numerical.cutoff=NULL,
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  numerical.cutoff=NULL,
                                  greater=FALSE, less=TRUE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,
                                  bep = 'BEP', bep.indicator=1)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  numerical.cutoff=NULL,
                                  greater=TRUE, less=TRUE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  numerical.cutoff=NULL,
                                  greater=TRUE, less=TRUE, greater.by.less = TRUE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  numerical.cutoff=NULL,
                                  greater=FALSE, less=FALSE,
                                  within.bin=TRUE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=c(50,100),
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=c(50,100),
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,
                                  bep = 'BEP', bep.indicator=1, covariate="Age")

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=c(50,100),
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,
                                  bep = 'BEP', bep.indicator=1, strata="Sex")

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=c(50,100),
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,
                                  bep = 'BEP', bep.indicator=1,covariate="Age",strata="Sex")

## ---- fig.height=4,fig.width=6-------------------------------------------
PlotTabForestBiomarker(data=subset(input, Arm=="TRT"),
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=c(50,100),
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,covariate="Age",strata="Sex")

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  numerical.cutoff=NULL,
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,
                                  covariate=NULL, #Sex
                                  strata=NULL, #Age
                                  placebo.code='TRT',
                                  active.code='CTRL')

## ---- fig.height=12,fig.width=9------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric", var.name="KRAS exprs",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  numerical.cutoff=NULL,
                                  greater=TRUE, less=FALSE,
                                  within.bin=FALSE,
                                  show.itt=TRUE, show.bep=TRUE,
                       across.and.within = TRUE)

## ---- fig.height=5,fig.width=7-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="binary",
                                  outcome.var=c("Response"),
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical")

## ---- fig.height=10,fig.width=9------------------------------------------
PlotTabForestBiomarker(data=input,                                  
                                  outcome.class="binary",
                                  outcome.var=c("Response"),
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical",
                                  tabforest = TRUE)

## ---- fig.height=5,fig.width=7-------------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="continuous",
                                  outcome.var=c("Lab_ontrt"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical")

## ---- fig.height=10,fig.width=10-----------------------------------------
PlotTabForestBiomarker(data=input,
                                  outcome.class="continuous",
                                  outcome.var=c("Lab_ontrt"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical",tabforest=T)

## ---- fig.height=12,fig.width=9------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1,
                   compare.bep.itt=TRUE
                   )

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant"
                   )

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant", show.itt=TRUE
                   )

## ---- fig.height=12,fig.width=9------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant", show.itt=TRUE, show.bep=TRUE
                   )

## ---- fig.height=12,fig.width=9------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Age"), compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant", show.itt=TRUE, show.bep=TRUE
                   )

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestMulti(data=subset(input,Arm=="TRT"),
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=TRUE
                   )

## ---- fig.height=12,fig.width=9------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm", percentile.cutoff=c(.33,.66),
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant"
                   )

## ---- fig.height=12,fig.width=9------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm", percentile.cutoff=c(.33,.66),less=TRUE,greater=FALSE,
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant"
                   )

## ---- fig.height=8,fig.width=9-------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm", percentile.cutoff=c(.33,.66),within.bin=TRUE,
                                  var=c("Sex","Age"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=FALSE,
                   compare.subgroup=TRUE, subgroup="KRAS.mutant"
                   )

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   compare.bep.itt=FALSE, compare.subgroup=FALSE,itt.name=""
                   )

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestMulti(data=input,
                                                                   outcome.class="binary",
                                  outcome.var=c("Response"),
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   compare.bep.itt=FALSE, compare.subgroup=FALSE,itt.name=""
                   )

## ---- fig.height=10,fig.width=10-----------------------------------------
PlotTabForestMulti(data=input,
                                                                   outcome.class="binary",
                                  outcome.var=c("Response"),
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   compare.bep.itt=FALSE, compare.subgroup=FALSE,itt.name="", tabforest=T
                   )

## ---- fig.height=6,fig.width=9-------------------------------------------
PlotTabForestMulti(data=input,
                               outcome.class="continuous",
                                  outcome.var=c("Lab_ontrt"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   compare.bep.itt=FALSE, compare.subgroup=FALSE,itt.name=""
                   )

## ---- fig.height=10,fig.width=10-----------------------------------------
PlotTabForestMulti(data=input,
                 outcome.class="continuous",
                                  outcome.var=c("Lab_ontrt"),
                                  trt="Arm",
                                  var=c("Sex","Age"),
                   compare.bep.itt=FALSE, compare.subgroup=FALSE,itt.name="", tabforest=T
                   )

## ------------------------------------------------------------------------
PlotSTEPP(data = input,
          outcome.var = c("PFS", "PFS.event"),
          outcome.class = "survival",
          trt = "Arm",
          var = "KRAS.exprs",
          placebo.code = "CTRL",
          active.code = "TRT",
          csv.name = NULL,
          pdf.name = NULL
) 

## ------------------------------------------------------------------------
PlotSTEPP(data = input,
          outcome.var = "Baseline.SLD",
          outcome.class = "continuous",
          trt = "Arm",
          var = "KRAS.exprs",
          covariate= "Sex",
          placebo.code = "CTRL",
          active.code = "TRT",
          csv.name = NULL,
          pdf.name = NULL
)

## ------------------------------------------------------------------------
PlotSTEPP(data = input,
          outcome.var = "ECOG",
          outcome.class = "binary",
          trt = "Arm",
          var = "KRAS.exprs",
          placebo.code = "CTRL",
          active.code = "TRT",
          csv.name = NULL,
          pdf.name = NULL
)

## ------------------------------------------------------------------------
CoxTab(data=sample.data, tte="OS", cens="OS.event",bep='BEP', var='Sex' )



## ------------------------------------------------------------------------
kable(
CoxTab(data=sample.data,tte="OS", cens="OS.event",  var=c('Sex',"Country","Age"), 
       var.class=c("categorical","categorical","numeric"))
)

## ------------------------------------------------------------------------
kable(
CoxTab(data=sample.data,tte="OS", cens="OS.event",  var=c('Sex',"Country","Age"))
)


## ------------------------------------------------------------------------
kable(
CoxTab(data=sample.data,tte="OS", cens="OS.event",  var=c('Sex',"Country","Age"), 
       var.class=c("categorical","categorical","numeric"), bep="BEP")
)

## ------------------------------------------------------------------------
kable(
CoxTab(data=sample.data,tte="OS", cens="OS.event",  var=c('Sex',"Country","Age"), 
       var.class=c("ordered.factor","categorical","numeric"), 
       ordered.factor.levels.list=list(Sex=c("M","F")),bep="BEP")
)

## ------------------------------------------------------------------------
kable(
CoxTab(data=sample.data,tte="OS", cens="OS.event",  var=c('Sex',"Country","Age"),
       additive=FALSE)
)


## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input$OS, treatment.var = input$Arm, 
            placebo.code = "CTRL", active.code = "TRT", 
            outcome.class = "continuous")

## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input$OS, treatment.var = input$Arm, 
            placebo.code = "CTRL", active.code = "TRT", 
            outcome.class = "continuous", covariate.var = input$Sex)

## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input$OS, treatment.var = input$Arm, 
            placebo.code = "CTRL", active.code = "TRT", 
            outcome.class = "continuous", covariate.var = input$Sex, 
            return.fit = TRUE)

## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input$OS.event, treatment.var = input$Arm, 
            placebo.code = "CTRL", active.code = "TRT", 
            outcome.class = "binary")

## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input[, c("OS", "OS.event")], 
            treatment.var = input$Arm, placebo.code = "CTRL", 
            active.code = "TRT", outcome.class = "survival")

## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input[, c("OS", "OS.event")], 
            treatment.var = input$Arm, placebo.code = "CTRL", 
            active.code = "TRT", outcome.class = "survival", 
            covariate.var = input$Sex)

## ------------------------------------------------------------------------
SummaryTwoGroups(outcome.var = input[, c("OS", "OS.event")], 
            treatment.var = input$Arm, placebo.code = "CTRL", 
            active.code = "TRT", outcome.class = "survival", 
            covariate.var = input$Sex, strat.factor.var = input$Age)

## ------------------------------------------------------------------------
  kable(
    LogRankTab(data=input,tte="PFS",cens="PFS.event",var="Arm")
  )

## ------------------------------------------------------------------------
example <- data.frame( y=c(rnorm(30)+10, rnorm(4)+20, rnorm(15)+15, NA), 
 							time=c(rep("t2", 30), rep("t4",4), rep("t1", 15), "t3"),
 							grp=sample(1:3, 50, TRUE), sex=sample(1:2, 50, TRUE))
head(example)
str(example)

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(list(a=rnorm(50,,3), b=rnorm(25,1,4), c=rnorm(75,2,1)))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(list(a=rnorm(50,,3), b=rnorm(25,1,4), c=rnorm(75,2,1)), horizontal=TRUE,
 			 Xaxis=list(las=2, hadj=2), Xaxis2=list(las=2, hadj=-.25))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(matrix(rnorm(100), ncol=4, dimnames=list(NULL, LETTERS[1:4])) )

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(rnorm(50,,3), rnorm(25,1,2), rnorm(75,2,1))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(rnorm(50,,3), rnorm(25,1,2), rnorm(75,2,1), horizontal=TRUE)

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example, y~time, sc.pch=16)

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example, y~time, sc.pch=16, box.type="bp")

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example, y~time, sc.pch=16, box.type="bp", Title=list(main="Custom Main Title"))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example, y~time, trend="median", ylab="Y-Axis Label")

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example, y~time, trend="median", ylab="Y-Axis Label",
           sc.col=c("red", "blue", "green")[example$grp] )

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example, y~time, trend="median", ylab="Y-Axis Label",
           sc.col=c("red", "blue", "green")[example$grp],
 		     sc.pch=c(5,10)[example$sex] )
 
legend("bottomright", fill=c("red", "green", "blue", "white", "white"), 
           legend=c("Stage I", "Stage II", "Stage III", "Female", "Male"), 
           pch=c(-1, -1, -1, 5, 10), border=NA)

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
example2 <- data.frame(y1=12+rnorm(15), y2=15+rnorm(15), y3=25+rnorm(15))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example2, var=c("y1", "y2", "y3"), Grid=TRUE, trend="mean")

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(example2, var=c("y1", "y2", "y3"), box.type="bp", Grid=TRUE)

## ----include=TRUE,fig.width=7, fig.height=4------------------------------

BoxPlot(rnorm(1856, 5), runif(1245, 2,10), exp(rnorm(2311)), sc.col=as.rgb("black", .05), 
            box.type="bp", Xaxis=list(labels=c("~N(5,1)", "~Unif(2,10)", "~exp(N(0,1))")), ylim=c(0,15))     

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
    mat <- matrix(c(rep(10,45), rep(15,45), rep(35,45))+rnorm(135), ncol=3)
    BoxPlot(mat, trend="mean", trend.col="red", Ylabel=list(text="Example Measurements"))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
data(mtcars)
BoxPlot(mtcars, mpg~gear:cyl)  

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(mtcars, mpg~cyl %in% gear)     

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
dat <- mtcars
dat$cyl <- factor(dat$cyl, levels=c(4,6,8), labels=c("4Cyl", "6Cyl", "8Cyl"))
dat$gear <- factor(dat$gear, levels=c(3,4,5), labels=c("3Gear", "4Gear", "5Gear")) 
BoxPlot(dat, mpg~cyl %in% gear)

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(dat, mpg~cyl:gear, XaxisTab=list(),mar=c(8,3,5,1))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(mtcars, mpg~cyl:gear, XaxisTab=list(), mar=c(5,8,5,4), 
 			 horizontal=TRUE, Ylabel=list(text="Y-axis label now appearing on X-axis"))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(dat, mpg~cyl:gear, XaxisTab=list(font=2, col="darkblue", cex=1.25), mar=c(5,3,5,1))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(dat, mpg~cyl:gear, XaxisTab=list(Label=list(font=2, col="darkblue", cex=1.25), 
 			 Text=list(col="red")), mar=c(5,3,5,1))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(dat, mpg~cyl:gear:vs, XaxisTab=list(), mar=c(5,3,5,1))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(mtcars, mpg~cyl %in% gear, 
 			 Xaxis=list(labels=paste(rep(c("4Cyl", "6Cyl", "8Cyl"),3), 
 						c(rep("3Gear",3), rep("4Gear",3), rep("5Gear",3)), sep=".")))

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(dat, mpg~cyl %in% gear, Title=list(main="Miles per Gallon by Number of Gears", 
 			 col.main="Green", cex.main=2.5), vline=c(3.5, 6.5), vl.lty=2, vl.col="gray", vl.lwd=2,
           Xaxis=list(labels=NA, at=1:9, tick=TRUE), col=c(rep("blue", 3), rep("red", 3), rep("green", 3)),
           Xaxis2=list(tick=FALSE), Yaxis=list(at=seq(10,34,2)), Grid=list(x=1:9, y=seq(10,34,2)),
           Xlabel=list(text=paste(rep(c("4Cyl", "6Cyl", "8Cyl"),3), c(rep("3Gear",3), rep("4Gear",3), rep("5Gear",3))),
                       at=1:9, las=2, adj=1, line=0.75, col=c(rep("blue", 3), rep("red", 3), rep("green", 3))), 
           mean.col=c(rep("cyan", 3), rep("orange", 3), rep("magenta", 3)), Box=FALSE, trend="mean",mar=c(6,4,5,2) )

## ----include=TRUE,fig.width=7, fig.height=4------------------------------
BoxPlot(mtcars, mpg~cyl %in% gear, Title=list(main="Miles per Gallon by Number of Gears", col.main="#84d52b", cex.main=1.5), 
 		vline=c(3.5, 6.5), vl.lty=2, vl.col="gray", vl.lwd=2,
 		Xaxis2=list(tick=FALSE, las=2, hadj=-.25), Yaxis=list(at=seq(10,34,2)), Grid=list(x=1:9, y=seq(10,34,2)),
 		mean.col=c(rep("cyan", 3), rep("orange", 3), rep("magenta", 3)), Box=FALSE, trend="mean",
 		mar=c(3, 7, 4, 4), horizontal=TRUE, sc.pch=c(0, 15)[dat$am+1], sc.col="wheat4",
 		XaxisTab=list(Text=list(col=c(rep("cyan", 3), rep("orange", 3), rep("magenta", 3)))) )
 
 legend(x="topright", pch=c(0, 15), legend=c("automatic", "manual"), box.lty=0, col="wheat4")

## ------------------------------------------------------------------------
PlotLong(longbmkr, aes(x=vm, y=ep),
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Mean +/- SEM')

## ------------------------------------------------------------------------
PlotLong(longbmkr, aes(x=vm, y=ep), fun.data = 'tukey',
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Tukey Whiskers and Hinges')

## ------------------------------------------------------------------------
PlotLong(longbmkr, aes(x=vm, y=ep, group=trt, color=trt, fill=trt), 
         fun.data = 'tukey', facet.fun = . ~ sex,
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Biomarker Timecourse',
         labs.caption = 'Ribbons represent Tukey hinges and whiskers')

## ---- fig.height = 3.5,eval=TRUE-----------------------------------------
library(dplyr)
PlotLong(longbmkr %>% filter(vm <= 24), 
         aes(x=vm, y=ep, group=trt, color=trt, fill=trt), 
         fun.data = 'tukey', facet.fun = . ~ sex,
         show.counts = 'table',
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Biomarker Timecourse',
         labs.caption = 'Ribbons represent Tukey hinges and whiskers')

## ---- fig.height = 3.5, eval=TRUE----------------------------------------
PlotLong(longbmkr %>% filter(vm <= 24), 
         aes(x=vm, y=ep, group=trt, color=trt, fill=trt), 
         fun.data = 'tukey', facet.fun = . ~ sex,
         show.counts = 'table',
         plot.style = 'errorbars',
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Biomarker Timecourse',
         labs.caption = 'Ribbons represent Tukey hinges and whiskers')

