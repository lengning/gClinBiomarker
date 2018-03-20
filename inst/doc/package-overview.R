## ----setup, include=FALSE------------------------------------------------
if (!require(gClinBiomarker)) { 
  install_github("lengning/gClinBiomarker")
  library(gClinBiomarker)
}

library(knitr)

## ---- warning=FALSE, message=FALSE, eval=FALSE---------------------------
#  library(devtools)
#  install_github("lengning/gClinBiomarker")

## ----warning=FALSE, message=FALSE, eval=FALSE----------------------------
#  library(gClinBiomarker)

## ------------------------------------------------------------------------
head(input)
str(input)

## ------------------------------------------------------------------------
  SummaryVars(data=input, trt='Arm', subgroup='BEP', var=c('Sex','Age'), 
       var.class=c('categorical','numeric'))

## ------------------------------------------------------------------------
kable(
  SummaryVars(data=input, trt='Arm', subgroup='BEP', var=c('Sex','Age'), 
       var.class=c('categorical','numeric'))
)

## ---- fig.height=5, fig.width=10-----------------------------------------
CompareKM(data=input,tte='PFS', cen='PFS.event',trt='Arm', bep='BEP')

## ---- fig.height=7,fig.width=10------------------------------------------
PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Weight"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=TRUE
                   )

## ------------------------------------------------------------------------
input.bep <- subset(input, BEP==1) # dataset with only BEP samples

kable(
  CoxTab(data=input,tte="PFS", cens="PFS.event",  var='Arm', 
       var.class="categorical"),
caption="full population, unadjusted"
)

kable(
  CoxTab(data=input,tte="PFS", cens="PFS.event",  var=c('Arm','Sex',"Weight"), 
       var.class=c("categorical","categorical","numeric")),
caption="full population, adjusted for Sex, Weight"
)

kable(
  CoxTab(data=input.bep,tte="PFS", cens="PFS.event",  var='Arm', 
       var.class="categorical"),
caption="BEP, unadjusted"
)
kable(
  CoxTab(data=input.bep,tte="PFS", cens="PFS.event",  var=c('Arm','Sex',"Weight"), 
       var.class=c("categorical","categorical","numeric")),
caption="BEP, adjusted for Sex, Weight"
)


## ----fig.width=12, fig.height=8------------------------------------------
PlotProperty(data=input, biomarker.var="KRAS.exprs", biomarker.class="numeric",
             var=c("Sex", "Country"), var.class=c("categorical", "categorical"),
             log2=TRUE, par.param = list(mfrow=c(2,3)))

## ----fig.width=8, fig.height=4-------------------------------------------
input.ctrl <- subset(input, Arm=="CTRL") ## Data with only ctrl samples
res.multicut.ctrl <- PlotTabForestBiomarker(data=input.ctrl,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), main.prefix="CTRL",
                                  greater=TRUE, less=FALSE)

input.trt <- subset(input, Arm=="TRT") ## Data with only trt samples
res.multicut.trt <- PlotTabForestBiomarker(data=input.trt,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), main.prefix="TRT", 
                                  greater=TRUE, less=FALSE)

## ---- fig.height=5,fig.width=10------------------------------------------
res.multicut <- PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  greater=TRUE, less=FALSE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=5,fig.width=10------------------------------------------
res.withinbin <- PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  within.bin=TRUE,
                                  show.itt=TRUE, show.bep=TRUE)

## ------------------------------------------------------------------------
stepp.out <- PlotSTEPP(data = input,
          outcome.var = c("PFS", "PFS.event"),
          outcome.class = "survival",
          trt = "Arm",
          var = "KRAS.exprs",
          placebo.code = "CTRL",
          active.code = "TRT"
) 

## ---- fig.height=5,fig.width=8-------------------------------------------
res.subgroups.cut <- PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=100, within.bin=TRUE,
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=5,fig.width=10------------------------------------------
input$KRAS.exprs.2group <- ifelse(input$KRAS.exprs >= 100, ">=100","<100")
input$KRAS.exprs.2group <- factor(input$KRAS.exprs.2group, levels=c(">=100","<100")) # ">=100" as Dx+
res.subgroups.bi <- PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.exprs.2group", 
                                  var.class="categorical",
                                  show.itt=TRUE, show.bep=TRUE)

## ----include=TRUE,fig.width=9, fig.height=9------------------------------
PlotKM(data=input, tte="PFS",cen="PFS.event", bep="BEP", 
             main="PFS BEP by treatment, by KRAS expression subgroups", 
             var=c("Arm","KRAS.exprs.2group"), legend.loc="topright",
       plot.median=TRUE)

## ------------------------------------------------------------------------
input.bep <- subset(input, BEP==1)
kable(
SummaryVars(data=input.bep,trt='Arm', subgroup='KRAS.exprs.2group', var=c('Sex'), 
       var.class="categorical", subgroup.indicator=">=100",compare.subgroup=TRUE, subgroup.name="exprs")
)

## ---- fig.height=7,fig.width=10------------------------------------------
res.subgroup.cov <- PlotTabForestMulti(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var=c("Sex","Weight"),
                   compare.bep.itt=FALSE,
                   compare.subgroup=TRUE,
                   subgroup="KRAS.exprs.2group"
                   )

## ------------------------------------------------------------------------
input.trt <- subset(input, Arm=="TRT")

## ------------------------------------------------------------------------
kable(
  SummaryVars(data=input.trt, trt=NULL, subgroup='BEP', var=c('Age','Sex'), 
       var.class=c('numeric','categorical'))
)

## ---- fig.height=5, fig.width=5------------------------------------------
CompareKM(data=input.trt,tte='PFS', cen='PFS.event',trt=NULL, bep='BEP')

## ---- fig.height=7,fig.width=10------------------------------------------
PlotTabForestMulti(data=input.trt,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var=c("Sex","Weight"),
                   bep="BEP",bep.indicator=1, compare.bep.itt=TRUE
                   )

## ------------------------------------------------------------------------
input.trt.bep <- subset(input.trt, BEP==1) # dataset with only BEP samples

kable(
  CoxTab(data=input.trt,tte="PFS", cens="PFS.event",  var=c('Sex',"Weight"), 
       var.class=c("categorical","numeric")),
caption="the full population, model of Sex, Weight"
)
kable(
  CoxTab(data=input.trt.bep,tte="PFS", cens="PFS.event",  var=c('Sex',"Weight"), 
       var.class=c("categorical","numeric")),
caption="BEP, model of Sex, Weight"
)

## ----fig.width=7, fig.height=4-------------------------------------------
PlotProperty(data=input.trt, biomarker.var="KRAS.exprs", biomarker.class="numeric",
             var=c("Sex", "Country"), var.class=c("categorical", "categorical"),
             log2=TRUE, par.param = list(mfrow=c(1,2)))

## ---- fig.height=5,fig.width=8-------------------------------------------
res.multicut <- PlotTabForestBiomarker(data=input.trt,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  greater=TRUE, less=FALSE,
                                  show.itt=FALSE, show.bep=FALSE)

## ---- fig.height=4,fig.width=8-------------------------------------------
res.subgroups.cut <- PlotTabForestBiomarker(data=input.trt,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=NULL,
                                  numerical.cutoff=100, 
                                  greater=TRUE,less=FALSE,
                                  show.itt=FALSE, show.bep=FALSE)

## ---- fig.height=4,fig.width=9-------------------------------------------
input.trt$KRAS.exprs.trt.2group <- ifelse(input.trt$KRAS.exprs >= 100, ">=100","<100")
input.trt$KRAS.exprs.2group <- factor(input.trt$KRAS.exprs.2group, levels=c(">=100","<100")) # ">=100" as Dx+

res.subgroups.bi <- PlotTabForestBiomarker(data=input.trt,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var="KRAS.exprs.trt.2group", 
                                  var.class="categorical",
                                  show.itt=FALSE, show.bep=FALSE)

## ----include=TRUE,fig.width=7, fig.height=7------------------------------
PlotKM(data=input.trt, tte="PFS",cen="PFS.event", bep="BEP", 
             main="PFS in BEP, by KRAS expression subgroups, single arm", 
             var=c("KRAS.exprs.trt.2group"), legend.loc="topright",
       plot.median=TRUE)

## ------------------------------------------------------------------------
input.trt.bep <- subset(input.trt, BEP==1)
kable(
SummaryVars(data=input.trt.bep,trt=NULL, subgroup='KRAS.exprs.2group', var=c('Sex'), 
       var.class="categorical", subgroup.indicator=">=100",compare.subgroup=TRUE, subgroup.name="exprs")
)

## ---- fig.height=4,fig.width=9-------------------------------------------
res.subgroups.cat <- PlotTabForestBiomarker(data=input,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt="Arm",
                                  var="KRAS.mutant", 
                                  var.class="categorical",
                                  show.itt=TRUE, show.bep=TRUE)

## ---- fig.height=4,fig.width=9-------------------------------------------
res.subgroups.cat <- PlotTabForestBiomarker(data=input.trt,
                                  outcome.class="survival",
                                  outcome.var=c("PFS","PFS.event"),
                                  trt=NULL,
                                  var="KRAS.mutant", 
                                  var.class="categorical")

## ---- fig.height=7,fig.width=7-------------------------------------------
Rsp.out <- PlotRspBar(data=input, outcome.var="Response", binary=FALSE,
rsp.levels=c("CR", "PR","SD","NON CR/PD", "PD","NE"),trt="Arm",
compare.bep.itt=TRUE, bep = "BEP")

## ------------------------------------------------------------------------
kable(Rsp.out$count,caption="count")

## ------------------------------------------------------------------------
kable(round(Rsp.out$perc,2), caption="percentage")

## ---- fig.height=7,fig.width=7-------------------------------------------
PlotRspBar(data=input, outcome.var="Response", binary=FALSE,
rsp.levels=c("CR", "PR","SD","NON CR/PD", "PD","NE"),trt="Arm",
compare.bep.itt=TRUE, bep = "BEP", plot.count=TRUE)

## ---- fig.height=7,fig.width=7-------------------------------------------
Rsp.out.2 <- PlotRspBar(data=input, outcome.var="Response", binary=TRUE,
rsp.response = c("CR","PR"),
rsp.nonresponse = c("SD", "PD","NON CR/PD","NE"),trt="Arm",
compare.bep.itt=TRUE, bep = "BEP")

## ------------------------------------------------------------------------
kable(Rsp.out.2$count,caption="count")

## ------------------------------------------------------------------------
kable(round(Rsp.out.2$perc,2), caption="percentage")

## ---- fig.height=7,fig.width=7-------------------------------------------
PlotRspBar(data=input, outcome.var="Response", binary=TRUE,
rsp.response = c("CR","PR"),
rsp.nonresponse = c("SD", "PD","NON CR/PD","NE"),trt="Arm",
compare.bep.itt=TRUE, bep = "BEP", plot.count = TRUE)

## ---- fig.height=6,fig.width=8-------------------------------------------
input.ctrl <- subset(input, Arm=="CTRL") ## Data with only ctrl samples
res.multicut.ctrl <- PlotTabForestBiomarker(data=input.ctrl,
                                  outcome.class=c("binary"),
                                  outcome.var=c("Response"),
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), main.prefix="CTRL",
                                  greater=TRUE, less=FALSE)
  

input.trt <- subset(input, Arm=="TRT") ## Data with only trt samples
res.multicut.trt <- PlotTabForestBiomarker(data=input.trt,
                                  outcome.class=c("binary"),
                                  outcome.var=c("Response"),
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), main.prefix="TRT",
                                  greater=TRUE, less=FALSE)

## ---- fig.height=6,fig.width=10------------------------------------------
res.multicut.2arm <- PlotTabForestBiomarker(data=input,
                                  outcome.class=c("binary"),
                                  outcome.var=c("Response"),
                                  trt = "Arm",
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75),
                                  greater=TRUE, less=FALSE)

## ---- fig.height=6,fig.width=10------------------------------------------
res.multicut.2arm <- PlotTabForestBiomarker(data=input,
                                  outcome.class=c("binary"),
                                  outcome.var=c("Response"),
                                  trt = "Arm",
                                  rsp.cat = TRUE,
                                  rsp.response = c("CR","PR"),
                                  rsp.nonresponse = c("SD", "PD","NON CR/PD","NE",NA),
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75),
                                  within.bin=TRUE)


## ---- fig.height=7,fig.width=7-------------------------------------------
input$KRAS.exprs.2group <- ifelse(input$KRAS.exprs >= 100, ">=100","<100")
input$KRAS.exprs.2group <- factor(input$KRAS.exprs.2group, levels=c(">=100","<100")) # ">=100" as Dx+
Rsp.out.2 <- PlotRspBar(data=input, outcome.var="Response", binary=TRUE,
rsp.response = c("CR","PR"),
rsp.nonresponse = c("SD", "PD","NON CR/PD","NE"),trt="Arm",
 bep = "BEP",compare.var=TRUE, var="KRAS.exprs.2group")

## ---- fig.height=5,fig.width=7-------------------------------------------

res.multicut <- PlotTabForestBiomarker(data=input,
                                  outcome.class="continuous",
                                  outcome.var="Lab_ontrt",
                                  trt="Arm",
                                  var="KRAS.exprs", 
                                  var.class="numeric",
                                  percentile.cutoff=c(.25,.5,.75), 
                                  greater=TRUE, less=FALSE,
                                  show.itt=TRUE, show.bep=TRUE)

## ------------------------------------------------------------------------
stepp.out <- PlotSTEPP(data = input,
          outcome.var = "Lab_ontrt",
          outcome.class = "continuous",
          trt = "Arm",
          var = "KRAS.exprs",
          placebo.code = "CTRL",
          active.code = "TRT"
) 

## ---- message=FALSE------------------------------------------------------
head(longbmkr)

## ----longbmkr_description, message=FALSE, echo=FALSE---------------------
kable(
    data.frame(
        `Column Name` = names(longbmkr),
        `English Name` = c("PatientID", "Treatment", "Sex", "Age", "Education", "Biomarker Measurement", "Visitation Month", "Endpoint Measurement"),
        `Description` = c(
            "Unique Patient Identifier",
            "Binary Classification of Treatment (1) or Control (0)",
            "Male (\"m\") or Female (\"f\") binary classification of sex",
            "Age of Patient",
            "Years of Education for the Patient",
            "Biomarker Measurement|Mock Biomarker Measurement Reading",
            "Timepoint that data was collected",
            "Mock Endpoint Measurement Reading"
        ),
        check.names = FALSE),
    caption = "Longitudinal Biomarker Sample Data"
)

## ------------------------------------------------------------------------
baseline_data <- longbmkr[longbmkr['vm'] == 0,]

kable(
    head(baseline_data),
    caption = 'longitudinal data subset by visit month to produce baseline data'
)

## ------------------------------------------------------------------------
PlotLong(longbmkr, x=vm, y=ep)

## ------------------------------------------------------------------------
PlotLong(longbmkr, x=vm, y=ep, fun.data = 'tukey',
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Tukey Whiskers and Hinges')

## ------------------------------------------------------------------------
PlotLong(longbmkr, x=vm, y=ep, color=trt, fill=trt, 
         fun.data = 'tukey',
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Biomarker Timecourse',
         labs.caption = 'Ribbons represent Tukey hinges and whiskers')

## ------------------------------------------------------------------------
PlotLong(longbmkr, x=vm, y=ep, group=trt, color=trt, fill=trt, 
         fun.data = 'tukey', facet.fun = . ~ sex,
         xlab = 'Visitation Month',
         ylab = 'Biomarker Endpoint',
         labs.title = 'Biomarker Timecourse',
         labs.caption = 'Ribbons represent Tukey hinges and whiskers')

## ------------------------------------------------------------------------
library(dplyr)

## ---- fig.height = 5-----------------------------------------------------
cutoff = 1.5

longbmkr %>%
    filter(vm <= 24) %>%
    mutate(bmkr = ifelse(bmkr > cutoff, 'Positive', 'Negative')) %>%
    mutate(trt  = ifelse(trt == 0, 'CTRL', 'TRT')) %>%
    
    PlotLong(x=vm, y=ep, color=trt, fill=trt, linetype=bmkr,
             fun.data = 'mean_se',
             show.counts = 'table',
             plot.style = 'errorbars',
             errorbar.linetype = 1,
             xlab = 'Visitation Month',
             ylab = 'Endpoint',
             labs.title = 'Biomarker Timecourse',
             labs.caption = '*ribbons represent mean value +/- standard error of the mean.')

## ---- fig.height = 5-----------------------------------------------------
cutoff = 1.5

longbmkr %>%
    filter(vm <= 24) %>%
    mutate(bmkr = ifelse(bmkr > cutoff, 'Positive', 'Negative')) %>%
    mutate(trt  = ifelse(trt == 0, 'CTRL', 'TRT')) %>%
    
    PlotLong(x=vm, y=ep, color=trt, fill=trt, linetype=bmkr,
             model = lm,
             model.formula = ep ~ sex + edu + age, 
             fun.data = 'mean_se',
             show.counts = 'table',
             plot.style = 'errorbars',
             errorbar.linetype = 1,
             xlab = 'Visitation Month',
             ylab = 'Endpoint',
             labs.title = 'Biomarker Timecourse',
             labs.caption = '*ribbons represent mean value +/- standard error of the mean.')

