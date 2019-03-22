

rbimodal <- function (n, cpct, mu1, mu2, sig1, sig2) {
  y0 <- rnorm(n,mean=mu1, sd = sig1)
  y1 <- rnorm(n,mean=mu2, sd = sig2)
  flag <- rbinom(n,size=1,prob=cpct)
  y <- y0*(1 - flag) + y1*flag
}

#' Simulated longitudinal biomarker data
#'
#' Contains the following variables
#'
#' @format A data frame with 8,500 rows and 8 variables:
#' \describe{
#'   \item{pid}{Patient Identifier}
#'   \item{trt}{Treatment Arm (1, 0)}
#'   \item{sex}{Patient Sex (m, f)}
#'   \item{age}{Patient Age}
#'   \item{edu}{Patient years of education}
#'   \item{bmkr}{Baseline biomarker reading}
#'   \item{vm}{Patient visit time in months}
#'   \item{ep}{Biomarker endpoint reading}
#' }
#'
"longbmkr"

longbmkr <-
  # build list of parameters
  list(
    patients_n = 1000,  # number of patients
    visits_n   = 10,    # number of visits
    visits_d   = 6,     # months of separation between visits
    endpt_b    = 200,   # lambda used for poisson distribution
    endpt_m    = 5,     # endpoint slope
    endpt_m_sd = 1.5    # endpoint slope standard deviation
  ) %>%

  # set up variables given initial parameter list
  (function(.) { data.frame(list(
    pid  = rep(1:.$patients_n, each=.$visits_n),
    trt  = rep(round(runif(1:.$patients_n)), each=.$visits_n),
    sex  = rep(round(runif(1:.$patients_n)), each=.$visits_n),
    age  = 18 + rep(runif(1:.$patients_n) * (80 - 18), each=.$visits_n),
    edu  = 8 + rep(rpois(.$patients_n, lambda=4), each=.$visits_n),
    bmkr = rep(rbimodal(.$patients_n, 0.5, 1, 2, 0.2, 0.5), each=.$visits_n),
    vm   = rep((1:.$visits_n - 1) * .$visits_d, .$patients_n),
    ep_b = rep(rpois(.$patients_n, lambda=.$endpt_b), each=.$visits_n),
    ep_m = rep(.$endpt_m, .$patients_n * .$visits_n),
    ep_s = rep(.$endpt_m_sd, .$patients_n * .$visits_n)
  )) }) %>%

  # build timecourse data and add noise
  group_by(pid) %>%
  #mutate(ep_m = rep(rnorm(1, first(ep_m) * (1 - first(trt)), first(ep_s)), n())) %>%
  ungroup() %>%
  mutate(ep = rnorm(n(), ep_b, ep_b/10) + vm *
           (1                          * rnorm(n(), ep_m, ep_s)
         + (sex - 0.5)                 * rnorm(n(), 1,    0.2)
         + (mean(edu) - edu)/mean(edu) * rnorm(n(), 0.14, 0.05)
         + (mean(age) - age)/mean(age) * rnorm(n(), 0.18, 0.05)
         + (bmkr/mean(bmkr) ^ 3 )      * rnorm(n(), 0.7,  0.05)
         + (0.5 - trt)                 * rnorm(n(), 2,    0.2)) ) %>%
  #mutate(ep = rnorm(n(), ep, ep_b * 0.01 * (vm * 0.5 + 12)))%>%
  mutate(sex = ifelse(sex == 1, "f", "m")) %>%
  select(-ep_m, -ep_s, -ep_b) %>%
  sample_frac(0.85) %>%
  arrange(pid, vm) %>%
  mutate(trt = as.factor(trt),
         sex = as.factor(sex))

usethis::use_data(longbmkr, overwrite = TRUE)
