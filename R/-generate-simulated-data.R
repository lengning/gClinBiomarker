# set random number seed
set.seed(1)

bmkr <-
  # build list of parameters
  list(
    patients_n = 1000,  # number of patients
    visits_n   = 10,    # number of visits
    visits_d   = 6,     # months of separation between visits
    endpt_b    = 200,   # lambda used for poisson distribution
    endpt_m    = 5,     # endpoint slope
    endpt_m_sd = 0.5    # endpoint slope standard deviation
  ) %>%

  # set up variables given initial parameter list
  (function(.) { data.frame(list(
    pid  = rep(1:.$patients_n, each=.$visits_n),
    trt  = rep(round(runif(1:.$patients_n)), each=.$visits_n),
    sex  = rep(round(runif(1:.$patients_n)), each=.$visits_n),
    vm   = rep((1:.$visits_n - 1) * .$visits_d, .$patients_n),
    ep_b = rep(rpois(.$patients_n, lambda=.$endpt_b), each=.$visits_n),
    ep_m = rep(.$endpt_m, .$patients_n * .$visits_n),
    ep_s = rep(.$endpt_m_sd, .$patients_n * .$visits_n)
  )) }) %>%

  # build timecourse data and add noise
  group_by(pid) %>%
  mutate(ep_m = rep(rnorm(1, first(ep_m) * first(trt), first(ep_s)), n())) %>%
  ungroup() %>%
  mutate(ep = ep_b + ep_m * vm * ((sex - 0.5) * 0.2 + 1)) %>%
  mutate(ep = rnorm(n(), ep, ep_b * 0.01 * (vm * 0.5 + 12)))%>%
  mutate(sex = ifelse(sex == 1, "f", "m")) %>%
  select(-ep_m, -ep_s) %>%
  sample_frac(0.85) %>%
  arrange(pid, vm) %>%
  mutate(trt = as.factor(trt),
         sex = as.factor(sex))

save(bmkr, file = file.path(getwd(), "data", "bmkr.rda"))

# reset random number seed based on system time
rm(.Random.seed, envir=globalenv())
