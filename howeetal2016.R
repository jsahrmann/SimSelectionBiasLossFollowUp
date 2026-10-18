# Replicate simulations in Howe et al. (2016)
#
# References:
# [1] [2016 Howe, Cole, et al. Selection Bias Due to Loss to Follow Up in Cohort Studies. Epidemiology](doi.org/10.1097/eDe.0000000000000409)
# [2] [Use dynamic name for new column/variable in `dplyr`])https://stackoverflow.com/a/26003971)
# [3] [checking length of (possibly NULL) data-masked argument for flow control in R, dplyr](https://stackoverflow.com/a/79793270)
# [4] [What is data-masking and why do I need `{{`?](https://rlang.r-lib.org/reference/topic-data-mask.html)
# [5] [Defusing R expressions](https://rlang.r-lib.org/reference/topic-defuse.html)
#
# Date created: 20250606
# Last modified: 20260901
# John M Sahrmann


# * Setup

library(dplyr)
library(furrr)
library(ggplot2)
library(patchwork)
library(progressr)
library(logbin)
library(rlang)
library(survival)
library(tidyr)
library(webpea)


# * Constant definitions

m <- 500         # number of simulations
n <- 1000        # observations per simulation


# ** Figure 1 / Figures A1.1-A1.5

mu_t <- 200.9    # mean of exponential distribution for death times
mu_c <- 66.9     # mean of exponential distribution for censor times
prob_idu <- 0.5  # probability of injection drug use (exposure, A)
prob_hau <- 0.5  # probability of heavy alcohol use (covariate, L)


# ** Figure 2 / Figure A1.6

# proportion AA (exposure, A) in UNC HIV cohort (Table p 6/96)
prob_aa <- 0.66
# probability of composite (covariate, L) given AA, a "single binary
# composite variable (1=yes, 0=no) representing the net effect of CD4
# cell count, a prior AIDS-defining illness diagnosis, HIV RNA level,
# and insurance status at study entry on time to death" (p S5); other
# than "African Americans were assumed to be twice as likely to have a
# value of 1 for the composite variable" (p S5), no information is
# given on its expected frequency; this value (and the resulting
# vector of joint probabilities) produces a weighted survival curve
# that looks reasonably close to that of Figure A1.6 (p S13) when also
# assuming that what's given as 12.18 in the text should actually be
# 121.8
prob_co_given_aa <- 0.20
# proportion of each exposure-covariate combination
prop_aa0_co0 <- (1 - prob_aa) * (1 - (prob_co_given_aa / 2))
prop_aa0_co1 <- (1 - prob_aa) * (prob_co_given_aa / 2)
prop_aa1_co0 <- prob_aa * (1 - prob_co_given_aa)
prop_aa1_co1 <- prob_aa * prob_co_given_aa
# Weibull distribution parameters for death times ('sh' for shape,
# 'sc' for scale)
sh_t <- 1.25
sc_t_aa0_co0 <- 121.8
sc_t_aa0_co1 <- sc_t_aa0_co0 * 0.57
sc_t_aa1_co0 <- sc_t_aa0_co0 * 0.49
sc_t_aa1_co1 <- sc_t_aa0_co0 * 0.57 * 0.49
# exponential distribution parameters for censoring times
mu_c_aa0_co0 <- 36.7
mu_c_aa0_co1 <- mu_c_aa0_co0 * 0.50
mu_c_aa1_co0 <- mu_c_aa0_co0 * 0.58
mu_c_aa1_co1 <- mu_c_aa0_co0 * 0.58 * 0.50


# ** Graphical constants

label_visit <- "Visit"
label_surv <- "Survival"
label_sqerr_surv <- "Survival Mean Squared Error"
label_rd <- "Risk Difference (RD)"
label_sqerr_rd <- "RD Mean Squared Error"
label_lrr <- "Log Risk Ratio (RR)"
label_sqerr_lrr <- "Log RR Mean Squared Error"
label_lor <- "Log Odds Ratio (OR)"
label_sqerr_lor <- "Log OR Mean Squared Error"
breaks_visit <- seq(0, 36, 6)
breaks_surv <- seq(0, 1, 0.2)
breaks_sqerr <- c(
  5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 5
)
breaks_rd <- seq(-0.1, 0.5, 0.1)
breaks_lr <- c(-1, 0, 1, 2, 3)
labels_sqerr <- c(
  "5E-8", "5E-7", "5E-6", "5E-5", "5E-4", "5E-3", "5E-2", "5E-1", "5"
)
limits_surv <- c(0, 1)
limits_sqerr <- c(5e-8, 5)
limits_rd <- c(-0.1, 0.5)
limits_lr <- c(-1, 3)
# line type for crude, nonparametric-weighted, weighted, and true
# curves
linetypes <- c("solid", "dotted", "dashed", "solid")
line_alpha <- 0.7
# colors for crude, nonparametric-weighted, weighted, and true curves
colours <- c("gray", "#FF653F", "black", "black")
theme_base_size <- 18
theme_base_line_size_mult <- 1 / 16
legend_coord_topright <- c(.74, .88)
legend_coord_topleft <- c(.28, .905)
legend_text_size <- 13
legend_key_height <- unit(8, "pt")
x_axis_title_top_padding_pts <- 10
y_axis_title_right_padding_pts <- 10

wi_label_y <- "Proportion at risk with L = 1"
wi_linetypes <- rep(c("solid", "dashed"), times = 3)
wi_colours <- rep(c("#760031", "#121358", "#8C5A3C"), each = 2)
wi_legend_coord_bottomleft <- c(.11, .18)


# ** Simulation scenarios

# *** Figure 1 / Figure A1.1-A1.5

truth_fig1 <- tibble(
  time = 1:36,
  surv_pooled = (
    0.25 * exp(-time / mu_t)
    + 0.25 * exp(-time / (mu_t * 0.3))
    + 0.25 * exp(-time / (mu_t * 0.41))
    + 0.25 * exp(-time / (mu_t * 0.3 * 0.41))
  ),
  surv_idu0 = (
    0.5 * exp(-time / mu_t)
    + 0.5 * exp(-time / (mu_t * 0.41))
  ),
  surv_idu1 = (
    0.5 * exp(-time / (mu_t * 0.3))
    + 0.5 * exp(-time / (mu_t * 0.3 * 0.41))
  ),
  rd = (1 - surv_idu1) - (1 - surv_idu0),
  lrr = log((1 - surv_idu1) / (1 - surv_idu0)),
  lor = log(
    ((1 - surv_idu1) / surv_idu1)
    / ((1 - surv_idu0) / surv_idu0)
  )
)

# *** Figure 2 / Figure A1.6

truth_fig2 <- tibble(
  time = 1:36,
  surv_aa0_co0 = exp(-(time / sc_t_aa0_co0)^sh_t),
  surv_aa0_co1 = exp(-(time / sc_t_aa0_co1)^sh_t),
  surv_aa1_co0 = exp(-(time / sc_t_aa1_co0)^sh_t),
  surv_aa1_co1 = exp(-(time / sc_t_aa1_co1)^sh_t),
  surv_pooled = (
    prop_aa0_co0 * surv_aa0_co0
    + prop_aa0_co1 * surv_aa0_co1
    + prop_aa1_co0 * surv_aa1_co0
    + prop_aa1_co1 * surv_aa1_co1
  ),
  surv_aa0 = (
    prop_aa0_co0 / (prop_aa0_co0 + prop_aa0_co1) * surv_aa0_co0
    + prop_aa0_co1 / (prop_aa0_co0 + prop_aa0_co1) * surv_aa0_co1
  ),
  surv_aa1 = (
    prop_aa1_co0 / (prop_aa1_co0 + prop_aa1_co1) * surv_aa1_co0
    + prop_aa1_co1 / (prop_aa1_co0 + prop_aa1_co1) * surv_aa1_co1
  ),
  rd = (1 - surv_aa1) - (1 - surv_aa0),
  lrr = log((1 - surv_aa1) / (1 - surv_aa0)),
  lor = log(
    ((1 - surv_aa1) / surv_aa1)
    / ((1 - surv_aa0) / surv_aa0)
  )
)


# * Function definitions


# ** Helper functions

tt <- function() {
  format(Sys.time(), format = "%Y%m%dT%H%M%S")
}

RD <- function(p0, p1) {
  p1 - p0
}

RR <- function(p0, p1) {
  dplyr::if_else(p0 > 0, p1 / p0, NA_real_)
}

OR <- function(p0, p1) {
  dplyr::if_else(
    p1 < 1 & p0 > 0 & p0 < 1,
    (p1 / (1 - p1)) / (p0 / (1 - p0)),
    NA_real_
  )
}


# ** Data generating processes

# *** Figure 1, Diagram I / Figure A1.1

make_data_f1d1 <- function(n = 1000) {
  tibble::tibble(
    id = seq_len(n),
    a = rbinom(n, 1, prob_idu),
    l = rbinom(n, 1, prob_hau),
    death_time_cont = dplyr::case_when(
      a == 1 & l == 1 ~ rexp(n, 1 / (mu_t * 0.3 * 0.41)),
      a == 1 & l == 0 ~ rexp(n, 1 / (mu_t * 0.3)),
      a == 0 & l == 1 ~ rexp(n, 1 / (mu_t * 0.41)),
      a == 0 & l == 0 ~ rexp(n, 1 / mu_t),
      .default = NA_real_
    ),
    censor_time_cont = rexp(n, 1 / mu_c),
    death_time = as.integer(ceiling(death_time_cont)),
    censor_time = as.integer(ceiling(censor_time_cont)),
    study_end = rep(36L, n),
    observed_time = pmin(death_time, censor_time, study_end),
    lost = ifelse(censor_time_cont < observed_time, 1L, 0L),
    died = ifelse(death_time_cont < observed_time, 1L, 0L)
  ) |>
    dplyr::select(!dplyr::ends_with("_cont"))
}


# *** Figure 1, Diagram II / Figure A1.2

make_data_f1d2 <- function(n = 1000) {
  tibble::tibble(
    id = seq_len(n),
    a = rbinom(n, 1, prob_idu),
    l = rbinom(n, 1, prob_hau),
    death_time_cont = dplyr::case_when(
      a == 1 & l == 1 ~ rexp(n, 1 / (mu_t * 0.3 * 0.41)),
      a == 1 & l == 0 ~ rexp(n, 1 / (mu_t * 0.3)),
      a == 0 & l == 1 ~ rexp(n, 1 / (mu_t * 0.41)),
      a == 0 & l == 0 ~ rexp(n, 1 / mu_t),
      .default = NA_real_
    ),
    censor_time_cont = ifelse(
      l == 1,
      rexp(n, 1 / (mu_c * 0.27)),
      rexp(n, 1 / mu_c)
    ),
    death_time = as.integer(ceiling(death_time_cont)),
    censor_time = as.integer(ceiling(censor_time_cont)),
    study_end = rep(36L, n),
    observed_time = pmin(death_time, censor_time, study_end),
    lost = ifelse(censor_time_cont < observed_time, 1L, 0L),
    died = ifelse(death_time_cont < observed_time, 1L, 0L)
  ) |>
    dplyr::select(!dplyr::ends_with("_cont"))
}

# *** Figure 1, Diagram III / Figure A1.3

make_data_f1d3 <- function(n = 1000) {
  tibble::tibble(
    id = seq_len(n),
    a = rbinom(n, 1, prob_idu),
    l = rbinom(n, 1, prob_hau),
    death_time_cont = dplyr::case_when(
      a == 1 & l == 1 ~ rexp(n, 1 / (mu_t * 0.3 * 0.41)),
      a == 1 & l == 0 ~ rexp(n, 1 / (mu_t * 0.3)),
      a == 0 & l == 1 ~ rexp(n, 1 / (mu_t * 0.41)),
      a == 0 & l == 0 ~ rexp(n, 1 / mu_t),
      .default = NA_real_
    ),
    censor_time_cont = ifelse(
      a == 1,
      rexp(n, 1 / (mu_c * 0.27)),
      rexp(n, 1 / mu_c)
    ),
    death_time = as.integer(ceiling(death_time_cont)),
    censor_time = as.integer(ceiling(censor_time_cont)),
    study_end = rep(36L, n),
    observed_time = pmin(death_time, censor_time, study_end),
    lost = ifelse(censor_time_cont < observed_time, 1L, 0L),
    died = ifelse(death_time_cont < observed_time, 1L, 0L)
  ) |>
    dplyr::select(!dplyr::ends_with("_cont"))
}

# *** Figure 1, Diagram IV / Figure A1.4

make_data_f1d4 <- function(n = 1000) {
  tibble::tibble(
    id = seq_len(n),
    a = rbinom(n, 1, prob_idu),
    l = rbinom(n, 1, prob_hau),
    death_time_cont = dplyr::case_when(
      a == 1 & l == 1 ~ rexp(n, 1 / (mu_t * 0.3 * 0.41)),
      a == 1 & l == 0 ~ rexp(n, 1 / (mu_t * 0.3)),
      a == 0 & l == 1 ~ rexp(n, 1 / (mu_t * 0.41)),
      a == 0 & l == 0 ~ rexp(n, 1 / mu_t),
      .default = NA_real_
    ),
    censor_time_cont = dplyr::case_when(
      a == 1 & l == 1 ~ rexp(n, 1 / (mu_c * 0.27 * 0.41)),
      a == 1 & l == 0 ~ rexp(n, 1 / (mu_c * 0.27)),
      a == 0 & l == 1 ~ rexp(n, 1 / (mu_c * 0.41)),
      a == 0 & l == 0 ~ rexp(n, 1 / mu_c),
      .default = NA_real_
    ),
    death_time = as.integer(ceiling(death_time_cont)),
    censor_time = as.integer(ceiling(censor_time_cont)),
    study_end = rep(36L, n),
    observed_time = pmin(death_time, censor_time, study_end),
    lost = ifelse(censor_time_cont < observed_time, 1L, 0L),
    died = ifelse(death_time_cont < observed_time, 1L, 0L)
  ) |>
    dplyr::select(!dplyr::ends_with("_cont"))
}

# *** Figure 1, Diagram V / Figure A1.5

make_data_f1d5 <- function(n = 1000) {
  tibble::tibble(
    id = seq_len(n),
    q = rbinom(n, 1, 0.50),
    a = rbinom(n, 1, ifelse(q == 1, 0.6, 0.4)),
    z = rbinom(n, 1, 0.50),
    l = dplyr::case_when(
      q == 1 & z == 1 ~ rbinom(n, 1, 0.10 * 3 * 3),
      q == 1 & z == 0 ~ rbinom(n, 1, 0.10 * 3),
      q == 0 & z == 1 ~ rbinom(n, 1, 0.10 * 3),
      q == 0 & z == 0 ~ rbinom(n, 1, 0.10)
    ),
    death_time_cont = dplyr::case_when(
      a == 1 & z == 1 ~ rexp(n, 1 / (mu_t * 0.27 * 0.27)),
      a == 1 & z == 0 ~ rexp(n, 1 / (mu_t * 0.27)),
      a == 0 & z == 1 ~ rexp(n, 1 / (mu_t * 0.27)),
      a == 0 & z == 0 ~ rexp(n, 1 / mu_t),
      .default = NA_real_
    ),
    censor_time_cont = dplyr::case_when(
      a == 1 & l == 1 ~ rexp(n, 1 / (mu_c * 0.37 * 0.37)),
      a == 1 & l == 0 ~ rexp(n, 1 / (mu_c * 0.37)),
      a == 0 & l == 1 ~ rexp(n, 1 / (mu_c * 0.37)),
      a == 0 & l == 0 ~ rexp(n, 1 / mu_c),
      .default = NA_real_
    ),
    death_time = as.integer(ceiling(death_time_cont)),
    censor_time = as.integer(ceiling(censor_time_cont)),
    study_end = rep(36L, n),
    observed_time = pmin(death_time, censor_time, study_end),
    lost = ifelse(censor_time_cont < observed_time, 1L, 0L),
    died = ifelse(death_time_cont < observed_time, 1L, 0L)
  ) |>
    dplyr::select(!dplyr::ends_with("_cont"))
}

# *** Figure 2, Diagram I / Figure A1.6

make_data_f2d1 <- function(n = 1000) {
  tibble::tibble(
    id = seq_len(n),
    a = rbinom(n, 1, prob_aa),
    l = rbinom(
      n, 1,
      prob = ifelse(a == 1L, prob_co_given_aa, 0.5 * prob_co_given_aa)
    ),
    death_time_cont = dplyr::case_when(
      a == 0 & l == 0 ~ rweibull(n, sh_t, sc_t_aa0_co0),
      a == 0 & l == 1 ~ rweibull(n, sh_t, sc_t_aa0_co1),
      a == 1 & l == 0 ~ rweibull(n, sh_t, sc_t_aa1_co0),
      a == 1 & l == 1 ~ rweibull(n, sh_t, sc_t_aa1_co1),
      .default = NA_real_
    ),
    censor_time_cont = dplyr::case_when(
      a == 0 & l == 0 ~ rexp(n, 1 / mu_c_aa0_co0),
      a == 0 & l == 1 ~ rexp(n, 1 / mu_c_aa0_co1),
      a == 1 & l == 0 ~ rexp(n, 1 / mu_c_aa1_co0),
      a == 1 & l == 1 ~ rexp(n, 1 / mu_c_aa1_co1),
      .default = NA_real_
    ),
    death_time = as.integer(ceiling(death_time_cont)),
    censor_time = as.integer(ceiling(censor_time_cont)),
    study_end = rep(36L, n),
    observed_time = pmin(death_time, censor_time, study_end),
    lost = ifelse(censor_time_cont < observed_time, 1L, 0L),
    died = ifelse(death_time_cont < observed_time, 1L, 0L)
  ) |>
    dplyr::select(!dplyr::ends_with("_cont"))
}


# ** Discrete survival analysis functions

# *** lambda_hat, estimate the discrete hazard function

# Estimate the hazard of `event` at each time point in the data,
# optionally stratified by `stratum` and incorporating `weights`.  All
# arguments following `data` are data-masked, meaning they can be
# referred to by name as in `base::subset()` or `dplyr::mutate()`.
lambda_hat <- function(
  data, time, time2 = NULL, event, stratum = NULL, weights = NULL,
  var = time_
) {
  # If `data` is empty, return a zero-row output data set with the
  # correct columns.  (This is needed to handle the toy data set of p
  # S25.)
  if (nrow(data) == 0 && rlang::quo_is_null(rlang::enquo(stratum))) {
    return(
      tibble::tibble(
        {{var}} := numeric(0),
        lambda_pooled = numeric(0)
      )
    )
  } else if (nrow(data) == 0 && !rlang::quo_is_null(rlang::enquo(stratum))) {
    return(
      tibble::tibble(
        {{var}} := numeric(0),
        "lambda_{{stratum}}0" := numeric(0),
        "lambda_{{stratum}}1" := numeric(0)
      )
    )
  }
  # To avoid having to repeatedly check if `weights` is `NULL`, assign
  # a value if it wasn't specified.  Because `weights` is a
  # data-masked argument, we need to create it as a new quosure,
  # adding it to the same environment as the other data-masked
  # arguments.  Compare `print(rlang::enquo(weights))` before and
  # after this block when `weights` is `NULL`.  See [3] for the origin
  # of the condition for checking for a null data-masked argument and
  # [4] and [5] for further details on data-masking.
  if (rlang::quo_is_null(rlang::enquo(weights))) {
    data$weights <- 1
    weights <- rlang::new_quosure(
      rlang::expr(weights),
      env = rlang::env_parent()
    )
  }
  if (rlang::quo_is_null(rlang::enquo(time2))) {
    # data$time2 <- dplyr::pull(data, {{time}})
    data <- data |>
      dplyr::mutate(
        time2 = {{time}},
        {{time}} := 1L
      )
    time2 <- rlang::new_quosure(
      rlang::expr(time2),
      env = rlang::env_parent()
    )
  }
  # Get the sequence of times at which to estimate the hazard.
  times <- seq_len(max(pull(data, {{time2}})))
  # Compute the hazard estimate at time t by dividing the sum of the
  # sample weights among observations with the event at time t by the
  # sum of the sample weights among observations at risk for the event
  # at time t.  For each t, the result is a one-row tibble that gets
  # added to a list.
  result_list <- vector(mode = "list", length = length(times))
  if (rlang::quo_is_null(rlang::enquo(stratum))) {
    # If `stratum` is `NULL`, estimate the hazard in the complete
    # sample, creating the variable `lambda_pooled`.
    for (t in times) {
      # Get the weights of records with the event at time t.
      weights_event_at_t <- data |>
        dplyr::filter({{time2}} == t & {{event}} == 1) |>
        dplyr::pull({{weights}})
      # Get the weights of records at risk of the event at time t.  To
      # avoid counting the same individual multiple times in a data
      # set organized into intervals, select only the interval
      # overlapping t.  `dplyr::between()` would be the natural choice
      # here, but it requires that its first argument be the same
      # length as the boundaries or for the boundaries to be
      # recyclable to the length of the first argument; since `t` is a
      # scalar, this won't work and leads to an error.  Logical
      # expressions will do the trick.
      weights_at_risk_at_t <- data |>
        dplyr::filter({{time}} <= t & t <= {{time2}}) |>
        dplyr::pull({{weights}})
      result_list[[t]] <- tibble(
        {{var}} := t,
        lambda_pooled =
          sum(weights_event_at_t) / sum(weights_at_risk_at_t)
      )
    }
  } else {
    # If `stratum` is not `NULL`, compute the hazard estimate at time
    # t in each stratum separately, creating the variables
    # `lambda_{stratum}0` and `lambda_{stratum}1`.
    for (t in times) {
      weights_event_at_t_strat0 <- data |>
        dplyr::filter(
          {{stratum}} == 0 & {{time2}} == t & {{event}} == 1
        ) |>
        dplyr::pull({{weights}})
      weights_event_at_t_strat1 <- data |>
        dplyr::filter(
          {{stratum}} == 1 & {{time2}} == t & {{event}} == 1
        ) |>
        dplyr::pull({{weights}})
      weights_at_risk_at_t_strat0 <- data |>
        dplyr::filter(
          {{stratum}} == 0 & {{time}} <= t & t <= {{time2}}
        ) |>
        dplyr::pull({{weights}})
      weights_at_risk_at_t_strat1 <- data |>
        dplyr::filter(
          {{stratum}} == 1 & {{time}} <= t & t <= {{time2}}
        ) |>
        dplyr::pull({{weights}})
      result_list[[t]] <- tibble::tibble(
        {{var}} := t,
        "lambda_{{stratum}}0" :=
          sum(weights_event_at_t_strat0) / sum(weights_at_risk_at_t_strat0),
        "lambda_{{stratum}}1" :=
          sum(weights_event_at_t_strat1) / sum(weights_at_risk_at_t_strat1)
      )
    }
  }
  # Bind the one-row data sets in `result_list` into a tibble.
  dplyr::bind_rows(result_list)
}


# *** S_hat, estimate survival from the discrete hazard function

S_hat <- function(data, stratum = NULL) {
  if (rlang::quo_is_null(rlang::enquo(stratum))) {
    data |>
      dplyr::mutate(
        surv_pooled = cumprod(1 - lambda_pooled),
        .keep = "unused"
      )
  } else {
    stratum_str <- rlang::as_name(rlang::enquo(stratum))
    lambda_strat0_sym <- rlang::sym(paste0("lambda_", stratum_str, "0"))
    lambda_strat1_sym <- rlang::sym(paste0("lambda_", stratum_str, "1"))
    data |>
      dplyr::mutate(
        "surv_{{stratum}}0" := cumprod(1 - !!lambda_strat0_sym),
        "surv_{{stratum}}1" := cumprod(1 - !!lambda_strat1_sym),
        .keep = "unused"
      )
  }
}


# *** estimate_effects, compute the RD, RR, and OR

estimate_effects <- function(data, stratum) {
  stratum_str <- rlang::as_name(rlang::enquo(stratum))
  surv_strat0_sym <- rlang::sym(paste0("surv_", stratum_str, "0"))
  surv_strat1_sym <- rlang::sym(paste0("surv_", stratum_str, "1"))
  risk_strat0_sym <- rlang::sym(paste0("risk_", stratum_str, "0"))
  risk_strat1_sym <- rlang::sym(paste0("risk_", stratum_str, "1"))
  data |>
    dplyr::mutate(
      "risk_{{stratum}}0" := 1 - !!surv_strat0_sym,
      "risk_{{stratum}}1" := 1 - !!surv_strat1_sym,
      rd = RD(!!risk_strat0_sym, !!risk_strat1_sym),
      rr = RR(!!risk_strat0_sym, !!risk_strat1_sym),
      lrr = dplyr::if_else(!is.na(rr), log(rr), NA_real_),
      or = OR(!!risk_strat0_sym, !!risk_strat1_sym),
      lor = dplyr::if_else(!is.na(or), log(or), NA_real_)
    ) |>
    dplyr::select(!dplyr::starts_with(c("surv_", "risk_")), -rr, -or)
    # dplyr::select(!dplyr::starts_with(c("surv_")))
}


# ** Inverse probability weighting functions

# *** expand_singly, expand a data set to one row per time point

expand_singly <- function(
  data,
  time,
  time2 = NULL,
  invariants = NULL,
  events = NULL,
  var = time_
) {
  # If `invariants` has a value, we assume that columns other than
  # `{{time}}`, `{{time2}}` and those columns in `invariants` and
  # `events` should not be included in the output.  Otherwise, all
  # columns other than `{{time}}`, `{{time2}}` and those columns in
  # `events` are treated as invariants.
  if (!rlang::quo_is_null(rlang::enquo(invariants))) {
    data <- data |>
      dplyr::select(
        c(
          {{time}},
          {{time2}},
          !!!rlang::enquos(invariants),
          !!!rlang::enquos(events)
        )
      )
  }
  data |>
    dplyr::rowwise() |>
    dplyr::mutate(
      {{var}} := 
        if (rlang::quo_is_null(rlang::enquo(time2))) {
          list(seq(1L, {{time}}, by = 1))
        } else {
          list(seq({{time}}, {{time2}}, by = 1))
        },
      dplyr::across(
        !!!rlang::enquos(events),
        \(x) list(c(rep(0L, times = length({{var}}) - 1), x))
      )
    ) |>
    tidyr::unnest(
      c({{var}}, !!!rlang::enquos(invariants), !!!rlang::enquos(events))
    )
}

# *** collapse_intervals, collapse a data set of individual time increments

collapse_intervals <- function(
  data, interval, events = NULL, invariants = NULL
) {
  data |>
    dplyr::summarise(
      "{{interval}}_start" := min({{interval}}),
      "{{interval}}_end" := max({{interval}}),
      dplyr::across(!!!rlang::enquos(events), max),
      .by = c(!!!rlang::enquos(invariants))
    )
}

# *** estimate_np_ipcw, estimate nonparametric (un)stabilized inverse probability of censoring weights

estimate_np_ipcw <- function(
  data, time, censored, event, exposure, stratum, id
) {
  data_marg_0 <- data |> dplyr::filter({{exposure}} == 0)
  data_marg_1 <- data |> dplyr::filter({{exposure}} == 1)
  data_cond_00 <- data |>
    dplyr::filter({{exposure}} == 0 & {{stratum}} == 0)
  data_cond_01 <- data |>
    dplyr::filter({{exposure}} == 0 & {{stratum}} == 1)
  data_cond_10 <- data |>
    dplyr::filter({{exposure}} == 1 & {{stratum}} == 0)
  data_cond_11 <- data |>
    dplyr::filter({{exposure}} == 1 & {{stratum}} == 1)
  lambda_marg_0 <- data_marg_0 |>
    lambda_hat(time = {{time}}, event = {{censored}}, var = {{time}}) |>
    dplyr::mutate(
      {{exposure}} := 0,
      marg_prob_nltfu = 1 - lambda_pooled
    )
  lambda_marg_1 <- data_marg_1 |>
    lambda_hat(time = {{time}}, event = {{censored}}, var = {{time}}) |>
    dplyr::mutate(
      {{exposure}} := 1,
      marg_prob_nltfu = 1 - lambda_pooled
    )
  lambda_cond_00 <- data_cond_00 |>
    lambda_hat(time = {{time}}, event = {{censored}}, var = {{time}}) |>
    dplyr::mutate(
      {{exposure}} := 0,
      {{stratum}} := 0,
      cond_prob_nltfu = 1 - lambda_pooled
    )
  lambda_cond_01 <- data_cond_01 |>
    lambda_hat(time = {{time}}, event = {{censored}}, var = {{time}}) |>
    dplyr::mutate(
      {{exposure}} := 0,
      {{stratum}} := 1,
      cond_prob_nltfu = 1 - lambda_pooled
    )
  lambda_cond_10 <- data_cond_10 |>
    lambda_hat(time = {{time}}, event = {{censored}}, var = {{time}}) |>
    dplyr::mutate(
      {{exposure}} := 1,
      {{stratum}} := 0,
      cond_prob_nltfu = 1 - lambda_pooled
    )
  lambda_cond_11 <- data_cond_11 |>
    lambda_hat(time = {{time}}, event = {{censored}}, var = {{time}}) |>
    dplyr::mutate(
      {{exposure}} := 1,
      {{stratum}} := 1,
      cond_prob_nltfu = 1 - lambda_pooled
    )
  lambda_marg <- dplyr::bind_rows(lambda_marg_0, lambda_marg_1)
  lambda_cond <- dplyr::bind_rows(
    lambda_cond_00, lambda_cond_01, lambda_cond_10, lambda_cond_11
  )
  data |>
    dplyr::rename(time_ := {{time}}) |>
    expand_singly(
      time = time_,
      invariants = c({{id}}, {{exposure}}, {{stratum}}),
      events = c({{censored}}, {{event}}),
      var = {{time}}
    ) |>
    dplyr::inner_join(
      lambda_marg, dplyr::join_by({{time}}, {{exposure}})
    ) |>
    dplyr::inner_join(
      lambda_cond, dplyr::join_by({{time}}, {{exposure}}, {{stratum}})
    ) |>
    dplyr::mutate(
      wt = dplyr::if_else(
        {{event}} == 1, 0, cumprod(1 / cond_prob_nltfu)
      ),
      swt = cumprod(marg_prob_nltfu / cond_prob_nltfu),
      .by = {{id}}
    ) |>
    collapse_intervals(
      interval = {{time}},
      events = c({{censored}}, {{event}}),
      invariants = c({{id}}, {{exposure}}, {{stratum}}, wt, swt)
    )
}

# *** estimate_ipcw, estimate parametric (un)stabilized inverse probability of censoring weights

estimate_ipcw <- function(data, formula_numer, formula_denom, id, censored) {
  model_numer <- glm(formula_numer, family = binomial, data = data)
  model_denom <- glm(formula_denom, family = binomial, data = data)
  data |>
    dplyr::mutate(
      pred_marg = predict(model_numer, type = "response"),
      pred_cond = predict(model_denom, type = "response"),
      # Per Howe et al. (2016, p S15), unstabilized weights are zero
      # at and after loss to follow up.
      wt_ = dplyr::if_else({{censored}} == 1, 0, 1 / pred_cond),
      swt_ = pred_marg / pred_cond
    ) |>
    dplyr::mutate(
      wt = round(cumprod(wt_), 4),
      swt = round(cumprod(swt_), 4),
      .by = {{id}},
      .keep = "unused"
    )
}


# ** Simulation runners and summarizers

# *** sim_1, run one simulation

sim_1 <- function(data_generator, n = n) {
  fake_data <- data_generator(n)
  # Run a naive analysis.
  naive_pooled_survival <- fake_data |>
    lambda_hat(observed_time, event = died, var = observed_time) |>
    S_hat()
  naive_stratified_survival <- fake_data |>
    lambda_hat(
      observed_time, event = died, stratum = a, var = observed_time
    ) |>
    S_hat(a) |>
    estimate_effects(a)
  naive_estimates <- naive_pooled_survival |>
    dplyr::inner_join(naive_stratified_survival, by = "observed_time") |>
    dplyr::rename_with(
      \(x) ifelse(x == "observed_time", x, paste0(x, "_naive"))
    )
  # Run a nonparametric IPC-weighted analysis.
  np_fake_data <- fake_data |>
    estimate_np_ipcw(
      time = observed_time,
      censored = lost,
      event = died,
      exposure = a,
      stratum = l,
      id = id
    )
  np_pooled_survival <- np_fake_data |>
    lambda_hat(
      time = observed_time_start,
      time2 = observed_time_end,
      event = died,
      weights = swt,
      var = observed_time
    ) |>
    S_hat()
  np_stratified_survival <- np_fake_data |>
    lambda_hat(
      time = observed_time_start,
      time2 = observed_time_end,
      event = died,
      stratum = a,
      weights = swt,
      var = observed_time
    ) |>
    S_hat(a) |>
    estimate_effects(a)
  np_estimates <- np_pooled_survival |>
    dplyr::inner_join(
      np_stratified_survival,
      by = "observed_time"
    ) |>
    dplyr::rename_with(
      \(x) ifelse(x == "observed_time", x, paste0(x, "_np"))
    )
  # Run an IPC-weighted analysis.
  weighted_fake_data <- fake_data |>
    expand_singly(
      observed_time,
      invariants = c(id, a, l),
      events = c(lost, died),
      var = observed_time
    ) |>
    estimate_ipcw(
      formula_numer = lost == 0 ~ 0 + factor(observed_time) * a,
      formula_denom = lost == 0 ~ 0 + factor(observed_time) * a * l,
      id = id,
      censored = lost
    ) |>
    collapse_intervals(
      interval = observed_time,
      events = c(lost, died),
      invariants = c(id, a, l, wt, swt)
    )
  ipcw_pooled_survival <- weighted_fake_data |>
    lambda_hat(
      time = observed_time_start,
      time2 = observed_time_end,
      event = died,
      weights = swt,
      var = observed_time
    ) |>
    S_hat()
  ipcw_stratified_survival <- weighted_fake_data |>
    lambda_hat(
      time = observed_time_start,
      time2 = observed_time_end,
      event = died,
      stratum = a,
      weights = swt,
      var = observed_time
    ) |>
    S_hat(a) |>
    estimate_effects(a)
  ipcw_estimates <- ipcw_pooled_survival |>
    dplyr::inner_join(
      ipcw_stratified_survival,
      by = "observed_time"
    ) |>
    dplyr::rename_with(
      \(x) ifelse(x == "observed_time", x, paste0(x, "_ipcw"))
    )
  purrr::reduce(
    list(naive_estimates, np_estimates, ipcw_estimates),
    \(x, y)
      dplyr::inner_join(x, y, by = dplyr::join_by(observed_time))
  )
}

# *** sim, run `n_sim` simulations

sim <- function(data_generator, n_sim = m, n = n, seed = FALSE) {
  p <- progressr::progressor(steps = n_sim)
  sim_list <- furrr::future_map(
    seq_len(n_sim),
    \(i) {
      p()
      sim_1(data_generator, n = n) |>
        dplyr::mutate(sim_id = i, .before = 1)
    },
    .options = furrr::furrr_options(seed = seed)
  )
  dplyr::bind_rows(sim_list)
}

# *** summarize_sims, summarize simulation results

# Summarize variables in a set of simulation results by another
# factor, e.g., visit number.  For each column in `...`, calculate the
# mean and variance (required for the subsequent calculation of
# simulation mean squared error) at each level of `.by` across
# simulations.  The resulting tibble will contain `.by` along with a
# pair of columns for each variable in `...`---one for the mean and
# one for the variance.  Because dplyr verbs like `summarise()` allow
# later columns to use results of computations on earlier columns and
# because the mean columns reuse the variable names, we need to
# calculate the variances first, otherwise we'll be taking the
# variance of a single number, which returns `NA`.
summarize_sims <- function(sim_results, ..., .by) {
  dots <- rlang::enquos(...)
  if (!is.data.frame(sim_results))
    sim_results <- dplyr::bind_rows(sim_results)
  dots |>
    purrr::map(
      \(col) {
        sim_results |>
          dplyr::summarise(
            "var_{{col}}" := var({{col}}, na.rm = TRUE),
            {{col}} := mean({{col}}, na.rm = TRUE),
            .by = {{.by}}
          )
      }
    )|>
    purrr::reduce(
      \(x, y) {
        dplyr::inner_join(x, y, by = rlang::as_label(rlang::enquo(.by)))
      }
    )
}


# *** compute_sim_stats, compute simulation statistics

compute_sim_stats <- function(data, truth) {
  truth |>
    dplyr::select(1, surv_pooled, rd, lrr, lor) |>
    # Rename the column `time` in `truth` to match the name of the
    # corresponding column in `data`.
    dplyr::rename_with(
      .fn = \(x) ifelse(x == "time", names(data)[[1]], paste0(x, "_true"))
    ) |>
    dplyr::inner_join(
      data, by = names(data)[[1]]
    ) |>
    dplyr::mutate(
      sqerr_surv_pooled_naive =
        (surv_pooled_true - surv_pooled_naive)^2 + var_surv_pooled_naive,
      sqerr_rd_naive = (rd_true - rd_naive)^2 + var_rd_naive,
      sqerr_lrr_naive = (lrr_true - lrr_naive)^2 + var_lrr_naive,
      sqerr_lor_naive = (lor_true - lor_naive)^2 + var_lor_naive,
      sqerr_surv_pooled_np =
        (surv_pooled_true - surv_pooled_np)^2 + var_surv_pooled_np,
      sqerr_rd_np = (rd_true - rd_np)^2 + var_rd_np,
      sqerr_lrr_np = (lrr_true - lrr_np)^2 + var_lrr_np,
      sqerr_lor_np = (lor_true - lor_np)^2 + var_lor_np,
      sqerr_surv_pooled_ipcw =
        (surv_pooled_true - surv_pooled_ipcw)^2 + var_surv_pooled_ipcw,
      sqerr_rd_ipcw = (rd_true - rd_ipcw)^2 + var_rd_ipcw,
      sqerr_lrr_ipcw = (lrr_true - lrr_ipcw)^2 + var_lrr_ipcw,
      sqerr_lor_ipcw = (lor_true - lor_ipcw)^2 + var_lor_ipcw
    ) |>
    dplyr::relocate(
      observed_time,
      surv_pooled_true,
      surv_pooled_naive, sqerr_surv_pooled_naive,
      surv_pooled_np, sqerr_surv_pooled_np,
      surv_pooled_ipcw, sqerr_surv_pooled_ipcw,
      rd_true,
      rd_naive, sqerr_rd_naive,
      rd_np, sqerr_rd_np,
      rd_ipcw, sqerr_rd_ipcw,
      lrr_true,
      lrr_naive, sqerr_lrr_naive,
      lrr_np, sqerr_lrr_np,
      lrr_ipcw, sqerr_lrr_ipcw,
      lor_true,
      lor_naive, sqerr_lor_naive,
      lor_np, sqerr_lor_np,
      lor_ipcw, sqerr_lor_ipcw
    )
}


# *** elongate_sim_stats, elongate simulation statistics to facilitate plotting

elongate_sim_stats <- function(sim_stats) {
  sim_stats |>
    tidyr::pivot_longer(
      cols = !observed_time,
      names_to = c(".value", "type"),
      names_pattern = "(.+)_(true|naive|np|ipcw)"
    ) |>
    dplyr::mutate(
      type = factor(type,
        levels = c("naive", "np", "ipcw", "true"),
        labels = c("Crude", "NP-weighted", "Weighted", "True")
      )
    )
}


# ** Simulation visualizers

# *** visualize_sim_stats, produce plots from simulation results

visualize_sim_stats_ <- function(elongated_sim_stats) {
  surv_plot <- elongated_sim_stats |>
    dplyr::filter_out(type == "NP-weighted") |>
    ggplot() +
    geom_line(
      aes(observed_time, surv_pooled, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_surv) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_surv, limits = limits_surv) +
    scale_linetype_manual(name = NULL, values = linetypes[-2]) +
    scale_colour_manual(name = NULL, values = colours[-2]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.77, .91),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_surv_plot <- elongated_sim_stats |>
    dplyr::filter(type %in% c("Crude", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_surv_pooled, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_surv) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[c(-2, -4)]) +
    scale_colour_manual(name = NULL, values = colours[c(-2, -4)]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.23, .93),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  rd_plot <- elongated_sim_stats |>
    dplyr::filter_out(type == "NP-weighted") |>
    ggplot() +
    geom_line(
      aes(observed_time, rd, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_rd) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_rd, limits = limits_rd) +
    scale_linetype_manual(name = NULL, values = linetypes[-2]) +
    scale_colour_manual(name = NULL, values = colours[-2]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.77, .91),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_rd_plot <- elongated_sim_stats |>
    filter(type %in% c("Crude", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_rd, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_rd) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[c(-2, -4)]) +
    scale_colour_manual(name = NULL, values = colours[c(-2, -4)]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.23, .93),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  lrr_plot <- elongated_sim_stats |>
    dplyr::filter_out(type == "NP-weighted") |>
    ggplot() +
    geom_line(
      aes(observed_time, lrr, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_lrr) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_lr, limits = limits_lr) +
    scale_linetype_manual(name = NULL, values = linetypes[-2]) +
    scale_colour_manual(name = NULL, values = colours[-2]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.77, .91),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_lrr_plot <- elongated_sim_stats |>
    filter(type %in% c("Crude", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_lrr, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_lrr) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[c(-2, -4)]) +
    scale_colour_manual(name = NULL, values = colours[c(-2, -4)]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.23, .93),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  lor_plot <- elongated_sim_stats |>
    dplyr::filter_out(type == "NP-weighted") |>
    ggplot() +
    geom_line(
      aes(observed_time, lor, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_lor) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_lr, limits = limits_lr) +
    scale_linetype_manual(name = NULL, values = linetypes[-2]) +
    scale_colour_manual(name = NULL, values = colours[-2]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.77, .91),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_lor_plot <- elongated_sim_stats |>
    filter(type %in% c("Crude", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_lor, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_lor) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[c(-2, -4)]) +
    scale_colour_manual(name = NULL, values = colours[c(-2, -4)]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = c(.23, .93),
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  list(
    surv_plot = surv_plot,
    sqerr_surv_plot = sqerr_surv_plot,
    rd_plot = rd_plot,
    sqerr_rd_plot = sqerr_rd_plot,
    lrr_plot = lrr_plot,
    sqerr_lrr_plot = sqerr_lrr_plot,
    lor_plot = lor_plot,
    sqerr_lor_plot = sqerr_lor_plot
  )
}

visualize_sim_stats <- function(elongated_sim_stats) {
  surv_plot <- elongated_sim_stats |>
    ggplot() +
    geom_line(
      aes(observed_time, surv_pooled, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_surv) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_surv, limits = limits_surv) +
    scale_linetype_manual(name = NULL, values = linetypes) +
    scale_colour_manual(name = NULL, values = colours) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topright,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_surv_plot <- elongated_sim_stats |>
    dplyr::filter(type %in% c("Crude", "NP-weighted", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_surv_pooled, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_surv) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[-4]) +
    scale_colour_manual(name = NULL, values = colours[-4]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topleft,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  rd_plot <- elongated_sim_stats |>
    ggplot() +
    geom_line(
      aes(observed_time, rd, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_rd) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_rd, limits = limits_rd) +
    scale_linetype_manual(name = NULL, values = linetypes) +
    scale_colour_manual(name = NULL, values = colours) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topright,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_rd_plot <- elongated_sim_stats |>
    filter(type %in% c("Crude", "NP-weighted", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_rd, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_rd) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[-4]) +
    scale_colour_manual(name = NULL, values = colours[-4]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topleft,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  lrr_plot <- elongated_sim_stats |>
    ggplot() +
    geom_line(
      aes(observed_time, lrr, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_lrr) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_lr, limits = limits_lr) +
    scale_linetype_manual(name = NULL, values = linetypes) +
    scale_colour_manual(name = NULL, values = colours) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topright,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_lrr_plot <- elongated_sim_stats |>
    filter(type %in% c("Crude", "NP-weighted", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_lrr, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_lrr) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[-4]) +
    scale_colour_manual(name = NULL, values = colours[-4]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topleft,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  lor_plot <- elongated_sim_stats |>
    ggplot() +
    geom_line(
      aes(observed_time, lor, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_lor) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_continuous(breaks = breaks_lr, limits = limits_lr) +
    scale_linetype_manual(name = NULL, values = linetypes) +
    scale_colour_manual(name = NULL, values = colours) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topright,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  sqerr_lor_plot <- elongated_sim_stats |>
    filter(type %in% c("Crude", "NP-weighted", "Weighted")) |>
    ggplot() +
    geom_line(
      aes(observed_time, sqerr_lor, linetype = type, colour = type),
      alpha = line_alpha
    ) +
    xlab(label_visit) +
    ylab(label_sqerr_lor) +
    scale_x_continuous(breaks = breaks_visit) +
    scale_y_log10(
      breaks = breaks_sqerr, labels = labels_sqerr, limits = limits_sqerr
    ) +
    scale_linetype_manual(name = NULL, values = linetypes[-4]) +
    scale_colour_manual(name = NULL, values = colours[-4]) +
    theme_bw(
      base_size = theme_base_size,
      base_line_size = theme_base_size * theme_base_line_size_mult
    ) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(
        margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
      ),
      axis.title.y = element_text(
        margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
      ),
      legend.background = element_rect(fill = fill_alpha("white", 0)),
      legend.key.height = legend_key_height,
      legend.position = "inside",
      legend.position.inside = legend_coord_topleft,
      legend.text = element_text(size = legend_text_size),
      panel.grid = element_blank()
    )
  list(
    surv_plot = surv_plot,
    sqerr_surv_plot = sqerr_surv_plot,
    rd_plot = rd_plot,
    sqerr_rd_plot = sqerr_rd_plot,
    lrr_plot = lrr_plot,
    sqerr_lrr_plot = sqerr_lrr_plot,
    lor_plot = lor_plot,
    sqerr_lor_plot = sqerr_lor_plot
  )
}

# *** organize_plots, organize visualizations (via patchwork) to resemble Howe et al. Figure A1-A6

organize_plots <- function(plots_list) {
  with(plots_list, (
    (surv_plot | sqerr_surv_plot | rd_plot | sqerr_rd_plot)
    / (lrr_plot | sqerr_lrr_plot | lor_plot | sqerr_lor_plot)
  ))
}


# ** Weighting illustration

# *** proportion_l_end, frequency of L at each visit for data with one time variable

proportion_l_end <- function(data) {
  purrr::map_dbl(
    1:36,
    \(x) {
      data |>
        dplyr::filter(
          (observed_time == x & died == 0 & lost == 0)  # study end
          | observed_time > x
        ) |>
        dplyr::pull(l) |>
        mean()
    }
  )
}

# *** proportion_l_interval, frequency of L at each visit for data with start/end time variables

proportion_l_interval <- function(
  data,
  type = c("crude", "weighted", "unstabilized")
) {
  purrr::map_dbl(
    1:36,
    \(x) {
      data <- data |>
        dplyr::filter(
          observed_time_start <= x
          & (
            (observed_time_end == x & died == 0 & lost == 0)
            | observed_time_end > x)
          # observed_time_start <= x
        ) |>
        dplyr::mutate(wtl = wt * l, swtl = swt * l)
      if (type == "crude") {
        data |>
          dplyr::pull(l) |>
          mean()
      } else if (type == "weighted") {
        data |>
          dplyr::pull(swtl) |>
          mean()
      } else if (type == "unstabilized") {
        data |>
          dplyr::pull(wtl) |>
          mean()
      }
    }
  )
}


# * Data generating process summaries

set.seed(93792515)

f1d1_dat <- make_data_f1d1(1e6)
mean(f1d1_dat$lost)
#> [1] 0.322912

f1d2_dat <- make_data_f1d2(1e6)
mean(f1d2_dat$lost)
#> [1] 0.507858

f1d3_dat <- make_data_f1d3(1e6)
mean(f1d3_dat$lost)
#> [1] 0.503599

f1d4_dat <- make_data_f1d4(1e6)
mean(f1d4_dat$lost)
#> [1] 0.627226

f1d5_dat <- make_data_f1d5(1e6)
mod_hau <- logbin(l ~ q + z, data = f1d5_dat)
# summary(mod_hau)
exp(coef(mod_hau))
#> (Intercept)           q           z 
#>   0.1000437   2.9921532   3.0045772 
mean(f1d5_dat$lost)
#> [1] 0.532583


# * Weighting demonstration

set.seed(66500707)

f1d1_dat <- make_data_f1d1(20000)
f1d2_dat <- make_data_f1d2(20000)
f1d4_dat <- make_data_f1d4(20000)

f1d1_wdat <- f1d1_dat |>
  expand_singly(
    observed_time,
    invariants = c(id, a, l),
    events = c(lost, died),
    var = observed_time
  ) |>
  estimate_ipcw(
    formula_numer = lost == 0 ~ 0 + factor(observed_time) * a,
    formula_denom = lost == 0 ~ 0 + factor(observed_time) * a * l,
    id = id,
    censored = lost
  ) |>
  collapse_intervals(
    interval = observed_time,
    events = c(lost, died),
    invariants = c(id, a, l, wt, swt)
  )
f1d2_wdat <- f1d2_dat |>
  expand_singly(
    observed_time,
    invariants = c(id, a, l),
    events = c(lost, died),
    var = observed_time
  ) |>
  estimate_ipcw(
    formula_numer = lost == 0 ~ 0 + factor(observed_time) * a,
    formula_denom = lost == 0 ~ 0 + factor(observed_time) * a * l,
    id = id,
    censored = lost
  ) |>
  collapse_intervals(
    interval = observed_time,
    events = c(lost, died),
    invariants = c(id, a, l, wt, swt)
  )
f1d4_wdat <- f1d4_dat |>
  expand_singly(
    observed_time,
    invariants = c(id, a, l),
    events = c(lost, died),
    var = observed_time
  ) |>
  estimate_ipcw(
    formula_numer = lost == 0 ~ 0 + factor(observed_time) * a,
    formula_denom = lost == 0 ~ 0 + factor(observed_time) * a * l,
    id = id,
    censored = lost
  ) |>
  collapse_intervals(
    interval = observed_time,
    events = c(lost, died),
    invariants = c(id, a, l, wt, swt)
  )

proportion_l_end(f1d1_dat)
proportion_l_interval(f1d1_wdat, type = "crude")
proportion_l_interval(f1d1_wdat, type = "weighted")

proportion_l_end(f1d2_dat)
proportion_l_interval(f1d2_wdat, type = "crude")
proportion_l_interval(f1d2_wdat, type = "weighted")

proportion_l_interval(f1d4_wdat, type = "crude")
proportion_l_interval(f1d4_wdat, type = "weighted")

ldat <- bind_rows(
  tibble(
    type = "I, crude", visit = 1:36,
    prop_l = proportion_l_interval(f1d1_wdat, type = "crude")
  ),
  tibble(
    type = "I, weighted", visit = 1:36,
    prop_l = proportion_l_interval(f1d1_wdat, type = "weighted")
  ),
  tibble(
    type = "II, crude", visit = 1:36,
    prop_l = proportion_l_interval(f1d2_wdat, type = "crude")
  ),
  tibble(
    type = "II, weighted", visit = 1:36,
    prop_l = proportion_l_interval(f1d2_wdat, type = "weighted")
  ),
  tibble(
    type = "IV, crude", visit = 1:36,
    prop_l = proportion_l_interval(f1d4_wdat, type = "crude")
  ),
  tibble(
    type = "IV, weighted", visit = 1:36,
    prop_l = proportion_l_interval(f1d4_wdat, type = "weighted")
  )
) |>
  mutate(
    type = factor(type,
      levels = c(
        "I, crude", "I, weighted", "II, crude", "II, weighted",
        "IV, crude", "IV, weighted"
      )
    )
  )

fig_wt_illus <- ldat |>
  ggplot() +
  geom_line(
    aes(visit, prop_l, colour = type, linetype = type),
    alpha = line_alpha
  ) +
  xlab(label_visit) +
  ylab(wi_label_y) +
  scale_x_continuous(breaks = breaks_visit) +
  scale_linetype_manual(name = NULL, values = wi_linetypes) +
  scale_colour_manual(name = NULL, values = wi_colours) +
  theme_bw(
    base_size = theme_base_size,
    base_line_size = theme_base_size * theme_base_line_size_mult
  ) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title.x = element_text(
      margin = margin(t = x_axis_title_top_padding_pts, unit = "pt")
    ),
    axis.title.y = element_text(
      margin = margin(r = y_axis_title_right_padding_pts, unit = "pt")
    ),
    legend.background = element_rect(fill = fill_alpha("white", 0)),
    legend.key.height = legend_key_height,
    legend.position = "inside",
    legend.position.inside = wi_legend_coord_bottomleft,
    legend.text = element_text(size = legend_text_size),
    panel.grid = element_blank()
  )
webpea(
  "./img/figure_weight_illustration.webp", plot = fig_wt_illus,
  width = 600, height = 315, dpi = 72, units = "px"
)


# * Simulations


# ** Run and save simulations

# times based on runs on a desktop computer with the following specs:
# Processor	Intel(R) Core(TM) i7-14700   2.10 GHz
# Installed RAM	32.0 GB (31.7 GB usable)
# System type	64-bit operating system, x64-based processor

future::plan(multisession)
progressr::handlers(global = TRUE)

system.time({
sim_f1d1 <- sim(make_data_f1d1, n_sim = m, n = n, seed = 40010419)
})
#>    user  system elapsed 
#>   17.36    0.71 1258.41 
saveRDS(sim_f1d1, file = "./sim_f1d1_500_1000_20260509.rds")

system.time({
sim_f1d2 <- sim(make_data_f1d2, n_sim = m, n = n, seed = 77022729)
})
#>    user  system elapsed 
#>   12.98    0.56  973.11 
saveRDS(sim_f1d2, file = "./sim_f1d2_500_1000_20260509.rds")

system.time({
sim_f1d3 <- sim(make_data_f1d3, n_sim = m, n = n, seed = 28872964)
})
#>    user  system elapsed 
#>   13.39    0.48  981.32 
saveRDS(sim_f1d3, file = "./sim_f1d3_500_1000_20260509.rds")

system.time({
sim_f1d4 <- sim(make_data_f1d4, n_sim = m, n = n, seed = 66925650)
})
#>    user  system elapsed 
#>    9.91    0.46  795.35 
saveRDS(sim_f1d4, file = "./sim_f1d4_500_1000_20260509.rds")

system.time({
sim_f1d5 <- sim(make_data_f1d5, n_sim = m, n = n, seed = 11941628)
})
#>    user  system elapsed 
#>    9.86    0.53  825.72 
saveRDS(sim_f1d5, file = "./sim_f1d5_500_1000_20260509.rds")

system.time({
sim_f2d1 <- sim(make_data_f2d1, n_sim = m, n = n, seed = 41783356)
})
#>    user  system elapsed 
#>   13.42    0.63  778.05 
saveRDS(sim_f2d1, file = "./sim_f2d1_500_1000_20260509.rds")


system.time({
sim_f1d1x <- sim(make_data_f1d1, n_sim = m, n = 20000, seed = 75916690)
})
#>     user   system  elapsed 
#>    25.54    42.61 19705.93 
saveRDS(sim_f1d1x, file = "./sim_f1d1_500_20000_20260509.rds")

system.time({
sim_f1d2x <- sim(make_data_f1d2, n_sim = m, n = 20000, seed = 66267682)
})
#>     user   system  elapsed 
#>    31.86    40.90 14493.89 
saveRDS(sim_f1d2x, file = "./sim_f1d2_500_20000_20260509.rds")

system.time({
sim_f1d3x <- sim(make_data_f1d3, n_sim = m, n = 20000, seed = 39149152)
})
#>     user   system  elapsed 
#>    33.51    44.65 14942.38 
saveRDS(sim_f1d3x, file = "./sim_f1d3_500_20000_20260509.rds")

system.time({
sim_f1d4x <- sim(make_data_f1d4, n_sim = m, n = 20000, seed = 493020)
})
#>     user   system  elapsed 
#>    43.30    46.35 13663.97 
saveRDS(sim_f1d4x, file = "./sim_f1d4_500_20000_20260509.rds")

system.time({
sim_f1d5x <- sim(make_data_f1d5, n_sim = m, n = 20000, seed = 21236133)
})
#>     user   system  elapsed 
#>    40.03    41.61 13347.89 
saveRDS(sim_f1d5x, file = "./sim_f1d5_500_20000_20260509.rds")

system.time({
sim_f2d1x <- sim(make_data_f2d1, n_sim = m, n = 20000, seed = 92189932)
})
#>     user   system  elapsed 
#>    34.93    37.54 12300.17 
saveRDS(sim_f2d1x, file = "./sim_f2d1_500_20000_20260509.rds")


# ** Load simulations

sim_f1d1_1000 <- readRDS("./results/sim_f1d1_500_1000_20260509.rds")
sim_f1d2_1000 <- readRDS("./results/sim_f1d2_500_1000_20260509.rds")
sim_f1d3_1000 <- readRDS("./results/sim_f1d3_500_1000_20260509.rds")
sim_f1d4_1000 <- readRDS("./results/sim_f1d4_500_1000_20260509.rds")
sim_f1d5_1000 <- readRDS("./results/sim_f1d5_500_1000_20260509.rds")
sim_f2d1_1000 <- readRDS("./results/sim_f2d1_500_1000_20260509.rds")

sim_f1d1_20000 <- readRDS("./results/sim_f1d1_500_20000_20260509.rds")
sim_f1d2_20000 <- readRDS("./results/sim_f1d2_500_20000_20260509.rds")
sim_f1d3_20000 <- readRDS("./results/sim_f1d3_500_20000_20260509.rds")
sim_f1d4_20000 <- readRDS("./results/sim_f1d4_500_20000_20260509.rds")
sim_f1d5_20000 <- readRDS("./results/sim_f1d5_500_20000_20260509.rds")
sim_f2d1_20000 <- readRDS("./results/sim_f2d1_500_20000_20260509.rds")


# ** Summarize simulations

stats_f1d1_1000 <- sim_f1d1_1000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d2_1000 <- sim_f1d2_1000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d3_1000 <- sim_f1d3_1000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d4_1000 <- sim_f1d4_1000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d5_1000 <- sim_f1d5_1000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)

stats_f1d1_20000 <- sim_f1d1_20000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d2_20000 <- sim_f1d2_20000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d3_20000 <- sim_f1d3_20000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d4_20000 <- sim_f1d4_20000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)
stats_f1d5_20000 <- sim_f1d5_20000 |>
  summarize_sims(
    surv_pooled_naive, rd_naive, lrr_naive, lor_naive,
    surv_pooled_np, rd_np, lrr_np, lor_np,
    surv_pooled_ipcw, rd_ipcw, lrr_ipcw, lor_ipcw,
    .by = observed_time
  ) |>
  compute_sim_stats(truth_fig1)


# ** Visualize simulations

# *** N = 1000 simulations

toplot_f1d1_1000 <- stats_f1d1_1000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d1_1000 <- toplot_f1d1_1000 |>
  visualize_sim_stats()
fig_f1d1_1000 <- plots_f1d1_1000 |>
  organize_plots()
# fig_f1d1_1000
# ggsave(
#   paste0(
#     "FigureA_1_1000_",
#     tt(),
#     ".pdf"
#   ),
#   plot = fig_f1d1_1000, width = 11, height = 7, units = "in"
# )
webpea(
  # paste0("Figure1_4_1000_", tt(), ".webp"),
  "FigureA_1_1000.webp", plot = fig_f1d1_1000,
  width = 1200, height = 630, dpi = 72, units = "px"
)

plots_f1d1_1000x <- toplot_f1d1_1000 |>
  visualize_sim_stats_()
fig_f1d1_1000x <- plots_f1d1_1000x |>
  organize_plots()
webpea(
  "FigureA_1_1000_noNpWeights.webp", plot = fig_f1d1_1000x,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d2_1000 <- stats_f1d2_1000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d2_1000 <- toplot_f1d2_1000 |>
  visualize_sim_stats()
fig_f1d2_1000 <- plots_f1d2_1000 |>
  organize_plots()
webpea(
  "FigureA_2_1000.webp", plot = fig_f1d2_1000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d2_1000_ <- toplot_f1d2_1000 |>
  visualize_sim_stats_()
fig_f1d2_1000_ <- plots_f1d2_1000_ |>
  organize_plots()
webpea(
  "./img/FigureA_2_1000_noNpWeights.webp", plot = fig_f1d2_1000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d3_1000 <- stats_f1d3_1000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d3_1000 <- toplot_f1d3_1000 |>
  visualize_sim_stats()
fig_f1d3_1000 <- plots_f1d3_1000 |>
  organize_plots()
webpea(
  "./img/FigureA_3_1000.webp", plot = fig_f1d3_1000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d3_1000_ <- toplot_f1d3_1000 |>
  visualize_sim_stats_()
fig_f1d3_1000_ <- plots_f1d3_1000_ |>
  organize_plots()
webpea(
  "./img/FigureA_3_1000_noNpWeights.webp", plot = fig_f1d3_1000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d4_1000 <- stats_f1d4_1000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d4_1000 <- toplot_f1d4_1000 |>
  visualize_sim_stats()
fig_f1d4_1000 <- plots_f1d4_1000 |>
  organize_plots()
webpea(
  "./img/FigureA_4_1000.webp", plot = fig_f1d4_1000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d4_1000_ <- toplot_f1d4_1000 |>
  visualize_sim_stats_()
fig_f1d4_1000_ <- plots_f1d4_1000_ |>
  organize_plots()
webpea(
  "./img/FigureA_4_1000_noNpWeights.webp", plot = fig_f1d4_1000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d5_1000 <- stats_f1d5_1000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d5_1000 <- toplot_f1d5_1000 |>
  visualize_sim_stats()
fig_f1d5_1000 <- plots_f1d5_1000 |>
  organize_plots()
webpea(
  "./img/FigureA_5_1000.webp", plot = fig_f1d5_1000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d5_1000_ <- toplot_f1d5_1000 |>
  visualize_sim_stats_()
fig_f1d5_1000_ <- plots_f1d5_1000_ |>
  organize_plots()
webpea(
  "./img/FigureA_5_1000_noNpWeights.webp", plot = fig_f1d5_1000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

# *** N = 20000 simulations

toplot_f1d1_20000 <- stats_f1d1_20000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d1_20000 <- toplot_f1d1_20000 |>
  visualize_sim_stats()
fig_f1d1_20000 <- plots_f1d1_20000 |>
  organize_plots()
webpea(
  "./img/FigureA_1_20000.webp", plot = fig_f1d1_20000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d1_20000_ <- toplot_f1d1_20000 |>
  visualize_sim_stats_()
fig_f1d1_20000_ <- plots_f1d1_20000_ |>
  organize_plots()
webpea(
  "./img/FigureA_1_20000_noNpWeights.webp", plot = fig_f1d1_20000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d2_20000 <- stats_f1d2_20000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d2_20000 <- toplot_f1d2_20000 |>
  visualize_sim_stats()
fig_f1d2_20000 <- plots_f1d2_20000 |>
  organize_plots()
webpea(
  "./img/FigureA_2_20000.webp", plot = fig_f1d2_20000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d2_20000_ <- toplot_f1d2_20000 |>
  visualize_sim_stats_()
fig_f1d2_20000_ <- plots_f1d2_20000_ |>
  organize_plots()
webpea(
  "./img/FigureA_2_20000_noNpWeights.webp", plot = fig_f1d2_20000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d3_20000 <- stats_f1d3_20000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d3_20000 <- toplot_f1d3_20000 |>
  visualize_sim_stats()
fig_f1d3_20000 <- plots_f1d3_20000 |>
  organize_plots()
webpea(
  "./img/FigureA_3_20000.webp", plot = fig_f1d3_20000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d3_20000_ <- toplot_f1d3_20000 |>
  visualize_sim_stats_()
fig_f1d3_20000_ <- plots_f1d3_20000_ |>
  organize_plots()
webpea(
  "./img/FigureA_3_20000_noNpWeights.webp", plot = fig_f1d3_20000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d4_20000 <- stats_f1d4_20000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d4_20000 <- toplot_f1d4_20000 |>
  visualize_sim_stats()
fig_f1d4_20000 <- plots_f1d4_20000 |>
  organize_plots()
webpea(
  "./img/FigureA_4_20000.webp", plot = fig_f1d4_20000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d4_20000_ <- toplot_f1d4_20000 |>
  visualize_sim_stats_()
fig_f1d4_20000_ <- plots_f1d4_20000_ |>
  organize_plots()
webpea(
  "./img/FigureA_4_20000_noNpWeights.webp", plot = fig_f1d4_20000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)

toplot_f1d5_20000 <- stats_f1d5_20000 |>
  select(!starts_with("var_")) |>
  elongate_sim_stats()
plots_f1d5_20000 <- toplot_f1d5_20000 |>
  visualize_sim_stats()
fig_f1d5_20000 <- plots_f1d5_20000 |>
  organize_plots()
webpea(
  "./img/FigureA_5_20000.webp", plot = fig_f1d5_20000,
  width = 1200, height = 630, dpi = 72, units = "px"
)
plots_f1d5_20000_ <- toplot_f1d5_20000 |>
  visualize_sim_stats_()
fig_f1d5_20000_ <- plots_f1d5_20000_ |>
  organize_plots()
webpea(
  "./img/FigureA_5_20000_noNpWeights.webp", plot = fig_f1d5_20000_,
  width = 1200, height = 630, dpi = 72, units = "px"
)
