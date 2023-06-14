# -----------------------------------------------------------------------------
# MODULE
# Functions for defining complexity features in choices_13k data
# Author: Cassidy Shubatt <cshubatt@g.harvard.edu
# To Use: compl <- modules::use("helpers/complexity_feature_functions.R")
# -----------------------------------------------------------------------------

# Setup -----------------------------------------------------------------------
library(pracma)
library(stats)
library(testit)
library(magrittr)
library(dplyr) # mutate()
library(tidyr) # replae_na()
library(glue)
library(purrr) # map2
library(OneR) # bin()

l <- modules::use(here::here("lib", "lottery_funs.R"))

# Build Complexity Features ----------------------------------------------------
build_features <- function(df, incl_tiles = FALSE, lots = c("a", "b")){
  # input:
  # - df with columns x_a, x_b, p_a, p_b, x_a_cor, x_b_cor, p_ij_cor
  # - indicator incl_tiles; if TRUE, builds in quartile vars for continuous vals

  # output: df with complexity features

  # primitives <- c(
  #   glue("x_{lots[1]}"), glue("x_{lots[2]}"),
  #   glue("p_{lots[1]}"), glue("p_{lots[2]}"),
  #   glue("x_{lots[1]}_cor"), glue("x_{lots[2]}_cor"),
  #   glue("p_{lots[1]}{lots[2]}_cor")
  # )

  orig_cols <- names(df)
  if (!("compound" %in% names(df))) {
    df <- df %>%
      mutate(
        compound = 0,
        compound_range = 0
      )
  }

  for (i in 1:length(lots)) {
      lot_i <- lots[i]

      x_i <- df[[glue("x_{lot_i}")]]
      p_i <- df[[glue("p_{lot_i}")]]

      df[[glue("ev__{lot_i}")]] <- map2(
          .x = x_i, .y = p_i, ~l$lottery_value(x_list = .x, p_list = .y)
      ) %>% unlist
  }

  for (i in 1:(length(lots))) {
      lot_i <- lots[i]

      x_i <- df[[glue("x_{lot_i}")]]
      p_i <- df[[glue("p_{lot_i}")]]
      x_i_cor <- df[[glue("x_{lot_i}_cor")]]

      if (i < length(lots)) {
          for (j in (i + 1):length(lots)) {
              lot_j <- lots[j]

              x_j <- df[[glue("x_{lot_j}")]]
              p_j <- df[[glue("p_{lot_j}")]]
              x_j_cor <- df[[glue("x_{lot_j}_cor")]]
              p_ij_cor <- df[[glue("p_{lot_i}{lot_j}_cor")]]

              # ev diff
              ev_i <- df[[glue("ev__{lot_i}")]]
              ev_j <- df[[glue("ev__{lot_j}")]]
              ev_diff_ij <- glue("ev_diff__{lot_i}{lot_j}")
              df[[ev_diff_ij]] <- ev_i - ev_j

              # ev diff transformations
              abs_ev_diff <- glue("abs_{ev_diff_ij}")
              df[[abs_ev_diff]] <- abs(df[[ev_diff_ij]])

              abs_ev_diff_sq <- glue("{abs_ev_diff}_sq")
              df[[abs_ev_diff_sq]] <- df[[abs_ev_diff]]^2

              abs_ev_diff_cu <- glue("{abs_ev_diff}_cu")
              df[[abs_ev_diff_cu]] <- df[[abs_ev_diff]]^3

              sgn_ev_diff_sq <- glue("sgn_{ev_diff_ij}_sq")
              df[[sgn_ev_diff_sq]] <- df[[abs_ev_diff]] * df[[ev_diff_ij]]

              ev_diff_cu <- glue("{ev_diff_ij}_cu")
              df[[ev_diff_cu]] <- df[[ev_diff_ij]]^3

              # dominance
              df[[glue("dom__{lot_i}{lot_j}")]] <- pmap(
                  .l = list(x_i_cor, x_j_cor), .f = dominated
              ) %>% unlist

              # CDF diff
              # df[[glue("cdf_diff_sq__{lot_i}{lot_j}")]] <- pmap(
              #     .l = list(x_i_cor, x_j_cor, p_ij_cor), .f = CDF_diff_squared
              # ) %>% unlist

              cdf_diff_abs_val <- pmap(
                  .l = list(x_i_cor, x_j_cor, p_ij_cor), .f = CDF_diff_abs
              ) %>% unlist
              cdf_diff_abs <- glue("cdf_diff_abs__{lot_i}{lot_j}")
              df[[cdf_diff_abs]] <- round(
                cdf_diff_abs_val - df[[abs_ev_diff]], digits = 1
              )

              # BGS salience
              # df[[glue("bgs_sal1_var__{lot_i}{lot_j}")]]<- pmap(
              #     .l = list(x_i_cor, x_j_cor, p_ij_cor, 1),
              #     .f = salience_wts
              # ) %>% unlist

              # bushong et al range wts
              # df[[glue("range_wts_var__{lot_i}{lot_j}")]]<- pmap(
              #     .l = list(x_i_cor, x_j_cor, p_ij_cor, 1),
              #     .f = range_wts
              # ) %>% unlist

            if (incl_tiles) {
              tile_vars_multi <- c(
                # "cdf_diff_sq",
                "cdf_diff_abs"
              )

              for(tile_var in tile_vars_multi){
                source_var_i <- glue("{tile_var}__{lot_i}{lot_j}")
                tile_var_i <- glue("{source_var_i}_tile")

                df[[tile_var_i]] <- bin(
                  df[[source_var_i]], nbins = 4, labels = 1:4,
                  method = "content"
                )
              }
            }

            # add log versions of cardinal features ----------------------------
            card_feats <- c(
              # "cdf_diff_sq",
              "cdf_diff_abs"
            )

            for (feat in card_feats) {
              feat_i <- glue("{feat}__{lot_i}{lot_j}")
              log_i <- glue("ln_{feat_i}")
              sq_i <- glue("sq_{feat_i}")
              df[[log_i]] <- log(df[[feat_i]] + 1)
              df[[sq_i]] <- df[[feat_i]]^2
            }
          }
      }

      # generate indiv features for lottery i
      df[[glue("mixed__{lot_i}")]] <- map(
          x_i, ~mixed_sign(xs = .x)
      ) %>% unlist

      df[[glue("gains__{lot_i}")]] <- !df[[glue("mixed__{lot_i}")]] &  # nolint
          df[[glue("ev__{lot_i}")]] >= 0

      df[[glue("losses__{lot_i}")]] <- !df[[glue("mixed__{lot_i}")]] & # nolint
          df[[glue("ev__{lot_i}")]] < 0

      df[[glue("evil_probs__{lot_i}")]] <- map(
          p_i, ~evil_probabilities(probs = .x)
      ) %>% unlist

      df[[glue("scale__{lot_i}")]] <- map(
          x_i, ~x_scale(xs = .x)
      ) %>% unlist

      df[[glue("var__{lot_i}")]] <- map2(
          x_i, p_i, ~var_lottery(xs = .x, probs = .y)
      ) %>% unlist

      df[[glue("var_norm__{lot_i}")]] <- map2(
          x_i, p_i, ~var_norm(xs = .x, probs = .y)
      ) %>% unlist

      df[[glue("cdf_diff_self__{lot_i}")]] <- map2(
        x_i, p_i, ~CDF_diff_self(xs = .x, probs = .y)
      ) %>% unlist

      # df[[glue("var_1__{lot_i}")]] <- map2(
      #     x_i, p_i, ~var_1(xs = .x, probs = .y)
      # ) %>% unlist

      # df[[glue("var_2__{lot_i}")]] <- map2(
      #     x_i, p_i, ~var_2(xs = .x, probs = .y)
      # ) %>% unlist

      # df[[glue("var_3__{lot_i}")]] <- map2(
      #     x_i, p_i, ~var_3(xs = .x, probs = .y)
      # ) %>% unlist

      # df[[glue("var_4__{lot_i}")]] <- map2(
      #     x_i, p_i, ~var_4(xs = .x, probs = .y)
      # ) %>% unlist

      df[[glue("range__{lot_i}")]] <- map(
          x_i, ~x_range(xs = .x)
      ) %>% unlist

      df[[glue("payoff_var__{lot_i}")]] <- map(
              x_i, ~var_general(xs = .x)
          ) %>% unlist

      df[[glue("prob_var__{lot_i}")]] <- map(
              p_i, ~var_general(xs = .x)
          ) %>% unlist

      df[[glue("entropy__{lot_i}")]] <- map(
          p_i, ~entropy(probs = .x)
      ) %>% unlist

      df[[glue("adb__{lot_i}")]] <- map(
          p_i, ~avg_dist_from_boundary(probs = .x)
      ) %>% unlist

      df[[glue("payoff_wtd_db__{lot_i}")]] <- map2(
        x_i, p_i, ~payout_wtd_dist_from_boundary(xs = .x, probs = .y)
      ) %>% unlist

      df[[glue("nstates__{lot_i}")]] <- map(
              x_i, ~num_states(xs = .x)
          ) %>% unlist

      df[[glue("payoff_dispersion__{lot_i}")]] <- map(
              x_i, ~payout_dispersion(xs = .x)
          ) %>% unlist

      # df[[glue("max_x__{lot_i}")]] <- map(
      #     x_i, ~x_max(xs = .x)
      # )  %>% unlist

      # df[[glue("max_p__{lot_i}")]] <- map(
      #     p_i, ~x_max(xs = .x)
      # ) %>% unlist
      # df[[glue("min_x__{lot_i}")]] <- map(
      #     x_i, ~x_min(xs = .x)
      # ) %>% unlist

      # df[[glue("min_p__{lot_i}")]] <- map(
      #     p_i, ~x_min(xs = .x)
      # ) %>% unlist

      df[[glue("certain__{lot_i}")]] <- map(
              .x = p_i, ~any_certain(probs = .x)
          ) %>% unlist

      many_safe <- mean(df[[glue("certain__{lot_i}")]]) > 0.25
      if(incl_tiles & !many_safe){
        # generate quartile factor variables for cont. valued features
        # quartiles won't make sense when many lots are safe payments -- skip
        tile_vars <- c(
          "scale", "var",
          # glue("var_{1:4}"),
          "range", "payoff_var", "entropy", "adb"
        )

        for(tile_var in tile_vars){
          source_var_i <- glue("{tile_var}__{lot_i}")
          tile_var_i <- glue("{source_var_i}_tile")

          df[[tile_var_i]] <- bin(
            df[[source_var_i]], nbins = 4, labels = 1:4, method = "content"
          )
        }
      }
    # add log, square versions of cardinal features ----------------------------
    card_feats <- c(
      "scale", "range", "var",
      # "var_1", "var_2", "var_3", "var_4",
      "payoff_var", "nstates", "adb", "payoff_wtd_db", "entropy", "prob_var",
      "payoff_dispersion", "var_norm", "cdf_diff_self"
      # "min_x", "max_x", "min_p", "max_p"
    )

    for (feat in card_feats) {
      log_i <- glue("ln_{feat}__{lot_i}")
      sq_i <- glue("sq_{feat}__{lot_i}")
      feat_i <- glue("{feat}__{lot_i}")
      df[[log_i]] <- log(abs(df[[feat_i]]) + 1)
      df[[sq_i]] <- df[[feat_i]] ^2
    }
  }

  # add aves of feat_i, feat_j -------------------------------------------------
  ave_feats <- c(
    card_feats, glue("ln_{card_feats}"), glue("sq_{card_feats}"),
    "gains", "evil_probs"
  )
  lot_i <- lots[1]
  lot_j <- lots[2]

  for (ave_feat in ave_feats) {
    ave_name <- glue("ave_{ave_feat}__{lot_i}{lot_j}")
    feat_i <- glue("{ave_feat}__{lots[1]}")
    feat_j <- glue("{ave_feat}__{lots[2]}")

    df[[ave_name]] = map2(
      .x = df[[feat_i]], .y = df[[feat_j]], .f = ~mean(c(.x, .y))
    ) %>% unlist
  }

  # # add summary mixed/gain/loss vars -------------------------------------------
  # gains_i <- glue("gains__{lot_i}")
  # gains_j <- glue("gains__{lot_j}")

  # mixed_i <- glue("mixed__{lot_i}")
  # mixed_j <- glue("mixed__{lot_j}")

  # losses_i <- glue("losses__{lot_i}")
  # losses_j <- glue("losses__{lot_j}")

  # both_gains <- glue("smry_both_gains__{lot_i}{lot_j}")
  # df[[both_gains]] <- df[[gains_i]] & df[[gains_j]]

  # both_mixed <- glue("smry_both_mixed__{lot_i}{lot_j}")
  # df[[both_mixed]] <- df[[mixed_i]] & df[[mixed_j]]

  # both_losses <- glue("smry_both_losses__{lot_i}{lot_j}")
  # df[[both_losses]] <- df[[losses_i]] & df[[losses_j]]

  # gain_loss <- glue("smry_gain_loss__{lot_i}{lot_j}")
  # df[[gain_loss]] <- (df[[gains_i]] & df[[losses_j]]) |
  # (df[[losses_i]] & df[[gains_j]])

  # gain_mixed <- glue("smry_gain_mixed__{lot_i}{lot_j}")
  # df[[gain_mixed]] <- (df[[gains_i]] & df[[mixed_j]]) |
  # (df[[mixed_i]] & df[[gains_j]])

  # loss_mixed <- glue("smry_loss_mixed__{lot_i}{lot_j}")
  # df[[loss_mixed]] <- (df[[losses_i]] & df[[mixed_j]]) |
  # (df[[mixed_i]] & df[[losses_j]])

  ev_feats <- names(df)[grep("ev_", names(df))]
  lasso_feats <- names(df) %>%
    setdiff(orig_cols) %>%
    setdiff(ev_feats) %>%
    union(c("compound", "compound_range"))
  # for each feature, build a "dom" version and "non-dom"
  # these will be used in LASSO and are labeled as such
  dom <- df[[glue("dom__{lots[1]}{lots[2]}")]]
  nondom <- 1 - dom
  for (feat_name in lasso_feats) {
    lasso_feat_dom <- glue("lasso_{feat_name}_dom")
    df[[lasso_feat_dom]] <- dom * df[[feat_name]]

    lasso_feat_nondom <- glue("lasso_{feat_name}_nondom")
    df[[lasso_feat_nondom]] <- nondom * df[[feat_name]]
  }

  return(df)
}

outcomes_same_sign <- function(A_xs, B_xs){
  outcomes <- c(A_xs, B_xs) %>% unlist
  same_sign <- all(outcomes >= 0) | all(outcomes <= 0)
  return(same_sign)
}

correlate_states <- function(x_a, x_b, p_a, p_b, ...) {
  a_df <- data.frame(x = x_a, p = p_a)
  a_df <- a_df[order(a_df$x), ] %>%
    mutate(cdf = cumsum(p) %>% round(digits = 3))

  b_df <- data.frame(x = x_b, p = p_b)
  b_df <- b_df[order(b_df$x), ] %>%
    mutate(cdf = cumsum(p) %>% round(digits = 3))

  cdf_vals <- c(a_df$cdf, b_df$cdf) %>% round(digits = 4) %>% unique
  cdf_vals <- cdf_vals[order(cdf_vals)]

  pay_a <- c()
  pay_b <- c()

  for (i in 1:length(cdf_vals)) {
    cdf_i <- cdf_vals[i]

    a_better <- filter(a_df, cdf >= cdf_i)
    b_better <- filter(b_df, cdf >= cdf_i)
    assert(length(a_better$x) > 0 & length(b_better$x) > 0)
    assert(!(any(is.na(a_better$x))) & !(any(is.na(a_better$x))))

    pay_a_i <- min(a_better$x)
    pay_b_i <- min(b_better$x)

    pay_a <- c(pay_a, pay_a_i)
    pay_b <- c(pay_b, pay_b_i)
  }

  cor_states <- data.frame(state = 1:length(cdf_vals), cdf_vals, pay_a, pay_b) %>%
    mutate(
      p_ab_cor = c(cdf_vals[1], diff(cdf_vals))
    )
  return(cor_states)
}

dominated <- function(x_a_cor, x_b_cor, ...){
  if (any(is.na(c(x_a_cor, x_b_cor)))) {
    return(NA)
  }

  x_a <- unlist(x_a_cor)
  x_b <- unlist(x_b_cor)

  dom <- all(x_a >= x_b) |
    all(x_b >= x_a)

  equal <- all(x_a == x_b)
  dom <- dom & !equal

  return(dom)
}

salience_wts <- function(x_a_cor, x_b_cor, p_ij_cor, theta = 1, ...){
  # calculates bordalo et. al. salience of each state, returns normalized variance
  # high variance = large differences in salience
  # s(x1, x2) = |x1 - x2| / |x1| + |x2| + theta

  if(any(is.na(c(x_a_cor, x_b_cor, p_ij_cor)))){
    return(NA)
  }
  theta <- theta

  x_a <- unlist(x_a_cor)
  x_b <- unlist(x_b_cor)
  p <- unlist(p_ij_cor)

  df <- data.frame(x_a, x_b, p)

  if (length(x_a) == 1 & length(x_b) == 1) {
    var_sal <- 0
  } else {
    cor_states <- df %>%
      mutate(
        sal_num = abs(x_a - x_b),
        sal_denom = abs(x_a) + abs(x_b) + theta,
        sal = sal_num / sal_denom
      )

      var_sal <- var(cor_states$sal)
  }
  return(var_sal)
}

range_wts <- function(x_a_cor, x_b_cor, p_ij_cor, ...) {
  # calculates bushong et al (2021) range normalized weights of states, returns variance of these weights
  # high variance = large differences in range-normalized weights
  # w(x1, x2) = 1 / (1 + |x1 - x2|)

  if (any(is.na(c(x_a_cor, x_b_cor)))) {
    return(NA)
  }

  x_a <- unlist(x_a_cor)
  x_b <- unlist(x_b_cor)

  weights <- 1 / (1 + abs(x_a - x_b))
  var_wt <- var(weights)

  return(var_wt)
}

any_certain <- function(probs){
  if(any(is.na(probs))){
    return(NA)
  }
  certain <- any(unlist(probs) == 1)
  return(certain)
}

mixed_sign <- function(xs){
  if(any(is.na(xs))){
    return(NA)
  }
  xs <- unlist(xs)
  same_sign <- all(xs >= 0) | all(xs <= 0)
  mixed_sign <- !same_sign
  return(mixed_sign)
}

entropy <- function(probs){
  if(any(is.na(probs))){
    return(NA)
  }
  probs <- unlist(probs) %>%
    setdiff(c(0))
  entropy <- (probs * (-log(probs))) %>% sum
  return(entropy)
}

x_sd <- function(xs){
  if(any(is.na(xs))){
    return(NA)
  }
  xs <- unlist(xs)
  return(sd(xs))
}

x_scale <- function(xs){
  if(any(is.na(xs))){
    return(NA)
  }
  xs <- unlist(xs)
  scale <- mean(abs(xs))
  return(scale)
}

num_states <- function(xs){
  xs <- unlist(xs)

  if(all(is.na(xs))){
    return(NA)
  }
  num_states <- length(xs)
  return(num_states)
}

avg_dist_from_boundary <- function(probs){
  probs <- unlist(probs)

  if (all(is.na(probs))) {
    return(NA)
  }

  total_dist <- 0
  n <- length(probs[!is.na(probs)])
  for(i in 1:n){
    total_dist <- total_dist + min(abs(probs[i]), abs(1 - probs[i]))
  }

  avg_dist <- total_dist / n
  return(avg_dist)
}

payout_wtd_dist_from_boundary <- function(xs, probs){
  probs <- unlist(probs)
  xs <- unlist(xs)

  if (all(is.na(probs))) {
    return(NA)
  }
  total_wtd_dist <- 0

  n <- length(probs[!is.na(probs)])
  for(i in seq(n)){
    dist_i <- xs[i] * min(abs(probs[i]), abs(1 - probs[i]))
    total_wtd_dist <- total_wtd_dist + dist_i
  }

  avg_dist <- total_wtd_dist / n
  return(avg_dist)
}

# pweighted_sd <- function(xs, probs){
#   if(any(is.na(c(xs, probs)))){
#     return(NA)
#   }
#   xs <- unlist(xs)
#   probs <- unlist(probs)
  
#   rep_times <- round(probs*100)
#   pweighted_outcomes <- rep(xs, rep_times)
#   len_diff <- abs(length(pweighted_outcomes) - 100)
#   assert("Len pweighted outcomes ~ 100", len_diff < 8)
#   pweighted_sd <- sd(pweighted_outcomes)
  
#   return(pweighted_sd)
# }

var_1 <- function(xs, probs){
  if(any(is.na(c(xs, probs)))){
    return(NA)
  }
  xs <- unlist(xs)
  probs <- unlist(probs)

  xbar <- mean(xs)
  sqs <- (xs - xbar)^2
  sum_sqs <- probs %*% sqs
  var <- sum_sqs# / length(xs)

  return(var)
}

var_2 <- function(xs, probs){
  if(any(is.na(c(xs, probs)))){
    return(NA)
  }
  xs <- unlist(xs)
  probs <- unlist(probs)

  px <- xs * probs
  px_bar <- mean(px)

  sum_sqs <- sum((px - px_bar)^2)
  var <- sum_sqs#/length(xs)

  return(var)
}

var_3 <- function(xs, probs){
  if(any(is.na(c(xs, probs)))){
    return(NA)
  }
  xs <- unlist(xs)
  probs <- unlist(probs)

  pbar <- mean(probs)
  sqs_p <- (probs - pbar)^2

  xbar <- mean(xs)
  sqs_x <- (xs - xbar)^2

  sum_sqs <- sqs_p %*% sqs_x
  var <- sum_sqs# / length(xs)

  return(var)
}

var_4 <- function(xs, probs){
  if(any(is.na(c(xs, probs)))){
    return(NA)
  }
  xs <- unlist(xs)
  probs <- unlist(probs)

  pbar <- mean(probs)
  sqs_p <- (probs - pbar)^2

  sum_sqs <- abs(xs) %*% sqs_p
  var <- sum_sqs# / length(xs)

  return(var)
}

var_norm <- function(xs, probs){
  if (any(is.na(c(xs, probs)))) {
    return(NA)
  }
  xs <- unlist(xs)
  probs <- unlist(probs)

  scale_x <- mean(abs(xs))

  if (scale_x == 0) {
    var_norm <- 0
  } else {
    exp_of_square <- probs %*% xs^2
    square_of_exp <- (probs %*% xs)^2

    var_norm <- (exp_of_square - square_of_exp) / scale_x^2
  }

  return(var_norm)
}

var_general <- function(xs) {
  if (any(is.na(xs))) {
    return(NA)
  }
  xs <- unlist(xs)

  var <- sum((xs - mean(xs))^2 / length(xs))

  return(var)
}

var_lottery <- function(xs, probs) {
  if (any(is.na(c(xs, probs)))) {
    return(NA)
  }
  xs <- unlist(xs)
  probs <- unlist(probs)

  exp_of_square <- probs %*% xs^2
  square_of_exp <- (probs %*% xs)^2

  var <- exp_of_square - square_of_exp
  return(var)
}

CDF_diff_abs <- function(x_a_cor, x_b_cor, p_ij_cor, ...){
  # calculates total absolute area of difference between CDFs for two lotteries

  if (any(is.na(c(x_a_cor, x_b_cor, p_ij_cor)))) {
    return(NA)
  }
  x_a <- unlist(x_a_cor)
  x_b <- unlist(x_b_cor)
  p <- unlist(p_ij_cor)

  df <- data.frame(x_a, x_b, p) %>%
    mutate(
      val_diff = abs(x_a - x_b),
      area_i = val_diff * p,
    )

  area <- sum(df$area_i)

  return(area)
}

CDF_diff_squared <- function(x_a_cor, x_b_cor, p_ij_cor, ...){
  # calculates total absolute area of difference between CDFs for two lotteries
  
  if (any(is.na(c(x_a_cor, x_b_cor, p_ij_cor)))) {
    return(NA)
  }
  x_a <- unlist(x_a_cor)
  x_b <- unlist(x_b_cor)
  p <- unlist(p_ij_cor)

  df <- data.frame(x_a, x_b, p) %>%
    mutate(
      val_diff = abs(x_a - x_b),
      area_i = val_diff * p^2
    )
  
  area <- sum(df$area_i)

  return(area)
}

CDF_diff_self <- function(xs, probs, ...) {
  x <- unlist(xs)
  p <- unlist(probs)

  ev <- x %*% p

  cdf_diff <- 0

  for (i in seq_len(length(x))) {
    x_i <- x[i]
    p_i <- p[i]

    cdf_diff_i <- abs(x_i - ev) * p_i
    cdf_diff <- cdf_diff + cdf_diff_i
  }

  return(cdf_diff)
}

# hellinger_dist <- function(A_xs, B_xs, A_probabilities, B_probabilities, ...){
#   A_xs <- unlist(A_xs)
#   B_xs <- unlist(B_xs)
#   A_probs <- unlist(A_probabilities)
#   B_probs <- unlist(B_probabilities)

#   A_df <- data.frame(xs = A_xs, A_probs)
#   B_df <- data.frame(xs = B_xs, B_probs)

#   joined_df <- full_join(A_df, B_df, by = "xs") %>%
#     mutate(
#       A_probs = replace_na(A_probs, 0),
#       B_probs = replace_na(B_probs, 0)
#       )
  
#   # note: for any lotteries with no payoffs in common, this measure will be equal
#   adj_diffs <- (sqrt(joined_df$A_probs) - sqrt(joined_df$B_probs))^2
#   hellinger <- sum(adj_diffs)/sqrt(2)

#   return(hellinger)
# }

# pweighted_sd_norm <- function(xs, probs){
#   if(any(is.na(c(xs, probs)))){
#     return(NA)
#   }
#   xs <- unlist(xs)
#   if(length(xs) == 1){
#     pweighted_sd <- 0
#   }else{
#     payoff_scalar <- mean(abs(xs))
#     if(payoff_scalar == 0){
#       warning("Mean absolute payoff of 0, returning NA")
#       return(NA)
#       }
#     xs <- xs/payoff_scalar
    
#     probs <- unlist(probs)
    
#     rep_times <- round(probs*1000)
#     pweighted_outcomes <- rep(xs, rep_times)
#     len_diff <- abs(length(pweighted_outcomes) - 1000)
#     assert("Len pweighted outcomes ~ 100", len_diff < 8)
#     pweighted_sd <- sd(pweighted_outcomes)
#   }
  
#   return(pweighted_sd)
# }


prob_concentration <- function(probs){
  if(any(is.na(c(probs)))){
    return(NA)
  }
  probs <- unlist(probs)
  concentration <- prod(probs)
  k <- length(probs)
  conc_root <- nthroot(x = concentration, n = k)
  
  ratio <- conc_root/(1/k)
  return(ratio)
}

x_range <- function(xs){
  if(any(is.na(xs))){
    return(NA)
  }
  xs <- unlist(xs)
  range <- max(xs) - min(xs)

  return(range)
}

x_min <- function(xs){
  if(any(is.na(xs))){
    return(NA)
  }
  xs <- unlist(xs)
  return(min(xs))
}

x_max <- function(xs){
  if(any(is.na(xs))){
    return(NA)
  }
  xs <- unlist(xs)
  return(max(xs))
}

payout_dispersion <- function(xs) {
  if (any(is.na(xs))) {
    return(NA)
  }
  xs <- unlist(xs)

  n <- length(xs)
  xbar <- mean(abs(xs))

  if (xbar == 0) {
    return(0)
  }

  norm_dispersion <- abs(xs - xbar) / xbar
  avg_dispersion <- (1 / n) * sum(norm_dispersion)

  return(avg_dispersion)
}

evil_probabilities <- function(probs) {
  probs <- unlist(probs)

  nice_probs <- c(0, 0.01, seq(0.05, 0.95, 0.05), 0.99, 1) %>% round(2)

  probs_nice <- all(round(probs, 2) %in% nice_probs)

  evil <- !probs_nice

  return(evil)
}