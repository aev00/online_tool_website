# -----------------------------------------------------------------------------
# MODULE
# Functions for generating lotteries
# Author: Cassidy Shubatt <cshubatt@g.harvard.edu
# -----------------------------------------------------------------------------

# Setup -----------------------------------------------------------------------
library(testit) # assert
library(purrr) # map_2
library(stats) # pnorm
library(dplyr) # mutate()
library(magrittr) # %>%
library(glue)

l <- modules::use(here::here("lib", "lottery_funs.R"))

# Global Variables -------------------------------------------------------------
nice_ps <- c(0.01, seq(0.05, 0.95, 0.05), 0.99)

# Functions -------------------------------------------------------------------

base_lottery <- function(evil_probs = 0, lo_x_wt = 0.3, corr_prob = 0.5){
    low_indic <- rbinom(2, 1, lo_x_wt)
    low_exp_draw <- -rexp(2, 0.04)
    high_exp_draw <- rexp(2, 0.04)

    base_xs <- low_indic * low_exp_draw +
        (1 - low_indic) * high_exp_draw

    corr_indic <- rbinom(1, 1, corr_prob)
    if (corr_indic == 1) {
        base_xs[2] <- rnorm(1, base_xs[1], 20)
    }

    if (evil_probs == 1) {
        base_p1 <- 0
        while(base_p1 < 0.01 | base_p1 > 0.99){
            base_p1 <- rnorm(1, mean = 0.5, sd = 0.25) %>% round(digits = 2)
        }
    } else {
        base_p1 <- sample(nice_ps, 1)
    }

    base_p2 <- 1 - base_p1
    base_ps <- c(base_p1, base_p2) %>% round(digits = 2)

    # if probs 0.01/0.99, reclassify as "evil"
    evil_probs <- ((evil_probs == 1)| any(abs(base_ps - 0.01) < 0.001)) %>% as.integer

    base_lottery <- list(x = base_xs, p = base_ps, evil = evil_probs)

    return(base_lottery)
}

base_s_lottery <- function(ev_s, poss_ps = nice_ps, evil_probs = 0) {
    base_x1 <- rnorm(1, ev_s, sd = 20) %>% round

    if (evil_probs == 1) {
        base_p1 <- 0
        while (base_p1 < 0.01 | base_p1 > 0.99) {
            base_p1 <- rnorm(1, mean = 0.5, sd = 0.25) %>% round(digits = 2)
        }
    } else {
        base_p1 <- sample(poss_ps, 1)
    }

    base_p2 <- 1 - base_p1
    base_ps <- c(base_p1, base_p2)

    base_x2 <- ((ev_s - base_x1 * base_p1) / base_p2) %>% round
    base_xs <- c(base_x1, base_x2)

    # if probs 0.01/0.99, reclassify as "evil"
    evil_probs <- ((evil_probs == 1)| any(abs(base_ps - 0.01) < 0.001)) %>% as.integer

    base_lottery <- list(x = base_xs, p = base_ps, evil = evil_probs)

    return(base_lottery)
}

split_states <- function(x_cur, p_cur, states_cur, q_unif = 0.5, evil_probs = 0){
    # q_unif: in [0,1], wt on unif state sampling draw

    range <- max(x_cur) - min(x_cur)

    # choose state for splitting
    splittable_states <- (1:states_cur)[which(p_cur > 0.05)]
    prob_weights <- p_cur[splittable_states]

    if (length(splittable_states) == 1) {
        # sample function can't be given a single value to sample from
        split_index <- splittable_states
    } else {
        # mixture btw uniform sampling probs, p_cur sampling probs
        # weight on uniform sampling probs
        q_unif <- 0.5
        unif_wts_indic <- rbinom(1, 1, q_unif)

        unif_wt_draw <- sample(splittable_states, size = 1)
        prob_wt_draw <- sample(splittable_states, 1, prob = prob_weights)

        split_index <- unif_wts_indic * unif_wt_draw +
            (1 - unif_wts_indic) * prob_wt_draw
    }

    # split state
    x_split <- x_cur[split_index]
    p_split <- p_cur[split_index]

    split_val <- x_split * p_split

    if (evil_probs == 1) {
        perturb_vals <- seq(-0.04, 0.04, 0.01)
        p_perturb <- sample(perturb_vals, 1)

        new_x1 <- rnorm(1, mean = x_split, sd = 0.25 * range) %>% round
        new_p1 <- ((100 * p_split / 2) %>% round) / 100 + p_perturb %>%
            max(0.01) %>%
            round(digits = 2)

        # ensure positive probability left for new_p2
        new_p1 <- ifelse(p_split - new_p1 < 0.01, p_split - 0.01, new_p1)

        new_p2 <- p_split - new_p1
        # hold EV constant
        new_x2 <- ((split_val - new_x1 * new_p1) / new_p2) %>% round
    } else {
        x_perturb_mean <- 0.25 * range
        x_perturb_pos <- rexp(1, rate = 1 / x_perturb_mean)
        if (abs((p_split * 10) - round(p_split * 10)) < 0.01) {
            # if multiple of 0.10, split in half
            new_p1 <- p_split / 2
            new_p2 <- p_split / 2

            new_x1 <- (x_split + x_perturb_pos) %>% round
            new_x2 <- (x_split - x_perturb_pos) %>% round
        } else {
            # multiple of 0.05, round one half up and one half down (no .025s)
            new_p1 <- (p_split / 2) - 0.025
            new_p2 <- (p_split / 2) + 0.025

            rand_sign <- 2 * rbinom(1,1,0.5) - 1
            new_x1 <- (x_split + rand_sign * x_perturb_pos) %>% round
            # hold EV constant
            new_x2 <- ((split_val - new_x1 * new_p1) / new_p2) %>% round
            assert(abs((new_p1 * 20) - round(new_p1 * 20)) < 0.001)
        }
    }
    assert("Split probabilities both positive", (new_p2 > 0) & (new_p1 > 0))

    x_cur[split_index] <- new_x1
    p_cur[split_index] <- new_p1

    x_cur <- append(x_cur, new_x2, after = split_index)
    p_cur <- append(p_cur, new_p2, after = split_index)

    split_lottery <- list(x = x_cur, p = p_cur)

    return(split_lottery)
}

draw_nstates_s <- function(nstates_c, certain_prob = 0.6) {
    certain <- rbinom(1, 1, certain_prob)

    if (nstates_c == 2) {
        nstates_s_uncertain <- 2
    } else {
        nstates_s_uncertain <- sample(2:nstates_c, 1)
    }

    nstates_s <- certain * 1 + (1 - certain) * nstates_s_uncertain

    return(nstates_s)
}

draw_nstates_multi <- function(nstates_c, any_certain) {
    # certain <- ifelse(any_certain, 0, rbinom(1, 1, 0.5))
    certain <- 0

    if (nstates_c == 2) {
        nstates_s_uncertain <- 2
    } else {
        nstates_s_uncertain <- sample(2:3, size = 1, prob = c(0.7, 0.3))
    }

    nstates_s <- certain * 1 + (1 - certain) * nstates_s_uncertain

    return(nstates_s)
}

any_repeats <- function(v1, v2) {
    v1 <- unlist(v1)
    v2 <- unlist(v2)

    v1_repeats <- any(duplicated(v1))
    v2_repeats <- any(duplicated(v2))

    repeats <- v1_repeats | v2_repeats
    
    return(repeats)
}

any_evil <- function(p1, p2) {
    p1 <- unlist(p1)
    p2 <- unlist(p2)

    nice_probs <- c(0, 0.01, seq(0.05, 0.95, 0.05), 0.99, 1) %>% round(2)

    p1_nice <- all(round(p1, 2) %in% nice_probs)
    p2_nice <- all(round(p2, 2) %in% nice_probs)

    nice <- p1_nice & p2_nice

    evil <- !nice
    
    return(evil)
}

dominance_lot <- function(x_c, p_c, c_dominant, evil_probs = 0) {
    if (c_dominant == 1) {
        # get a value less than min(x_c) with prob Q
        # get an interior value less than the Qth quantile with probability 1-Q
        p_min <- p_c[which(x_c == min(x_c))]
        if (evil_probs == 1) {
            q <- runif(1, min = p_min + 0.01, max = 0.99)
        } else {
            acceptable_qs <- nice_ps[nice_ps >= p_min]
            if (length(acceptable_qs) == 0) {
                acceptable_qs <- setdiff(acceptable_qs, 0.99)
            }
            q <- sample(acceptable_qs, 1)
        }

        e1 <- sample(1:5, 1)

        q_x <- l$quantile(round(q, 2), x_c, p_c, min_geq = FALSE)
        max_dist <- q_x - min(x_c)
        assert(max_dist >= 1)

        max_subtr <- min(5, max_dist)

        e2 <- sample(1:max_subtr, 1)

        x1 <- min(x_c) - e1
        x2 <- q_x - e2

        x <- c(x1, x2)
        p <- c(q, 1 - q)
    } else {
        # get a value greater than max(x_c) with prob Q
        # get an interior value greater than the (1-Q)th quantile with probability 1-Q
        p_max <- p_c[which(x_c == max(x_c))]
        if (evil_probs == 1) {
            q <- runif(1, min = p_max, max = 1)
        } else {
            acceptable_qs <- nice_ps[nice_ps >= p_max]
            if (length(acceptable_qs) > 1) {
                acceptable_qs <- setdiff(acceptable_qs, 0.99)
            }
            q <- sample(acceptable_qs, 1)
        }

        e1 <- sample(1:5, 1)

        # 1 - Qth quantile
        q_inv_x <- l$quantile(round(1 - q, 2), x_c, p_c, min_geq = TRUE)
        max_dist <- max(x_c) - q_inv_x

        assert(max_dist >= 1)

        max_add <- min(5, max_dist)

        e2 <- sample(1:max_add, 1)

        x1 <- max(x_c) + e1
        x2 <- q_inv_x + e2

        x <- c(x1, x2)
        p <- c(q, 1 - q)
    }
    dom_lottery <- list(x = x, p = p)

    return(dom_lottery)
}

# Unpack payoffs/probabilities for lottery B into separate features ------------
unpack_column <- function(
    unpack_col, name_scheme, df, ncol = 7
    ){
    n <- nrow(df)
    mat <- matrix(NA, nrow = n, ncol = ncol)
    for (i in 1:n){
        ktemp <- length(unlist(df[[unpack_col]][i]))
        mat[i, 1:ktemp] <- unlist(df[[unpack_col]][i])
    }
    df <- data.frame(mat)
    names(df) <- glue("{name_scheme}_{1:ncol}")

    return(df)
}

# Repack payoffs/probabilities from wide to long
repack_column <- function(
    repacked_col, name_scheme, df, ncol = 7
) {
    df <- df %>%
        mutate(
            repacked = pmap(
                select(., starts_with(name_scheme)), c
            )
        )
    df[[repacked_col]] <- df$repacked

    drop_cols <- c(
        "repacked", glue("{name_scheme}{1:ncol}")
    )

    df <- df %>% select(-all_of(drop_cols))
    return(df)
}

# sort payoff/probability vectors ----------------------------------------------
sort_lottery <- function(to_sort, sort_by) {
    sorted <- unlist(to_sort)[order(unlist(sort_by), decreasing = TRUE)]

    return(list(sorted))
}

# ev difference between best, second best
ev_diff_sb <- function(ev_a, ev_b, ev_c, ev_d, ev_e) {
    evs <- c(ev_a, ev_b, ev_c, ev_d, ev_e) %>% sort %>% rev

    ev_sb <- evs[1] - evs[2]
    return(ev_sb)
}