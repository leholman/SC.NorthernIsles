library(DescTools)

# https://github.com/hyu-ub/perk
# devtools::install_github("hyu-ub/perk/perk")
library(perk)
library(tidyverse)
library(furrr)


dmg_fwd_CCC <- function(df, smp, ci = "z-transform", nperm = NULL, nproc = 1) {
    # Register a parallel backend to use
    options(future.globals.maxSize = 2291289600)
    plan(multisession, workers = nproc) # Choose the number of workers based on your CPU

    # Parallelized code
    results <- future_map_dfr(smp, function(X) {
        dist_data <- df |>
            filter(label == X) |>
            select(tax_name, label, starts_with("fwf")) |>
            pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -label)) |>
            mutate(x = gsub("fwf", "", type)) |>
            select(-type) |>
            mutate(x = as.numeric(x))
        dist_fit <- df |>
            filter(label == X) |>
            select(tax_name, label, matches("^fwdx\\d+")) |>
            pivot_longer(names_to = "type", values_to = "Dx_fwd", c(-tax_name, -label)) |>
            mutate(x = gsub("fwdx", "", type)) |>
            select(-type) |>
            mutate(x = as.numeric(x))
        data <- dist_data |>
            inner_join(dist_fit, by = join_by(tax_name, label, x)) |>
            group_by(tax_name) |>
            arrange(x, .by_group = TRUE)
        if (nrow(dist_data) == 0) {
            return(NULL)
        } else {
            fits <- data |>
                do({
                    ccc <- DescTools::CCC(.$f_fwd, .$Dx_fwd, ci = ci)
                    if (!is.null(nperm)) {
                        ccc_perm <- perk::perk_test(.$f_fwd, .$Dx_fwd, B = nperm, method = "ccc", alternative = "greater")
                        tibble(rho_c = ccc$rho.c$est, rho_lwr_ci = ccc$rho.c$lwr.ci, rho_upr_ci = ccc$rho.c$upr.ci, C_b = ccc$C.b, l_shift = ccc$l.shift, s_shift = ccc$s.shift, label = X, rho_c_perm = ccc_perm$estimate, rho_c_perm_pval = ccc_perm$p.value)
                    } else {
                        tibble(rho_c = ccc$rho.c$est, rho_lwr_ci = ccc$rho.c$lwr.ci, rho_upr_ci = ccc$rho.c$upr.ci, C_b = ccc$C.b, l_shift = ccc$l.shift, s_shift = ccc$s.shift, label = X)
                    }
                })
        }
    }, .progress = TRUE, .options = furrr_options(seed = T)) # Optional progress bar

    # Reset the default plan if needed
    plan(sequential)
    return(results)
}


get_dmg_decay_fit <- function(df, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
    df_dx_fwd <- df |>
        select(tax_name, label, matches("^fwdx\\d+")) |>
        pivot_longer(names_to = "type", values_to = "Dx_fwd", c(-tax_name, -label)) |>
        mutate(x = gsub("fwdx", "", type)) |>
        select(-type)

    df_dx_rev <- df |>
        select(tax_name, label, matches("^bwdx\\d+")) |>
        pivot_longer(names_to = "type", values_to = "Dx_rev", c(-tax_name, -label)) |>
        mutate(x = gsub("bwdx", "", type)) |>
        select(-type)

    df_dx_conf_fwd <- df |>
        select(tax_name, label, starts_with("fwdxConf")) |>
        pivot_longer(names_to = "type", values_to = "Dx_fwd_conf", c(-tax_name, -label)) |>
        mutate(x = gsub("fwdxConf", "", type)) |>
        select(-type)

    df_dx_conf_rev <- df |>
        select(tax_name, label, starts_with("bwdxConf")) |>
        pivot_longer(names_to = "type", values_to = "Dx_rev_conf", c(-tax_name, -label)) |>
        mutate(x = gsub("bwdxConf", "", type)) |>
        select(-type)


    df_fit_fwd <- df |>
        select(tax_name, label, starts_with("fwf")) |>
        pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -label)) |>
        mutate(x = gsub("fwf", "", type)) |>
        select(-type)

    df_fit_rev <- df |>
        select(tax_name, label, starts_with("bwf")) |>
        pivot_longer(names_to = "type", values_to = "f_rev", c(-tax_name, -label)) |>
        mutate(x = gsub("bwf", "", type)) |>
        select(-type)

    dat <- df_dx_fwd |>
        inner_join(df_dx_rev, by = join_by(tax_name, label, x)) |>
        inner_join(df_dx_conf_fwd, by = join_by(tax_name, label, x)) |>
        inner_join(df_dx_conf_rev, by = join_by(tax_name, label, x)) |>
        inner_join(df_fit_fwd, by = join_by(tax_name, label, x)) |>
        inner_join(df_fit_rev, by = join_by(tax_name, label, x)) |>
        inner_join(df |> select(label) |> distinct(), by = join_by(label)) |>
        mutate(x = as.numeric(x)) |>
        filter(x <= pos) |>
        # inner_join(cdata |> select(label, member_unit)) |>
        rowwise() |>
        mutate(
            Dx_fwd_min = f_fwd - Dx_fwd_conf,
            Dx_fwd_max = f_fwd + Dx_fwd_conf,
            Dx_rev_min = f_rev - Dx_rev_conf,
            Dx_rev_max = f_rev + Dx_rev_conf
        )

    # fwd_max <- dat |>
    #     group_by(as.character(x)) |>
    #     summarise(val = mean(Dx_std_fwd) + sd(Dx_std_fwd)) |>
    #     pull(val) |>
    #     max()
    # fwd_min <- dat |>
    #     group_by(as.character(x)) |>
    #     summarise(val = mean(Dx_std_fwd) - sd(Dx_std_fwd)) |>
    #     pull(val) |>
    #     min()
    # rev_max <- dat |>
    #     group_by(as.character(x)) |>
    #     summarise(val = mean(Dx_std_rev) + sd(Dx_std_rev)) |>
    #     pull(val) |>
    #     max()
    # rev_min <- dat |>
    #     group_by(as.character(x)) |>
    #     summarise(val = mean(Dx_std_rev) - sd(Dx_std_rev)) |>
    #     pull(val) |>
    #     min()


    samples <- dat |>
        select(label, tax_name) |>
        distinct() |>
        ungroup()

    smooth_it <- function(x, y, n = 1000, method = "natural") {
        t <- seq_along(x)
        new_t <- seq(min(t), max(t), length.out = n)
        new_x <- spline(t, x, xout = new_t, method = method)$y
        new_y <- spline(t, y, xout = new_t, method = method)$y
        data.frame(t = new_t, x = new_x, y = new_y)
    }


    if (orient == "fwd") {
        dat1 <- pmap_dfr(samples, function(...) {
            current <- tibble(...)
            d <- dat |>
                filter(label == current$label, tax_name == current$tax_name)
            d <- smooth_it(d$x, d$Dx_fwd, n = 1000, method = "natural")
            d <- d |> mutate(
                x = as.numeric(x),
                label = current$label,
                tax_name = current$tax_name,
            )
        })
        ggplot() +
            geom_ribbon(data = dat, aes(x, ymin = Dx_fwd_min, ymax = Dx_fwd_max, group = interaction(label, tax_name)), alpha = 0.1, fill = "#284B63") +
            geom_path(data = dat, aes(x, f_fwd, group = interaction(label, tax_name)), color = "#284B63", linewidth = 1.5) +
            geom_point(data = dat, aes(x, Dx_fwd), alpha = .3, size = 2, fill = "black") +
            # stat_summary(data = dat, aes(x, f_fwd), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 1, fill = "black") +
            # stat_summary(data = dat, aes(x, f_fwd), fun.data = mean_sd, geom = "errorbar", alpha = .3, size = 1, fill = "black") +
            # stat_summary(data = dat, aes(x, f_fwd), fun = mean, geom = "line", size = 0.8, color = "black") +
            xlab("Position") +
            ylab("Frequency") +
            theme_bw() +
            # scale_y_continuous(limits = c(y_min, y_max), breaks = p_breaks) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(), axis.line = element_line(),
                text = element_text(size = 12),
            )
    } else {
        dat1 <- pmap_dfr(samples, function(...) {
            current <- tibble(...)
            d <- dat |>
                filter(label == current$label, tax_name == current$tax_name)
            d <- smooth_it(d$x, d$Dx_rev, n = 1000, method = "natural")
            d <- d |> mutate(
                x = as.numeric(x),
                label = current$label,
                tax_name = current$tax_name,
            )
        })
        ggplot() +
            geom_ribbon(data = dat, aes(x, ymin = Dx_rev_min, ymax = Dx_rev_max, group = interaction(label, tax_name)), alpha = 0.1, fill = "#B04035") +
            geom_path(data = dat, aes(x, f_rev, group = interaction(label, tax_name)), color = "#B04035", linewidth = 1.5, linejoin = "round") + # geom_point(data = dat, aes(x, f_rev), alpha = .3, size = 2, fill = "black") +
            # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 1, fill = "black") +
            # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "errorbar", alpha = .3, size = 1, fill = "black") +
            # stat_summary(data = dat, aes(x, f_rev), fun = mean, geom = "line", size = 0.8, color = "black") +
            xlab("Position") +
            ylab("Frequency") +
            theme_bw() +
            scale_x_continuous(trans = "reverse") +
            scale_y_continuous(position = "right") +
            # scale_y_continuous(limits = c(y_min, y_max), position = "right", breaks = p_breaks) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 12),
            )
    }
}


calculate_plot_grid <- function(n) {
    cols <- ceiling(sqrt(n))
    rows <- ceiling(n / cols)

    list(rows = rows, cols = cols)
}

penalized_weighted_median <- function(values, weights, ratios) {
    # Combine values, weights, and ratios into a single dataframe
    data <- data.frame(values, weights, ratios)

    # Remove rows with NA in any column
    data <- na.omit(data)

    # Create a new weighting scheme that penalizes the values
    # Here we multiply each value by its square and by the number of reads
    data$penalized_weights <- data$ratios^2 * data$weights

    # Order by values
    data <- data[order(data$values), ]

    # Calculate the cumulative penalized weights
    data$cum_weights <- cumsum(data$penalized_weights)

    # Total sum of penalized weights
    total_weights <- sum(data$penalized_weights)

    # Find the median based on the cumulative penalized weights
    median_index <- min(which(data$cum_weights >= total_weights / 2))
    return(data$values[median_index])
}
