#----- build M predictor matrix for QAP 
# create matrix with 1 for mother-daughter dyads, 0 otherwise.
#' Build a symmetric mother-daughter predictor matrix.
#'
#' @param association_matrix Association matrix whose row names define the bats.
#' @param edges Edge data frame containing `year`, `pair_type`, `from`, and `to`.
#' @param the_year Year for which to select mother-daughter edges.
#' @return A square matrix with 1 for mother-daughter dyads and 0 otherwise.
build_m_matrix <- function(association_matrix, edges, the_year) {

  ids <- rownames(association_matrix)
  #edges for this year + the desired pair_type
  edges_subset <- subset(edges, year == the_year & pair_type == "mother-daughter")
  
  #M matrix of zeros with same row/col
  M <- matrix(0, nrow = length(ids), ncol = length(ids),
              dimnames = list(ids, ids))
  
  #M matrix with 1 for mother–daughter
  for (r in seq_len(nrow(edges_subset))) {
    bat_i <- edges_subset$from[r]
    bat_j <- edges_subset$to[r]
    #only if both bats
    if (bat_i %in% ids && bat_j %in% ids) {
      M[bat_i, bat_j] <- 1
      M[bat_j, bat_i] <- 1
    }
  }
  
  return(M)
}


#-----  run QAP for all years 
# separately for each year - 10 000 permutations
#' Run a QAP correlation test for each network year.
#'
#' @param net_list Named list of association matrices, one per year.
#' @param edges Edge data frame used to build mother-daughter predictor matrices.
#' @return A data frame with the correlation and upper- and lower-tail p-values
#'   for each year.
run_qap <- function(net_list,edges){
  results_list <- list()
  
  #for each year
  for (the_year in names(net_list)) {
    A_yr <- net_list[[the_year]]
  
    #build predictor matrix
    M <- build_m_matrix(A_yr, edges, the_year) 
    g <- list(A_yr, M)
    
    #run QAP
    qap_res <- qaptest(g,FUN = gcor,g1  = 1,g2  = 2,reps = 10000)
    
    #extract statistics
    s <- summary(qap_res)
    test_val <- s$testval
    p_greq   <- s$pgreq
    p_leeq   <- s$pleeq
    
  #create result table 
    results_df <- data.frame(
      year        = the_year,
      correlation = test_val,
      pval_upper  = p_greq,
      pval_lower  = p_leeq,
      stringsAsFactors = FALSE
    )
    
    results_list[[the_year]] <- results_df
  }
  
  all_results_df <- do.call(rbind, results_list)
  rownames(all_results_df) <- NULL
  
  return(all_results_df)
}


# add_sig_stars add significance stars to QAP results based on p-values
add_sig_stars <- function(qap_df) {
  qap_df %>%
    mutate(star = case_when(
      pval_upper < 0.001 ~ "***",
      pval_upper < 0.01  ~ "**",
      pval_upper < 0.05  ~ "*",
      TRUE            ~ ""
    ))
}

# make_qap_plot creates a boxplot of association strengths by year and pair type, with significance stars added.
make_qap_plot <- function(edges, qap_df, colony_label) {
  edges_plot <- edges %>%
    mutate(year = as.character(year),
           pair_type = factor(pair_type,
                              levels = c("other", "mother-daughter")))
  stars <- add_sig_stars(qap_df) %>%
    mutate(year = as.character(year))

  ggplot(edges_plot, aes(x = year, y = weight, fill = pair_type)) +
    geom_boxplot(width = 0.6, outlier.size = 0.5, outlier.colour = "#c4c2c2") +
    geom_text(data = stars, aes(x = year, y = 0.97, label = star),
              inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c("mother-daughter" = "#402c3a",
                                 "other"           = "#c4c2c2")) +
    labs(title = colony_label, x = NULL,
         y = "Association strength (SRI)", fill = "pair type") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
}


# compute_perm_pvalues computes node-label permutation p-values (two-tailed)
compute_perm_pvalues <- function(obs_coef, perm_matrix) {
  sapply(names(obs_coef), function(nm)
    mean(abs(perm_matrix[, nm]) >= abs(obs_coef[nm]))
  )
}


# make_perm_table creates a table comparing observed and permuted distributions for permutation tests 
make_perm_table <- function(obs_coef, perm_matrix, p_values) {
  obs_df <- data.frame(term = names(obs_coef),
                       observed = as.numeric(obs_coef),
                       stringsAsFactors = FALSE)
  perm_summary <- as.data.frame(perm_matrix) %>%
    pivot_longer(everything(), names_to = "term", values_to = "value") %>%
    group_by(term) %>%
    summarise(perm_mean = mean(value), perm_sd = sd(value), .groups = "drop")
  obs_df %>%
    left_join(perm_summary, by = "term") %>%
    mutate(p_value = p_values[term]) %>%
    # node-label permutations only test terms involving pair type,
    # covariates (colony size, colony ID) are assessed with the model's Wald tests
    filter(grepl("pair_type", term)) %>%
    arrange(desc(observed))
}


# check_model performs diagnostic checks on a fitted model, including overdispersion, collinearity, and residual simulation.

check_model <- function(model, label) {
  cat("\n──", label, "──\n")
  print(check_overdispersion(model))
  print(check_collinearity(model))
  sim <- simulateResiduals(model, n = 1000)
  plot(sim, main = label)
}


# get_emmeans_contrasts computes pairwise contrasts of estimated marginal means for a fitted model, focusing on the `pair_type_2` factor.
get_emmeans_contrasts <- function(model) {

  emm <- emmeans(model, ~ pair_type_2)

  as.data.frame(
    pairs(emm, adjust = "tukey")
  ) %>%
    mutate(
      # Remove parentheses added by emmeans
      contrast_clean = gsub("\\(|\\)", "", contrast),

      # Extract the two original pair types
      level1 = trimws(sub(" - .*", "", contrast_clean)),
      level2 = trimws(sub(".* - ", "", contrast_clean)),

      # Rename them EXACTLY as they appear in the plot
      level1 = dplyr::recode(
        level1,
        "mother-daughter" = "mother-daughter",
        "mother-unrelated-juvenile" = "mother-unrelated juvenile",
        "daughter-daughter" = "juvenile-juvenile",
        "other" = "other"
      ),

      level2 = dplyr::recode(
        level2,
        "mother-daughter" = "mother-daughter",
        "mother-unrelated-juvenile" = "mother-unrelated juvenile",
        "daughter-daughter" = "juvenile-juvenile",
        "other" = "other"
      ),

      # Tukey-adjusted p-value
      p_value = p.value
    ) %>%

    # Don't plot comparisons involving the reference "other"
    filter(
      level1 != "other",
      level2 != "other"
    ) %>%

    select(level1, level2, p_value)
}


# perm_contrast, compute permutation p-value for a difference between two coefficients
perm_contrast <- function(perm_matrix, obs_coef, term1, term2) {

  obs_diff <- obs_coef[term1] - obs_coef[term2]
  perm_diff <- perm_matrix[, term1] - perm_matrix[, term2]

  p_value <- mean(abs(perm_diff) >= abs(obs_diff))

  data.frame(
    contrast = paste(term1, "vs", term2),
    observed_diff = obs_diff,
    perm_mean = mean(perm_diff),
    perm_sd = sd(perm_diff),
    p_value = p_value
  )
}


# network plot, original simple version, circles only
plot_network_year <- function(edge_df, target_year,colony_name) {

  df_year <- edge_df %>% filter(year == target_year, weight > 0) # only pairs that co-roosted

  df_year$pair_type_2 <- factor(df_year$pair_type_2,
                               levels = c("other", "mother-unrelated-juvenile",
                                          "daughter-daughter", "mother-daughter"),
                               labels = c("unrelated females", "mother-unrelated juvenile",
                                          "juvenile-juvenile", "mother-daughter"))

  g <- graph_from_data_frame(df_year, directed = FALSE)
  E(g)$weight <- df_year$weight

  set.seed(23)
  layout <- layout_with_fr(g, weights = E(g)$weight)

  nodes <- data.frame(name = V(g)$name, x = layout[, 1], y = layout[, 2])

  el <- as_edgelist(g) %>% as.data.frame() %>% rename(from = V1, to = V2)
  edges_full <- el %>%
  mutate(pair_type_2 = df_year$pair_type_2,
         weight = df_year$weight) %>%
  left_join(nodes %>% rename(from = name, x1 = x, y1 = y), by = "from") %>%
  left_join(nodes %>% rename(to   = name, x2 = x, y2 = y), by = "to")

  e_other <- edges_full %>% filter(pair_type_2 == "unrelated females")
  e_murd  <- edges_full %>% filter(pair_type_2 == "mother-unrelated juvenile")
  e_juv   <- edges_full %>% filter(pair_type_2 == "juvenile-juvenile")
  e_md    <- edges_full %>% filter(pair_type_2 == "mother-daughter")

  ggplot() +
    geom_segment(data = e_other, aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = weight*2),
               colour = pair_colors["unrelated females"], alpha = 0.15) +
    geom_segment(data = e_murd,  aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = weight*2),
               colour = pair_colors["mother-unrelated juvenile"], alpha = 0.6) +
    geom_segment(data = e_juv,   aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = weight*2),
               colour = pair_colors["juvenile-juvenile"], alpha = 0.7) +
    geom_segment(data = e_md,    aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = weight*2),
               colour = pair_colors["mother-daughter"], alpha = 0.8) +
  scale_linewidth(range = c(0.2, 3), guide = "none") +
    labs(title = colony_name) +
  geom_point(data = nodes, aes(x = x, y = y), size = 1.3, colour = "#4f4f4f") +
  theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0, size = 9, face = "plain",
                                 margin = margin(b = 5)))
}

# coefficient forest plot, vs "other" reference
sig_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}


# plot_coef_panel creates a coefficient plot with confidence intervals and significance brackets for pairwise contrasts.
plot_coef_panel <- function(model, perm_results, contrasts_df) {

  tidy_model <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

  coef_df <- tidy_model %>%
    filter(grepl("pair_type_2", term)) %>%
    left_join(perm_results %>% select(term, p_value), by = "term") %>%
    mutate(pair_label = gsub("pair_type_2", "", term),
           pair_label = dplyr::recode(pair_label,
                                "mother-daughter" = "mother-daughter",
                                "mother-unrelated-juvenile" = "mother-unrelated juvenile",
                                "daughter-daughter" = "juvenile-juvenile"),
           sig = sig_stars(p_value))

  coef_df$pair_label <- factor(coef_df$pair_label,
                                levels = c("mother-unrelated juvenile",
                                           "juvenile-juvenile",
                                           "mother-daughter"))

  y_pos <- c("mother-unrelated juvenile" = 1,
              "juvenile-juvenile" = 2,
              "mother-daughter" = 3)

  bracket_x_start <- max(coef_df$conf.high) + 0.08

  # stagger brackets by span so they don't overlap
  # span 1 = adjacent levels, span 2 = skip one level
  brackets <- contrasts_df %>%
  mutate(y1 = y_pos[level1],
         y2 = y_pos[level2],
         sig = sig_stars(p_value),
         span = abs(y_pos[level1] - y_pos[level2]),
         ymid = (y_pos[level1] + y_pos[level2]) / 2) %>%
    filter(sig != "") %>%
    arrange(span) %>%
    mutate(x = bracket_x_start + (row_number() - 1) * 0.12)

  ggplot(coef_df, aes(x = estimate, y = pair_label, color = pair_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.9, linewidth = 1.1) +
    # bracket vertical lines
    geom_segment(data = brackets,
                 aes(x = x, xend = x, y = y1, yend = y2),
                 color = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
    # bracket ticks
    geom_segment(data = brackets,
                 aes(x = x - 0.02, xend = x, y = y1, yend = y1),
                 color = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
    geom_segment(data = brackets,
                 aes(x = x - 0.02, xend = x, y = y2, yend = y2),
                 color = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
    # stars at midpoint
    geom_text(data = brackets,
              aes(x = x + 0.02, y = ymid, label = sig),
              color = "black", size = 4, hjust = 0, inherit.aes = FALSE) +
    scale_color_manual(values = pair_colors) +
    labs(x = "coefficient (vs other pairs)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 60, 5, 5)
    )
}


# make_h2_plot function for plotting predicted association strength (SRI) against colony size, with confidence intervals and pair type differentiation.
make_h2_plot <- function(pred, title, show_y = TRUE) {
  ggplot(pred, aes(x = colony_size, y = SRI,
                      color = pair_type_2, fill = pair_type_2)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = pair_colours, name = "pair type") +
    scale_fill_manual(values  = pair_colours, guide = "none") +
    labs(title = title,
         x = "colony size",
         y = if (show_y) "association strength (SRI)" else NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}


# summarise_edges function to summarise edges by year and pair type, and colony size per year
summarise_edges <- function(edges, colony_label) {
  cat("\nSRI summary by year and pair type —", colony_label, "--\n")
  edges %>%
    group_by(year, pair_type_2) %>%
    summarise(mean_sri = round(mean(weight), 3),
              sd_sri   = round(sd(weight),   3),
              n_dyads  = n(), .groups = "drop") %>%
    print()

  cat("\nColony size per year —", colony_label, "--\n")
  edges %>%
    group_by(year) %>%
    summarise(colony_size = mean(colony_size, na.rm = TRUE), .groups = "drop") %>%
    print()

  cat("\nNumber of MD pairs per year —", colony_label, "--\n")
  edges %>%
    filter(pair_type_2 == "mother-daughter") %>%
    group_by(year) %>%
    summarise(n_md = n(), .groups = "drop") %>%
    print()
}



# count dyads per pair type per year, both colonies
count_pairs <- function(edge_df, colony_name) {
  edge_df %>%
    group_by(year, pair_type_2) %>%
    summarise(n_dyads = n_distinct(pairID), .groups = "drop") %>%
    mutate(colony = colony_name)
}


# mean_density computes the mean network density for a list of networks and prints the result with a label.
mean_density <- function(net_list, label) {
  d <- sapply(net_list, function(net) {
    g <- igraph::graph_from_adjacency_matrix(net, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
    igraph::edge_density(g, loops = FALSE)
  })
  cat("mean network density —", label, ":", round(mean(d), 4), "\n")
}


# prepare_edges function to filter and prepare edge data for analysis, including calculating pair origin, relation age, and beta transformation of weights.
prepare_edges <- function(edges_full, md_file, mom_juv, juv_pairs, colony_size_df) {

  # sorted pairID
  df <- edges_full %>%
    mutate(pairID = if_else(from < to,
                            paste(from, to, sep = "-"),
                            paste(to, from, sep = "-")))

  # remove pairs of individuals that were both adults in 2013 (relationship age unknown)
  ind_2013    <- unique(c(df$from[df$year == 2013], df$to[df$year == 2013]))
  adults_2013 <- setdiff(ind_2013, md_file$rfid_off[md_file$year == 2013])
  df <- df %>%
    filter(!(from %in% adults_2013 & to %in% adults_2013))

  # one pair type per pair
  muj_pairs <- unique(mom_juv$pairID[mom_juv$pair_type_2 == "mother-unrelated-daughter"])
  df <- df %>%
    mutate(pair_origin = case_when(
      pairID %in% md_file$pairID   ~ "mother-daughter",
      pairID %in% muj_pairs        ~ "mother-unrelated-juvenile",
      pairID %in% juv_pairs$pairID ~ "juvenile-juvenile",
      TRUE                         ~ "other"
    ))

  # daughter age for MD pairs
  df <- df %>%
    mutate(
      first_year_md = md_file$year[match(pairID, md_file$pairID)],
      daughter_age  = if_else(
        is.na(first_year_md), NA_real_, as.numeric(year - first_year_md + 1)
      )
    )

  # relation age variable
  df <- df %>%
    group_by(pairID) %>%
    mutate(first_year_overall = min(year)) %>%
    ungroup() %>%
    mutate(relation_age = if_else(
      pair_origin == "mother-daughter",
      daughter_age,
      as.numeric(year - first_year_overall + 1)
    ))

  # remove first summer of early-life pairs (juveniles marked in August)
  df <- df %>%
    filter(!(pair_origin %in% c("mother-daughter",
                               "juvenile-juvenile",
                               "mother-unrelated-juvenile") &
           relation_age == 1))

  # each pair only once per year
  stopifnot(!anyDuplicated(df[, c("pairID", "year")]))

  # colony size
  df <- df %>%
    left_join(colony_size_df, by = "year") %>%
    mutate(colony_size_c = as.numeric(
      scale(colony_size, center = TRUE, scale = FALSE)
    ))

  # beta transformation, after all filtering
  n <- nrow(df)
  df <- df %>%
    mutate(weight_beta = (weight * (n - 1) + 0.5) / n)

  # factor levels and centre relation_age
  df %>%
    mutate(
      pair_origin = factor(pair_origin,
                           levels = c("other",
                                      "mother-unrelated-juvenile",
                                      "juvenile-juvenile",
                                      "mother-daughter")),
      relation_age_c = relation_age - mean(relation_age, na.rm = TRUE)
    )
}

count_pairs2 <- function(edges_full, edges_final) {
  all_pairs <- edges_full %>%
    mutate(pairID = if_else(from < to, paste(from, to, sep = "-"), paste(to, from, sep = "-")))
  list(
    pairs_total   = n_distinct(all_pairs$pairID),
    pairs_kept    = n_distinct(edges_final$pairID),
    pairs_by_type = edges_final %>% distinct(pairID, pair_origin) %>% count(pair_origin)
  )
}


# run_permutations function to perform node-label permutations for a specified number of iterations, fitting two models (main and interaction) to each permuted dataset and storing the coefficients.
run_permutations <- function(data, formula_main, formula_int, n_perm, groups = "year") {

  
  perm_main <- matrix(NA, n_perm, length(coef_names),     dimnames = list(NULL, coef_names))
  perm_int  <- matrix(NA, n_perm, length(coef_names_int), dimnames = list(NULL, coef_names_int))

  for (i in seq_len(n_perm)) {
    d <- data %>%
      group_by(across(all_of(groups))) %>%
      mutate(pair_origin = sample(pair_origin)) %>%
      ungroup()

    m <- glmmTMB(formula_main, family = beta_family(), data = d)
    perm_main[i, ] <- fixef(m)$cond[coef_names]

    m_int <- glmmTMB(formula_int, family = beta_family(), data = d)
    perm_int[i, ] <- fixef(m_int)$cond[coef_names_int]
  }
  list(main = perm_main, int = perm_int)
}


# compute_perm_pvalues computes node-label permutation p-values (two-tailed) for observed coefficients against a matrix of permuted coefficients.
compute_perm_pvalues <- function(obs_coef, perm_matrix) {
  sapply(names(obs_coef), function(nm)
    mean(abs(perm_matrix[, nm]) >= abs(obs_coef[nm]), na.rm = TRUE)
  )
}


# make_perm_table creates a table comparing observed and permuted distributions for permutation tests, focusing on terms involving pair type.
make_perm_table <- function(obs_coef, perm_matrix, p_values) {
  data.frame(term = names(obs_coef), observed = as.numeric(obs_coef),
             stringsAsFactors = FALSE) %>%
    left_join(
      as.data.frame(perm_matrix) %>%
        pivot_longer(everything(), names_to = "term", values_to = "value") %>%
        group_by(term) %>%
        summarise(perm_mean = mean(value, na.rm = TRUE), perm_sd = sd(value, na.rm = TRUE), .groups = "drop"),
      by = "term"
    ) %>%
    mutate(p_value = p_values[term]) %>%
    # node-label permutations only test terms involving pair type,
    # the main effect of relationship age is assessed with the model's Wald test
    filter(grepl("pair_origin", term)) %>%
    arrange(desc(observed))
}