library(tidyverse)

## Given boolean 'df' calculate the jaccard index between columns 'i' and 'j'. 
.jaccard <- function(df, i, j) {
    s <- summarize(df,
                   union = sum(.data[[i]] | .data[[j]]),
                   inter = sum(.data[[i]] & .data[[j]]))
    s$inter / s$union
}

## Given boolean 'df' calculate the "asymmetric jaccard index" between columns
## 'i' and 'j'. The "asymmetric jaccard index" is like the jaccard index except
## the denominator is the cardinality of the set denoted by 'i'.
.asymmetric_jaccard <- function(df, i, j) {
    s <- summarize(df,
                   total = sum(.data[[i]]),
                   inter = sum(.data[[i]] & .data[[j]]))
    s$inter / s$total
}

jaccard <- Vectorize(.jaccard, vectorize.args = c("i", "j"))
ajaccard <- Vectorize(.asymmetric_jaccard, vectorize.args = c("i", "j"))

to_binary <- . %>% mutate(across(everything(), ~ !is.na(.x)))

## Given a vector 'xs' return a tibble with all combinations of the vector in
## columns 'var.x' and 'var.y'
cross_tibble <- function(xs) {
    .tbl <- tibble(var = xs) %>%
        mutate(.dummy = 1)
    full_join(.tbl, .tbl, by = ".dummy") %>%
        select(-.dummy)
}

## Given 'df' make a tile plot with varibles 'x' and 'y' on their respective
## axis with 'metric' as the fill
make_xy_tile_plot <- function(df, x, y, metric, xlab = NULL, ylab = NULL) {
    ggplot(df, aes_string(x, y, fill = metric)) +
        geom_tile() +
        xlab(xlab) +
        ylab(ylab) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = "bottom")
}

## Given 'df' return a correlation plot for all columns
make_cor_plot <- function(df) {
    df %>%
        cor() %>%
        as_tibble(rownames = "var1") %>%
        gather(-var1, key = var2, value = r) %>%
        make_xy_tile_plot("var1", "var2", "r")
}

## Given a 'df' and two variables 'var1' and 'var2' return the subsets whose
## self-cartesian products produce the combinations denoted by 'var1' and
## 'var2'.
##
## NOTE: this works by assuming that the input is the full self-cartesian
## product of one or more subsets concatenated together. It "might" work on
## incomplete self-cartesian products, but I haven't tested and don't feel like
## thinking about it as I will only be dealing with full sets.
##
## 1. Filter combinations having jaccard index of 1
## 2. Sort lexically
## 3. Remove all combinations where the second in the pair is "greater" than
##    the first (assumption being that this removes all redundant combinations,
##    including all where the pair contains two identical parameters)
## 4. Since the first column represents each unique subset, group by this and 
##    magically transform to a list
deconvolve_products <- function(df, var1, var2) {
    df %>%
        select(all_of(c(var1, var2))) %>%
        arrange(.data[[var1]], .data[[var2]]) %>%
        filter(.data[[var1]] < .data[[var2]]) %>%
        filter(!.data[[var1]] %in% .data[[var2]]) %>%
        group_split(.data[[var1]]) %>%
        map(~ c(.x[[var1]][1], .x[[var2]]))
}

perfect_overlapping_sets <- function(df_bool, combinations, var1, var2) {
    combinations %>%
        arrange(.data[[var1]], .data[[var2]]) %>%
        filter(.data[[var1]] < .data[[var2]]) %>%
        mutate(jaccard = jaccard(df_bool, .data[[var1]], .data[[var2]])) %>%
        filter(jaccard == 1) %>%
        filter(!.data[[var1]] %in% .data[[var2]]) %>%
        group_split(.data[[var1]]) %>%
        map(~ c(.x[[var1]][1], .x[[var2]]))
}
