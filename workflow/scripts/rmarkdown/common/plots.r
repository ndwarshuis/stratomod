## given a vector x, return a dataframe with transformed vectors depending on
## the domain
infer_transform <- function(x) {
    lower <- min(x)
    upper <- max(x)
    if (0 <= lower && upper <= 1) {
        tibble(identity = x)
    } else if (0 < lower) {
        tibble(identity = x, trans_log10 = log10(x))
    } else {
        tibble(identity = x, trans_arsinh = asinh(x - median(x)))
    }
}
