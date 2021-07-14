## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  message = FALSE, 
  warning = FALSE, 
  echo = TRUE
)

## -----------------------------------------------------------------------------
# Load Packages
library("dplyr")  # Data Manipulation and pipe (%>%) 
library("ffp")    # Fully-Flexible Probabilities
# Load Data
data("db_tbl")
db_tbl

## -----------------------------------------------------------------------------
# Inputs in Black-Scholes option pricing
invariants <- db_tbl %>% 
  dplyr::select(VIX, SWAP10YR, `S&P 500`) %>% 
  purrr::map_df(~diff(log(.))) # compute the continuous return for every column

# get last observations
last_obs <- db_tbl %>% 
  dplyr::slice_tail(n = 1) %>% 
  dplyr::select(VIX, SWAP10YR, `S&P 500`)

# Last observation for the underlying (SP500)
S_0   <- last_obs[["S&P 500"]]  
# Last observation for the implied-volatility
vol_0 <- last_obs[["VIX"]]
# Last observation for the risk-free rate
rf_0  <- last_obs[["SWAP10YR"]] 

# Gererate Paths from Historical Scenarios
paths <- purrr::map2_df(.x = last_obs, .y = invariants, .f = ~ .x * exp(.y)) %>%
  dplyr::rename(S_T = `S&P 500`, vol_T = `VIX`, rf_T = `SWAP10YR`) %>% 
  tibble::rowid_to_column(var = "id") 
paths

## -----------------------------------------------------------------------------
inflation <- db_tbl %>% 
  dplyr::select(`10YR Inflation Swap Rate`) %>% 
  dplyr::slice(1:(nrow(db_tbl) - 1))

## -----------------------------------------------------------------------------
call_price <- function(p, k, r, t, s) {
  d_1 <- log(p / k) + (r + s * s / 2) * t
  d_2 <- d_1 - s * sqrt(t)
  c <- p * stats::pnorm(d_1) - k * exp(-r * t) * stats::pnorm(d_2)
  c
}

## -----------------------------------------------------------------------------
N <- 20
# call parameters 
K      <- S_0 * (seq(0.8, 1.1, length = N))
Expiry <- (2:(N + 1)) / 252

## -----------------------------------------------------------------------------
# Market Pricing
pnl <- tibble::tibble(K = K, Expiry = Expiry, panel = as.factor(1:20), S_0, vol_0, rf_0) %>%
  dplyr::mutate(paths = list(paths)) %>%
  tidyr::unnest(paths) %>%
  dplyr::select(id, panel, dplyr::everything()) %>%
  # portfolio scenarios
  dplyr::mutate(
    # Pricing
    call_open  = call_price(p = S_0, k = K, r = rf_0, t = Expiry, s = vol_0),
    call_close = call_price(p = S_T, k = K, r = rf_T, t = Expiry - 1 / 252, s = vol_T),
    pnl        = call_close - call_open,
    # Units to buy and sell
    u          = rep(c(1, -1), each = nrow(paths) * 10)
  ) %>%
  # Aggregate by "day" (here represented by the "id" variable)
  dplyr::group_by(id) %>%
  dplyr::summarise(pnl_u = as.double(pnl %*% u))
pnl

## ---- fig.align='center', fig.height=5, fig.width=7---------------------------
pnl %>% 
  ggplot2::ggplot(ggplot2::aes(x = pnl_u)) + 
  ggdist::stat_histinterval(breaks = 100, outline_bars = TRUE) + 
  ggplot2::scale_x_continuous(labels = scales::dollar_format(prefix = "U$ ")) + 
  ggplot2::labs(title = "P&L Simulation",
                subtitle = "Long-Short Call Options Strategy", 
                x = NULL, 
                y = NULL) 

## ----message=FALSE, warning=FALSE---------------------------------------------
#### Full Information ####
# exponential-smoothing
fp_es1 <- exp_decay(invariants, 0.0166)
fp_es2 <- exp_decay(invariants, 0.0055)
# crisp-conditioning on inflation
fp_cc <- crisp(inflation, lgl = as.logical(inflation >= 2.8))
# normal kernel on inflation
fp_kd <- kernel_normal(inflation, mean = 3, sigma = var(diff(inflation[[1]])))
#### Partial Information ####
# entropy-pooling by kernel-dumping on inflation
fp_ekd <- kernel_entropy(inflation, mean = 3, sigma = var(diff(inflation[[1]])))
# entropy-pooling by moment-matching
fp_emc <- double_decay(invariants, slow = 0.0055, fast = 0.0166)

## ---- echo = FALSE, fig.align='center', fig.height=5, fig.width=7-------------
bind_probs(fp_es1, fp_es2, fp_cc, fp_kd, fp_ekd, fp_emc) %>% 
  ggplot2::ggplot(ggplot2::aes(x = rowid, y = probs, color = key)) + 
  ggplot2::geom_line(show.legend = FALSE) + 
  ggplot2::facet_wrap(~key, 
                      labeller = ggplot2::labeller(
                        key = c("1" = "Exp. Smoothing", "2" = "Exp. Smoothing", 
                                "3" = "Market-Conditioning", "4" = "Normal Kernel", 
                                "5" = "FFP Kernel", "6" = "FFP Double-Decay"))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format()) + 
  ggplot2::scale_x_continuous(labels = NULL, breaks = NULL) + 
  ggplot2::labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) 

## ---- fig.align='center', fig.height=5, fig.width=7---------------------------
scenario_density(pnl$pnl_u, fp_es2) + 
  ggplot2::scale_x_continuous(labels = scales::dollar_format(prefix = "U$ ")) +
  ggplot2::labs(title = "Scenario Analysis with Fully-Flexible Probabilities (FFP)", 
                subtitle = "Historical vs. Exponential Smoothing")

## -----------------------------------------------------------------------------
empirical_stats(pnl %>% dplyr::select(pnl_u), fp_es1)

