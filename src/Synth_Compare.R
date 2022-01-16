library(pacman)
p_load(data.table, patchwork, ggplot2, parallel, tictoc, Synth, CVXR)
theme_set(theme_minimal())

# %%
source("0CVXSynth.R")

# %% data prep
load('../data/ADH2015.RData')
ADH2015$country = factor(ADH2015$country)
pan = data.table(ADH2015)
pan[, treat := ifelse(country == "West Germany" & year >= 1990, 1, 0)]
treat_name = 'West Germany'
T0 = pan[country == "West Germany" & treat != 1, nunique(year)]
T1 = pan[country == "West Germany" & treat == 1, nunique(year)]

# %% # number of post-treatment periods reshape to wide
wide = pan[, .(country, year, gdp)] |> dcast(year ~ country, value.var = 'gdp')
setcolorder(wide, c('year', treat_name))
y_treat_pre = wide[1:T0, 2] |> as.matrix()
y_ctrl_pre  = wide[1:T0, -(1:2)] |> as.matrix()

# %% # Fit the SC algorithm
tic()
sc_fit = synth(X1 = y_treat_pre, X0 = y_ctrl_pre, Z1 = y_treat_pre, Z0 = y_ctrl_pre)
toc()

# %% synthetic control
tic()
ω_sc_MOSEK = sc_solve(y_treat_pre, y_ctrl_pre, solv = "MOSEK")
toc()
# %%
tic()
ω_sc_OSQP = sc_solve(y_treat_pre, y_ctrl_pre, solv = "OSQP")
toc()

# %%
tic()
ω_sc_ECOS = sc_solve(y_treat_pre, y_ctrl_pre, solv = "ECOS")
toc()
# %%
round(cbind(sc_fit$solution.w, ω_sc_MOSEK, ω_sc_OSQP, ω_sc_ECOS, ω_sc_SCS), 3) |>
  print()

# %%
RMSPE_CVXSynth = function(j, ypre, solv = "MOSEK"){
  # assign unit j as treated and other as as untreated, solve for weights
  y_j  = ypre[, j]; y_nj = ypre[, -j]
  wts = sc_solve(y_j, y_nj, solv = solv)
  if (is.null(wts)){ return(NULL) } else {
    # compute prediction error
    mse = (y_j -  y_nj  %*% wts)^2 |> mean() |> sqrt()
    return(mse)
  }
}


# %%
RMSPE_Synth = function(j, ypre){
  # assign unit j as treated and other as as untreated, solve for weights
  y_j  = ypre[, j] |> as.matrix(); y_nj = ypre[, -j] |> as.matrix()
  sc_fit = synth(X1 = y_j, X0 = y_nj, Z1 = y_j, Z0 = y_nj)
  if (is.null(wts)){ return(NULL) } else {
    # compute prediction error
    mse = (y_j -  y_nj  %*% sc_fit$solution.w)^2 |> mean() |> sqrt()
    return(mse)
  }
}
# %%
compare_RMSPE = function(x) c(RMSPE_CVXSynth(x, y_ctrl_pre), RMSPE_Synth(x, y_ctrl_pre))
compare_RMSPE(1)
RMSPE_comparison = lapply(1:16, compare_RMSPE)

# %%
RMSPE_comparison

# %%
