# remove all existing objects ----------------
rm(list = ls())

# 1. load data ------------------------------------------------------
load(file = "data/objects.Rdata", verbose = T)
load(file = "model/settings.Rdata", verbose = T)
load(file = "model/condition_mse_objects.Rdata", verbose = T)


load(file = "model/newBRPs.Rdata", verbose = T)

df1 <- readRDS(file = "model/sum_wErr_noAR.rds")
df2 <- readRDS(file = "model/sum_wErr_wAR.rds")



# 2. summary figs ------------------------------------------------------------

L <- as.list(c("wErr_noAR", "wErr_wAR"))
RES <- lapply(L, function(x){
  
  if(x == "wErr_noAR"){
    df.x <- get("df1")
  }
  if(x == "wErr_wAR"){
    df.x <- get("df2")
  }
  
  newdat.x <- tibble(data.frame(Ftarg = seq(min(df.x$Ftarg), max(df.x$Ftarg), by = 1e-3), scen = x))
  yvars <- c("landings", "ssb")
  for(i in seq(yvars)){
    yvar <- yvars[i]
    agg <- aggregate(data ~ Ftarg, data = df.x, subset = slot == yvar & year %in% ((fy-ey):fy), FUN = quantile, prob = c(0.05, 0.5, 0.95))
    agg <- as.data.frame(cbind(agg[,1], agg[,-1]))
    names(agg) <- c("Ftarg", "q05", "q50", "q95")
    spl05 <- smooth.spline(x = agg$Ftarg, y = agg$q05)
    spl50 <- smooth.spline(x = agg$Ftarg, y = agg$q50)
    spl95 <- smooth.spline(x = agg$Ftarg, y = agg$q95)
    
    var <- paste(yvar, "q05", sep = "_")
    newdat.x <- newdat.x |> mutate(!!var := predict(spl05, x = Ftarg)$y)
    var <- paste(yvar, "q50", sep = "_")
    newdat.x <- newdat.x |> mutate(!!var := predict(spl50, x = Ftarg)$y)
    var <- paste(yvar, "q95", sep = "_")
    newdat.x <- newdat.x |> mutate(!!var := predict(spl95, x = Ftarg)$y)
  }
  
  return(newdat.x)
})

res <- tibble(do.call("rbind", RES))
res <- res |> mutate(
  scenario = case_when(
    scen %in% c("wErr_noAR") ~ "wErr_noAR (Fmsy, Bmsy)",
    scen %in% c("wErr_wAR") ~ "wErr_wAR (Fpa)"
  )
)

# res <- res |> mutate(scenario = factor(scenario, 
#   levels = c("wErr_noAR (Fmsy, Bmsy)", "noErr_noAR (Btrig)", "wErr_wAR (Fpa)", "wErr_noAR_noLag (test)")))

Frps <- as.data.frame(t(newBRPs[, c("Fmsy", "Fpa")]))
Frps$var <- factor(rownames(Frps), levels = c("Fmsy", "Fpa"))
names(Frps) <- c("value", "var")
Flabel <- paste(Frps$var, sprintf("%.3f", Frps$value), sep = " = ", collapse = "\n")

Brps <- as.data.frame(t(newBRPs[, c("Bmsy", "Btrigger", "Blim", "Bpa")]))
Brps$var <- factor(rownames(Brps), levels = c("Bmsy", "Btrigger", "Blim", "Bpa"))
names(Brps) <- c("value", "var")
Blabel <- paste(Brps$var, sprintf("%.0f", Brps$value), sep = " = ", collapse = "\n")


p1 <- ggplot(res) + 
  aes(x = Ftarg, y = landings_q50, ymin = landings_q05, ymax = landings_q95, 
    group = scenario, col = scenario, fill = scenario) +
  geom_ribbon(alpha = 0.2, lty = 0) + 
  geom_line() + 
  geom_vline(data = Frps, mapping = aes(xintercept = value, linetype = var), linewidth = 0.2, show.legend = F) +
  annotate("label", x = Inf, y = Inf, label = Flabel, hjust = 1.1, vjust = 1.1, size = 3) +
  lims(y = c(0,NA)) +
  labs(x = expression(bar(F)), y = "Landings [t]") + 
  theme_bw()
# p1

p2 <- ggplot(res) + 
  aes(x = Ftarg, y = ssb_q50, ymin = ssb_q05, ymax = ssb_q95, 
    group = scenario, col = scenario, fill = scenario) +
  geom_ribbon(alpha = 0.2, lty = 0) + 
  geom_line() + 
  geom_vline(data = Frps, mapping = aes(xintercept = value, linetype = var), linewidth = 0.2) +
  geom_hline(data = Brps, mapping = aes(yintercept = value, linetype = var), linewidth = 0.2) +
  # scale_y_continuous(expand = c(0, 0.05)) +
  annotate("label", x = Inf, y = Inf, label = Blabel, hjust = 1.1, vjust = 1.1, size = 3) +
  lims(y = c(0,NA), ) +
  labs(x = expression(bar(F)), y = "SSB [t]") + 
  theme_bw()
# p2

p <- (p1 / p2) + plot_layout(guides = 'collect', axes = "collect")
# p

png("output/yield+ssb~Fscan.png", width = 6, height = 5, units = "in", res = 400)
print(p)
dev.off()

