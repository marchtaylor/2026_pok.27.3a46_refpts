# remove all existing objects ----------------
rm(list = ls())

# 1. load data ------------------------------------------------------
load(file = "data/objects.Rdata", verbose = T)
load(file = "model/settings.Rdata", verbose = T)
load(file = "model/condition_mse_objects.Rdata", verbose = T)


load(file = "model/newBRPs.Rdata", verbose = T)

df1 <- readRDS(file = "model/sum_wErr_noAR.rds")
df3 <- readRDS(file = "model/sum_wErr_wAR.rds")



# 2. summary figs ------------------------------------------------------------

L <- as.list(c("wErr_noAR", "wErr_wAR"))
RES <- lapply(L, function(x){
  
  if(x == "wErr_noAR"){
    df.x <- get("df1")
  }
  if(x == "wErr_wAR"){
    df.x <- get("df3")
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




# risk --------------------------------------------------------------------
# pcalc <- function(SSB, Blim){
#   # annual probability of SSB < Blim
#   annP <- apply(SSB, 1, function(x){mean(x < Blim)})
#   # Ps
#   P1 <- mean(annP)
#   P2 <- mean(apply(SSB, 2, function(x){max(x < Blim)}))
#   P3 <- max(annP)
#   return(list(P1 = P1, P2 = P2, P3 = P3))
# }

df_blim <- data.frame(blim = c(om@sr@params["b"]), iter = ac(seq(c(om@sr@params["b"]))))
df3a <- df3 |> left_join(df_blim)
df3a <- df3a |> mutate(unsafe = if_else(slot == "ssb", data < blim, NA))

tmp1 <- df3a |> filter(slot == "ssb", year %in% ((fy-ey):fy)) |> group_by(Ftarg, year) |> 
  summarise(annP = mean(data < blim)) |> ungroup() |> 
  group_by(Ftarg) |> summarise(prob1 = mean(annP), prob3 = max(annP))

tmp2 <- df3a |> filter(slot == "ssb", year %in% ((fy-ey):fy)) |> group_by(Ftarg, iter) |>
  summarise(prob = max(data < blim)) |> ungroup() |> 
  group_by(Ftarg) |> summarise(prob2 = mean(prob)) |> ungroup()
plot(tmp2)

tmp <- tmp1 |> left_join(tmp2)
tmp <- tmp |> pivot_longer(cols = -1, names_to = "risk_type", values_to = "value")
tmp


newdat <- data.frame(Ftarg = seq(0, max(df3$Ftarg), by = 0.001), risk_type = "prob3")
newdat$value <- with(subset(tmp, risk_type == "prob3"), approx(x = Ftarg, y = value, xout = newdat$Ftarg, rule = 2))$y
(Fpa <- newdat |> group_by(risk_type) |> summarise(Fpa = max(Ftarg[value < 0.05])) |> pull(Fpa))
Flabel <- paste("Fpa (Prob3)", sprintf("%.3f", Fpa), sep = " = ", collapse = "\n")


p <- ggplot(tmp) + aes(x = Ftarg, y = value, group = risk_type, colour = risk_type) + 
  geom_hline(yintercept = 0.05, lty = 3, col = 2) + 
  geom_vline(xintercept = Fpa, lty = 3) + 
  geom_line() +
  geom_point(size = 1) +
  # geom_line(data = newdat, lty = 2, color = 4) +
  annotate("label", x = Fpa, y = Inf, label = Flabel, hjust = 1.1, vjust = 1.1, size = 3) + 
  theme_bw()

png("output/risk_type~Fscan.png", width = 5, height = 3.5, units = "in", res = 400)
print(p)
dev.off()



# bootstrapped prob3 ------------------------------------------------------

dfboot <- df3a |> filter(slot == "ssb", year %in% ((fy-ey):fy)) 
newdat <- data.frame(Ftarg = seq(0,max(df3$Ftarg), by = 0.001), type = "p3")


# Split once outside loop
df_list <- split(dfboot, dfboot$iter)

nsamp <- 500
RES <- numeric(nsamp)
pb <- txtProgressBar(max = nsamp, style = 3)

for(i in seq_len(nsamp)) {

  samp_i <- sample(seq_len(it), size = it, replace = TRUE)

  # bind rows directly instead of join
  dfboot_i <- bind_rows(df_list[samp_i])

  tmp_i <- dfboot_i %>%
    group_by(Ftarg, year) %>%
    summarise(annP = mean(data < blim), .groups = "drop") %>%
    group_by(Ftarg) %>%
    summarise(prob3 = max(annP), .groups = "drop")

  newdat$prob3 <- approx(
    x = tmp_i$Ftarg,
    y = tmp_i$prob3,
    xout = newdat$Ftarg,
    rule = 2
  )$y

  RES[i] <- max(newdat$Ftarg[newdat$prob3 < 0.05])

  setTxtProgressBar(pb, i)
}

close(pb)
hist(RES, nclass = 20); abline(v = newBRPs$Fpa, col = 2)

sd(RES)/newBRPs$Fpa
sd(RES)/mean(RES)
