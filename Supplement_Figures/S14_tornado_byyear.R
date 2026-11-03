## Figure S14.1: the one-at-a-time sensitivity analysis (Figure 6) split by year.
## Reads Output/Table4_Sensitivity_OAT.csv, so it uses the same values as Table 4 and Figure 6.
## Run by make_supplement_figures.R, with the working directory set to the package folder.
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
d <- read.csv("Output/Table4_Sensitivity_OAT.csv", stringsAsFactors = FALSE)
gv <- function(run, col) d[d$run == run, col]

spec <- list(
  list(nm="Energy floor  E<sub>D</sub>",     A="4b", B="4a"),
  list(nm="Gonad allocation  E<sub>G</sub>", A="5a", B="5b"),
  list(nm="Pejus temperature",               A="2a", B="2b"),
  list(nm="Encounter velocity",              A="1",  B="0"),
  list(nm="Initial energy  E<sub>I</sub>",   A="6a", B="6b"),
  list(nm="Spawning cost  E<sub>S</sub>",    A="7b", B="7a"),
  list(nm="Bay/delta  E<sub>Mbd</sub>",      A="8",  B="0"),
  list(nm="Moving fraction",                 A="3a", B="3b"))
FITTED <- "Moving fraction"

## ---- factor order is taken from the POOLED spans, i.e. from Figure 6 -------------------
ord <- do.call(rbind, lapply(spec, function(s)
  data.frame(factor=s$nm, span=abs(gv(s$A,"overall") - gv(s$B,"overall")), stringsAsFactors=FALSE)))
ord <- ord[order(ord$span), ]                       # ascending -> widest at the top
lev <- ord$factor
stopifnot(nrow(ord) == 8, ord$factor[1] == FITTED)  # the fitted parameter is still shortest

YRS <- c("2022 (warm)"="y2022", "2023 (cool)"="y2023")
fac <- do.call(rbind, lapply(names(YRS), function(yl) {
  cl <- YRS[[yl]]
  z <- do.call(rbind, lapply(spec, function(s) {
    a <- gv(s$A, cl); b <- gv(s$B, cl)
    data.frame(Year=yl, factor=s$nm, lo=min(a,b), hi=max(a,b), stringsAsFactors=FALSE) }))
  z$base <- gv("0", cl); z }))
fac$span <- fac$hi - fac$lo
fac$y    <- match(fac$factor, lev)
fac$grp  <- ifelse(fac$factor == FITTED, "fitted", "measured")
bl <- unique(fac[, c("Year","base")])

h <- 0.32
rects <- do.call(rbind, lapply(seq_len(nrow(fac)), function(i) { r <- fac[i,]; o <- list()
  if (r$lo < r$base) o[[length(o)+1]] <- data.frame(Year=r$Year, y=r$y, xmin=r$lo,
                                                    xmax=min(r$base,r$hi), grp=r$grp, dir="lower")
  if (r$hi > r$base) o[[length(o)+1]] <- data.frame(Year=r$Year, y=r$y, xmin=max(r$base,r$lo),
                                                    xmax=r$hi, grp=r$grp, dir="raise")
  do.call(rbind, o) }))
rects$key <- paste(rects$grp, rects$dir, sep="_")

## Value labels only; the parameter levels are given in Figure 6 and Table 4.
labs <- rbind(
  data.frame(Year=fac$Year, y=fac$y, x=fac$lo, hj=1.15,
             lab=sprintf("%.1f%s", fac$lo, ifelse(fac$factor=="Pejus temperature","*",""))),
  data.frame(Year=fac$Year, y=fac$y, x=fac$hi, hj=-0.15, lab=sprintf("%.1f", fac$hi)))

pal <- c(measured_lower="#D55E00", measured_raise="#0072B2",
         fitted_lower  ="#E8A87C", fitted_raise  ="#7FB3D9")
p <- ggplot() +
  geom_vline(data=bl, aes(xintercept=base), linetype="dashed", colour="grey35", linewidth=0.5) +
  geom_rect(data=rects, aes(xmin=xmin, xmax=xmax, ymin=y-h, ymax=y+h, fill=key),
            colour="white", linewidth=0.3) +
  geom_text(data=labs, aes(x=x, y=y, label=lab, hjust=hj), size=2.6, colour="grey15") +
  geom_text(data=bl, aes(x=base, y=8.85, label=sprintf("baseline %.1f%%", base)),
            size=3.0, colour="grey25", fontface="bold") +
  facet_wrap(~Year, ncol=2) +
  scale_fill_manual(values=pal, breaks=c("measured_lower","measured_raise","fitted_lower"),
                    labels=c("lowers spawning success","raises spawning success",
                             "the moving fraction — the only fitted parameter"),
                    name=NULL) +
  scale_y_continuous(breaks=seq_along(lev), labels=lev,
                     expand=expansion(add=c(0.9, 1.15))) +
  scale_x_continuous(limits=c(-6, 54), breaks=seq(0, 50, 10), expand=expansion(mult=0)) +
  labs(x="Predicted spawning success (%)", y=NULL,
       caption="* at 1 SD below the central pejus estimate no fish is predicted to survive in either year; that bound is excluded by our own telemetry.") +
  theme_bw(base_size=11) +
  theme(axis.text.y=ggtext::element_markdown(size=9.5, colour="grey15"),
        axis.title.x=element_text(size=10.5, margin=margin(t=7)),
        strip.background=element_rect(fill="grey92", colour="grey30", linewidth=0.4),
        strip.text=element_text(size=10.5, face="bold", colour="grey15"),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="grey30", linewidth=0.4),
        panel.spacing=unit(9,"pt"),
        legend.position="top", legend.justification="left", legend.direction="horizontal",
        legend.key.size=unit(10,"pt"), legend.margin=margin(b=-4),
        legend.text=element_text(size=9),
        plot.caption=element_text(hjust=0, size=7.6, colour="grey30", margin=margin(t=6)),
        plot.margin=margin(6,10,6,6))

out <- "Supplement_Figures/FigS14-1_tornado_byyear"
ggsave(paste0(out,".jpeg"), p, width=10, height=5.4, dpi=300)
cat("wrote", paste0(out,".jpeg"), "\n\n")

## ---- what the figure is for: the ranking, per year ------------------------------------
for (yl in names(YRS)) {
  z <- fac[fac$Year==yl, c("factor","lo","hi","span","base")]
  z <- z[order(-z$span), ]
  z$rel <- round(z$span / z$base, 2)
  cat(yl, "\n"); print(z, row.names=FALSE); cat("\n")
}
r22 <- fac$factor[fac$Year=="2022 (warm)"][order(-fac$span[fac$Year=="2022 (warm)"])]
r23 <- fac$factor[fac$Year=="2023 (cool)"][order(-fac$span[fac$Year=="2023 (cool)"])]
cat("rank differences between years:", sum(r22 != r23), "of 8 positions\n")

cat("rank order 2022:", r22, "\nrank order 2023:", r23, "\n")
rel <- with(fac, tapply(span/base, list(factor, Year), sum))
cat("\nrelative span (span / that year's baseline):\n"); print(round(rel, 2))
ex <- rownames(rel)[rel[,1] <= rel[,2]]                # 2022 NOT more sensitive
cat("2022 not the more sensitive year for:", paste(ex, collapse=", "), "\n")
