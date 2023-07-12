library("pals")
library("ape")
library("BAMMtools")
library("phytools")
here::i_am("diversification/R/fossilbamm_results_noanc.R")

tree.pinnip <- read.tree(here::here("diversification/data/pinniped_median_noanc.tre"))

mcmcout <- read.csv(here::here("diversification/analyses/main_analysis/mcmc_out_pinnipedia_noanc.txt"), header = TRUE)
mcmcout <- mcmcout[floor(0.1 * nrow(mcmcout)):nrow(mcmcout),]
coda::effectiveSize(mcmcout$logLik)
coda::effectiveSize(mcmcout$N_shifts)

edata.pinnip <- getEventData(tree.pinnip, here::here("diversification/analyses/main_analysis/event_data_pinnipedia_noanc.txt"), burnin = 0.1)

post_probs <- table(mcmcout$N_shifts) / nrow(mcmcout)

summary(edata.pinnip)

png(here::here("diversification/outputs/figs/main_analysis/phylo_rates_pinnipedia_noanc_speciation.png"), width = 9, height = 9, units = "in", res = 300)
plot.bammdata(edata.pinnip, lwd = 2, legend = TRUE, xlim = c(0, 45), pal = parula(20)[-c(18:20)]);axisPhylo()
tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(edata.pinnip)
dev.off()

png(here::here("diversification/outputs/figs/main_analysis/phylo_rates_pinnipedia_noanc_extinction.png"), width = 9, height = 9, units = "in", res = 300)
plot.bammdata(edata.pinnip, lwd = 2, legend = TRUE, xlim = c(0, 45), pal = parula(20)[-c(18:20)], spex = "e");axisPhylo()
tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(edata.pinnip)
dev.off()

png(here::here("diversification/outputs/figs/main_analysis/phylo_rates_pinnipedia_noanc_netdiv.png"), width = 9, height = 9, units = "in", res = 300)
plot.bammdata(edata.pinnip, lwd = 2, legend = TRUE, xlim = c(0, 45), pal = parula(20)[-c(18:20)], spex = "netdiv");axisPhylo()
tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(edata.pinnip)
dev.off()

css <- credibleShiftSet(edata.pinnip, expectedNumberOfShifts = 1)
summary(css)

png(here::here("diversification/outputs/figs/main_analysis/css_pinnipedia_noanc.png"), width = 9, height = 9, units = "in", res = 300)
plot.credibleshiftset(css, lwd = 2, xlim = c(0, 45));axisPhylo()#;tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
dev.off()

best <- getBestShiftConfiguration(edata.pinnip, expectedNumberOfShifts = 1)
plot.bammdata(best, lwd = 2, legend = TRUE, xlim = c(0, 45));axisPhylo();tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(best, cex = 2.5)

#png(here::here("diversification/outputs/figs/main_analysis/rtt_speciation_pinnipedia_noanc.png"), width = 9, height = 3, units = "in", res = 300)
pdf(here::here("diversification/outputs/figs/main_analysis/rtt_speciation_pinnipedia_noanc_panel.pdf"), width = 9, height = 4)
plotRateThroughTime(edata.pinnip, ratetype = "speciation", intervalCol = "#CC6677", avgCol = "#CC6677", ylim = c(-0.1, 0.5))
abline(h = 0, lwd = 2, lty = 2, col = "grey")
dev.off()

#png(here::here("diversification/outputs/figs/main_analysis/rtt_extinction_pinnipedia_noanc.png"), width = 9, height = 3, units = "in", res = 300)
pdf(here::here("diversification/outputs/figs/main_analysis/rtt_extinction_pinnipedia_noanc_panel.pdf"), width = 9, height = 4)
plotRateThroughTime(edata.pinnip, ratetype = "extinction", intervalCol = "#CC6677", avgCol = "#CC6677", ylim = c(-0.1, 0.5))
abline(h = 0, lwd = 2, lty = 2, col = "grey")
dev.off()

#png(here::here("diversification/outputs/figs/main_analysis/rtt_netdiv_pinnipedia_noanc.png"), width = 9, height = 9, units = "in", res = 300)
pdf(here::here("diversification/outputs/figs/main_analysis/rtt_netdiv_pinnipedia_noanc_panel.pdf"), width = 9, height = 4)
plotRateThroughTime(edata.pinnip, ratetype = "netdiv", intervalCol = "#CC6677", avgCol = "#CC6677", ylim = c(-0.1, 0.5))
abline(h = 0, lwd = 2, lty = 2, col = "grey")
dev.off()


## Rates for Odobenidae

odob <- c("Odobenidae_gen_et_spec_indet_LACM_135920", "Odobenus_rosmarus")
odob.node <- getMRCA(tree.pinnip, odob)

## Rates for Otariidae

otar <- c("Eotaria_crypta", "Arctocephalus_australis")
otar.node <- getMRCA(tree.pinnip, otar)

## Rates for Phocinae

phoc <- c("Kawas_benegasorum", "Pusa_caspica")
phoc.node <- getMRCA(tree.pinnip, phoc)

## Rates for Monachinae

monac <- c("Frisiphoca_affine", "Hydrurga_leptonyx")
monac.node <- getMRCA(tree.pinnip, monac)

## Rates for Desmatophocidae

desmat <- c("Desmatophocidae_indet_USNM_335445", "Allodesmus_demerei")
desmat.node <- getMRCA(tree.pinnip, desmat)

## Rates for Stem
stem <- c("Odobenidae_gen_et_spec_indet_LACM_135920", "Hydrurga_leptonyx")
stem.node <- getMRCA(tree.pinnip, stem)

col.plots <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499")

#png(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_speciation_pinnipedia_per_family.png"), width = 9, height = 9, units = "in", res = 300)
pdf(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_speciation_pinnipedia_per_family.pdf"), width = 6, height = 4)
plotRateThroughTime(edata.pinnip, intervalCol = NULL, avgCol = NULL, ylim = c(-0.1, 0.4), lwd = 5)
plotRateThroughTime(edata.pinnip, node = stem.node, nodetype = "exclude", intervalCol = NULL, avgCol = col.plots[1], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = odob.node, nodetype = "include", intervalCol = NULL, avgCol = col.plots[2], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = otar.node, nodetype = "include", intervalCol = NULL, avgCol = col.plots[3], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = desmat.node, nodetype = "include", intervalCol = NULL, avgCol = col.plots[4], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = phoc.node, nodetype = "include", intervalCol = NULL, avgCol = col.plots[5], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = monac.node, nodetype = "include", intervalCol = NULL, avgCol = col.plots[6], add = TRUE, lwd = 4)
abline(h = 0, lwd = 2, lty = 2, col = "grey")
legend("topleft", legend = c("Stem", "Odobenidae", "Otariidae", "Desmatophocidae", "Phocinae", "Monachinae"), lwd = 4, col = col.plots[c(1, 2, 3, 4, 5, 6)])
dev.off()

#png(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_extinction_pinnipedia_per_family.png"), width = 9, height = 9, units = "in", res = 300)
pdf(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_extinction_pinnipedia_per_family.pdf"), width = 6, height = 4)
plotRateThroughTime(edata.pinnip, intervalCol = NULL, avgCol = NULL, ylim = c(-0.1, 0.4), lwd = 4)
plotRateThroughTime(edata.pinnip, node = stem.node, ratetype = "extinction", nodetype = "exclude", intervalCol = NULL, avgCol = col.plots[1], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = odob.node, ratetype = "extinction", nodetype = "include", intervalCol = NULL, avgCol = col.plots[2], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = otar.node, ratetype = "extinction", nodetype = "include", intervalCol = NULL, avgCol = col.plots[3], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = desmat.node, ratetype = "extinction", nodetype = "include", intervalCol = NULL, avgCol = col.plots[4], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = phoc.node, ratetype = "extinction", nodetype = "include", intervalCol = NULL, avgCol = col.plots[5], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = monac.node, ratetype = "extinction", nodetype = "include", intervalCol = NULL, avgCol = col.plots[6], add = TRUE, lwd = 4)
legend("topleft", legend = c("Stem", "Odobenidae", "Otariidae", "Desmatophocidae", "Phocinae", "Monachinae"), lwd = 4, col = col.plots[c(1, 2, 3, 4, 5, 6)])
abline(h = 0, lwd = 2, lty = 2, col = "grey")
dev.off()

#png(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_netdiv_pinnipedia_per_family.png"), width = 9, height = 9, units = "in", res = 300)
pdf(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_netdiv_pinnipedia_per_family.pdf"), width = 6, height = 4)
plotRateThroughTime(edata.pinnip, intervalCol = NULL, avgCol = NULL, ylim = c(-0.1, 0.4), lwd = 4)
plotRateThroughTime(edata.pinnip, node = stem.node, ratetype = "netdiv", nodetype = "exclude", intervalCol = NULL, avgCol = col.plots[1], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = odob.node, ratetype = "netdiv", nodetype = "include", intervalCol = NULL, avgCol = col.plots[2], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = otar.node, ratetype = "netdiv", nodetype = "include", intervalCol = NULL, avgCol = col.plots[3], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = desmat.node, ratetype = "netdiv", nodetype = "include", intervalCol = NULL, avgCol = col.plots[4], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = phoc.node, ratetype = "netdiv", nodetype = "include", intervalCol = NULL, avgCol = col.plots[5], add = TRUE, lwd = 4)
plotRateThroughTime(edata.pinnip, node = monac.node, ratetype = "netdiv", nodetype = "include", intervalCol = NULL, avgCol = col.plots[6], add = TRUE, lwd = 4)
legend("topleft", legend = c("Stem", "Odobenidae", "Otariidae", "Desmatophocidae", "Phocinae", "Monachinae"), lwd = 4, col = col.plots[c(1, 2, 3, 4, 5, 6)])
abline(h = 0, lwd = 2, lty = 2, col = "grey")
dev.off()
