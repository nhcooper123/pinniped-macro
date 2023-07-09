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

pdf(here::here("diversification/outputs/figs/main_analysis/phylo_rates_pinnipedia_noanc_speciation.pdf"), width = 9, height = 9)
plot.bammdata(edata.pinnip, lwd = 2, legend = TRUE, xlim = c(0, 45), pal = parula(20)[-c(18:20)]);axisPhylo()
tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(edata.pinnip)
dev.off()

pdf(here::here("diversification/outputs/figs/main_analysis/phylo_rates_pinnipedia_noanc_extinction.pdf"), width = 9, height = 9)
plot.bammdata(edata.pinnip, lwd = 2, legend = TRUE, xlim = c(0, 45), pal = parula(20)[-c(18:20)], spex = "e");axisPhylo()
tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(edata.pinnip)
dev.off()

pdf(here::here("diversification/outputs/figs/main_analysis/phylo_rates_pinnipedia_noanc_netdiv.pdf"), width = 9, height = 9)
plot.bammdata(edata.pinnip, lwd = 2, legend = TRUE, xlim = c(0, 45), pal = parula(20)[-c(18:20)], spex = "netdiv");axisPhylo()
tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(edata.pinnip)
dev.off()

css <- credibleShiftSet(edata.pinnip, expectedNumberOfShifts = 1)
summary(css)

pdf(here::here("diversification/outputs/figs/main_analysis/css_pinnipedia_noanc.pdf"), width = 9, height = 9)
plot.credibleshiftset(css, lwd = 2, xlim = c(0, 45));axisPhylo()#;tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
dev.off()

best <- getBestShiftConfiguration(edata.pinnip, expectedNumberOfShifts = 1)
plot.bammdata(best, lwd = 2, legend = TRUE, xlim = c(0, 45));axisPhylo();tiplabels(edata.pinnip$tip.label, frame = "none", bg = NULL, adj = c(-0.1, 0.5), cex = 0.75)
addBAMMshifts(best, cex = 2.5)

pdf(here::here("diversification/outputs/figs/main_analysis/rtt_speciation_pinnipedia_noanc.pdf"), width = 9, height = 9)
plotRateThroughTime(edata.pinnip, ratetype = "speciation", intervalCol = "#CC6677", avgCol = "#CC6677")
dev.off()

pdf(here::here("diversification/outputs/figs/main_analysis/rtt_extinction_pinnipedia_noanc.pdf"), width = 9, height = 9)
plotRateThroughTime(edata.pinnip, ratetype = "extinction", intervalCol = "#CC6677", avgCol = "#CC6677")
dev.off()

pdf(here::here("diversification/outputs/figs/main_analysis/rtt_netdiv_pinnipedia_noanc.pdf"), width = 9, height = 9)
plotRateThroughTime(edata.pinnip, ratetype = "netdiv", intervalCol = "#CC6677", avgCol = "#CC6677")
abline(h = 0, lwd = 2, lty = 2, col = "grey")
dev.off()


## Rates for Odobenidae

odob <- c("Odobenidae_gen_et_spec_indet_LACM_135920", "Odobenus_rosmarus")
odob.node <- getMRCA(tree.pinnip, odob)

## Rates for Otariidae

otar <- c("Eotaria_crypta", "Arctocephalus_australis")
otar.node <- getMRCA(tree.pinnip, otar)

## Rates for Phocidae

phoc <- c("Devinophoca_emryi", "Hydrurga_leptonyx")
phoc.node <- getMRCA(tree.pinnip, phoc)

## Rates for Desmatophocidae

desmat <- c("Desmatophocidae_indet_USNM_335445", "Allodesmus_demerei")
desmat.node <- getMRCA(tree.pinnip, desmat)

## Rates for Stem
stem <- c("Arctocephalus_australis", "Hydrurga_leptonyx")
stem.node <- getMRCA(tree.pinnip, stem)

col.plots <- c("#332288", "#44AA99", "#DDCC77", "#882255", "#CC6677")

pdf(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_speciation_pinnipedia_per_family.pdf"), width = 9, height = 9)
plotRateThroughTime(edata.pinnip, intervalCol = NA, avgCol = col.plots[5], ylim = c(0, 0.4), lwd = 5)
plotRateThroughTime(edata.pinnip, node = odob.node, nodetype = "include", intervalCol = NA, avgCol = col.plots[1], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = otar.node, nodetype = "include", intervalCol = NA, avgCol = col.plots[4], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = phoc.node, nodetype = "include", intervalCol = NA, avgCol = col.plots[3], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = desmat.node, nodetype = "include", intervalCol = NA, avgCol = col.plots[2], add = TRUE, lwd = 5)
legend("topright", legend = c("Odobenidae", "Otariidae", "Phocidae", "Desmatophocidae", "Full"), lwd = 5, col = col.plots[c(1, 4, 3, 2, 5)])
dev.off()


pdf(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_extiction_pinnipedia_per_family.pdf"), width = 9, height = 9)
plotRateThroughTime(edata.pinnip, intervalCol = NA, ratetype = "extinction", avgCol = col.plots[5], ylim = c(0, 0.4), lwd = 5)
plotRateThroughTime(edata.pinnip, node = odob.node, ratetype = "extinction", nodetype = "include", intervalCol = NA, avgCol = col.plots[1], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = otar.node, ratetype = "extinction", nodetype = "include", intervalCol = NA, avgCol = col.plots[4], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = phoc.node, ratetype = "extinction", nodetype = "include", intervalCol = NA, avgCol = col.plots[3], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = desmat.node, ratetype = "extinction", nodetype = "include", intervalCol = NA, avgCol = col.plots[2], add = TRUE, lwd = 5)
legend("topright", legend = c("Odobenidae", "Otariidae", "Phocidae", "Desmatophocidae", "Full"), lwd = 5, col = col.plots[c(1, 4, 3, 2, 5)])
dev.off()


pdf(here::here("diversification/outputs/figs/main_analysis/rtt_noanc_netdiv_pinnipedia_per_family.pdf"), width = 9, height = 9)
plotRateThroughTime(edata.pinnip, intervalCol = NA, ratetype = "netdiv", avgCol = col.plots[5], ylim = c(-0.1, 0.4), lwd = 5)
plotRateThroughTime(edata.pinnip, node = odob.node, ratetype = "netdiv", nodetype = "include", intervalCol = NA, avgCol = col.plots[1], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = otar.node, ratetype = "netdiv", nodetype = "include", intervalCol = NA, avgCol = col.plots[4], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = phoc.node, ratetype = "netdiv", nodetype = "include", intervalCol = NA, avgCol = col.plots[3], add = TRUE, lwd = 5)
plotRateThroughTime(edata.pinnip, node = desmat.node, ratetype = "netdiv", nodetype = "include", intervalCol = NA, avgCol = col.plots[2], add = TRUE, lwd = 5)
abline(h = 0, lty = 2, col = "lightgrey")
legend("topright", legend = c("Odobenidae", "Otariidae", "Phocidae", "Desmatophocidae", "Full"), lwd = 5, col = col.plots[c(1, 4, 3, 2, 5)])
dev.off()
