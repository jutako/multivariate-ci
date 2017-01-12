# Creates the toy example in manuscript intro (Fig 1)

source("experiments/custom_plots.R")
require(ggplot2)
require(RColorBrewer)

## General options
fig.save.path <- OPTS$dir$pubfig.dir
alpha <- 0.1
Nbase = 100
Nbump = 99
M = 80

## Make M dimensional data
set.seed(16) #to make reproducible
dmat <- make.toy.data2(Nbase = Nbase, Nbump = Nbump, M = M)


## Compute different CB
k <- round(alpha * nrow(dmat))
l <- 25
#cb_naive <- cb_quantile(dmat, k, quantile.fun = central_quantile)
cb_naive <- cb_quantile(dmat, k)

cb_l0 <- findcb_topdown(dmat, K = k, L = 0 )
cbm0 <- cbr_extract_cb(dmat, cb_l0)

cb_l <- findcb_topdown(dmat, K = k, L = l )
cbml <- cbr_extract_cb(dmat, cb_l)


total.contained(dmat, cb_naive$cb[1,], cb_naive$cb[2,], 0)
total.contained(dmat, cbm0[1,], cbm0[2,], 0)
total.contained(dmat, cbml[1,], cbml[2,], 0)

total.contained(dmat, cb_naive$cb[1,], cb_naive$cb[2,], l)
total.contained(dmat, cbm0[1,], cbm0[2,], l)
total.contained(dmat, cbml[1,], cbml[2,], l)


## Create ggplottable datasets
attr(dmat, 'groups')[cb_l$row.inc.idx]

norm.row <- grep('supp.normal', attr(dmat, 'groups'))[1]
locout.row <- cb_l$row.inc.idx[98] # grep('local.outlier', attr(dmat, 'groups'))[2]
globout.row <- grep('global.outlier', attr(dmat, 'groups'))[1]

Nsubset <- 17
subset.inds <- c(norm.row, locout.row, globout.row,
                 sample(1:(nrow(dmat)), Nsubset, replace = F))
subset.group <- attr(dmat, 'groups')[subset.inds]

t.inds <- 5:M
pd <- make.ggplot.df(dmat[subset.inds, t.inds])

pdn <- make.ggplot.df(cb_naive$cb[, t.inds])
pdn$type <- 'naive'

pd0 <- make.ggplot.df(cbm0[, t.inds])
pd0$type <- 'L=0'
pd0$row <- pd0$row + 2

pdl <- make.ggplot.df(cbml[, t.inds])
pdl$type <- sprintf('L=%d',l)
pdl$row <- pdl$row + 4

pdcb <- rbind(pdn, pd0, pdl)
pdcb$type <- ordered(pdcb$type, levels = c('naive','L=0',sprintf('L=%d',l)))


# ## Plot
# lw <- 0.8
# lw2 <- 0.8
# colors <- RColorBrewer::brewer.pal(3, 'Dark2')
# p <- ggplot(data = pd, mapping = aes(x = variable, y = value, group = row))
#
# #p <- p + geom_line(alpha = 0.1)
# p <- p + geom_line(color = "#DBDBDB")
#
# # p <- p + geom_line(data = pdcb, aes(x = variable, y = value, group = row, color = type),
# #                    linetype = 'solid' , size = lw2)
#
# p <- p + geom_line(data = pdn, color = colors[[3]], linetype = 'solid', size = lw2)
# p <- p + geom_line(data = pd0, color = colors[[1]], linetype = 'solid', size = lw2)
# p <- p + geom_line(data = pdl, color = colors[[2]], linetype = 'solid', size = lw2)
#
# # p <- p + geom_line(data = pdn, color = 'red', linetype = 'solid' , size = lw2)
# # p <- p + geom_line(data = pd0, color = 'blue', linetype = 'solid', size = lw2)
# # p <- p + geom_line(data = pdl, color = 'green',, linetype = 'solid', size = lw2)
#
#
# p <- p + geom_line(data = subset(pd, row == 1), size = lw, linetype = 'dotdash') #normal
# p <- p + geom_line(data = subset(pd, row == 2), size = lw, linetype = 'solid') #local outlier
# p <- p + geom_line(data = subset(pd, row == 3), size = lw, linetype = 'dashed') #global outier
#
# p <- p + theme_light()
# p <- p + theme( axis.ticks = element_blank(),
#                 axis.text.x = element_blank(),
#                 axis.text.y = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.grid.major = element_blank(),
#                 legend.position="bottom",
#                 panel.border = element_blank(),
#                 plot.margin = unit(c(1,1,1,1.5), "lines"))
# p <- p + labs(x = NULL, y = NULL)
#
# # annotate confidence bands
# datatext <- data.frame(x = c('V1','V1','V1','V1','V1','V1','V28'),
#                       y = c(0.08, 0.05, 0.035, -0.03, -0.05, -0.07, -0.15),
#                       label = c('MWE','N','MCI','MCI','N','MWE'),
#                       stringsAsFactors = F)
# p <- p + geom_text(mapping = aes(x=x, y=y, label=label, group=NULL, hjust = "right"),
#                     nudge_x = -1,
#                    data =  datatext,
#                    size = 3)
#
# # some code to override clipping
# # From: http://stackoverflow.com/questions/12409960/ggplot2-annotate-outside-of-plot
# gt <- ggplot_gtable(ggplot_build(p))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# #grid.draw(gt)




# original x = 28, y = 41
x.idx <- 41
y.idx <- 56 #41, 48

#source("init.R")
gt <- motivating.plot(pdn, pd, pd0, pd1,
                      sprintf('V%d',t.inds[x.idx]),
                      sprintf('V%d',t.inds[y.idx]))

savename <- sprintf('toyplot_N%d_Nb%d_M%d.pdf', Nbase, Nbump, M)
ggsave(filename = file.path(fig.save.path, savename),
       plot = gt,
       units = 'mm',
       width = 100,
       height = 100)

inch.in.mm <- 25.4
w = 150 #mm
h = 150 #mm
dmats <- scale(dmat)
savename <- sprintf('toyplot_cross_N%d_Nb%d_M%d.pdf', Nbase, Nbump, M)
pdf(file.path(fig.save.path, savename),
    width = w/inch.in.mm, height = h/inch.in.mm)
cross.plot(dmats[, x.idx], dmats[, y.idx], 0.1*nrow(dmats)) #41, 48
dev.off()



#
# ## Save
# w = 100 #mm
# h = 100 #mm
# savename <- sprintf('toyplot_N%d_Nb%d_M%d.png', Nbase, Nbump, M)
# ggsave(filename = file.path(fig.save.path, savename),
#        plot = gt,
#        units = 'mm',
#        width = 100,
#        height = 100)
#
# savename <- sprintf('toyplot_N%d_Nb%d_M%d.pdf', Nbase, Nbump, M)
# ggsave(filename = file.path(fig.save.path, savename),
#        plot = gt,
#        units = 'mm',
#        width = 100,
#        height = 100)
#
#
#
# ## Save
# dmats <- scale(dmat)
# #dmats <- dmat # does not play nice with the custom plot, hence use scaling
#
# w = 150 #mm
# h = 150 #mm
#
# savename <- sprintf('toyplot_cross_N%d_Nb%d_M%d.png', Nbase, Nbump, M)
# png(file.path(fig.save.path, savename), units = 'mm', res = 300,
#     width = w, height = h)
# cross.plot(dmats[, x.idx], dmats[, y.idx], 0.1*nrow(dmats)) #41, 48
# dev.off()
#
# inch.in.mm <- 25.4
# savename <- sprintf('toyplot_cross_N%d_Nb%d_M%d.pdf', Nbase, Nbump, M)
# pdf(file.path(fig.save.path, savename),
#     width = w/inch.in.mm, height = h/inch.in.mm)
# cross.plot(dmats[, x.idx], dmats[, y.idx], 0.1*nrow(dmats)) #41, 48
# dev.off()


#pdf('./figs/alternative_crossScatter.pdf')
#dmats <- scale(dmat)

#source("init.R")
#cross.plot(dmats[, x.idx], dmats[, y.idx], 0.1*nrow(dmats)) #41, 48

#dev.off()



