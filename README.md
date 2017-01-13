# multivariate-ci - multivariate confidence intervals for R

This repository contains R code to produce some of the figures in manuscript:

"Multivariate confidence intervals": Korpela J., Oikarinen E., Puolamäki K. and Ukkonen A., 2017, Proceedings of the 2017 SIAM International Conference on Data Mining (SDM 2017)

# Howto
To produce the figures use:
> source("run.R")

with the root of this repository as your working directory.

If some packages are missing, you can install them using:
> install.packages(c('grid','ggplot2','zoo','signal','RColorBrewer'))

# Contents

## Folders
Folder | Contents
------------ | -------------
experiments | Experiment related R code (usage examples)
rds | Real datasets in RDS format


## Files

File | Contents
------------ | -------------
mvci.R | Greedy algorithm for solving the MWE -problem
mvci_min_interval.R | Computation of minimum width intervals
mvci_tools.R | Naive bands, plotting, misc tools
experiments/mvci_data_loaders.R | Functions to create / load data
experiments/experiment_tools.R | Functions to create toy data, tools
experiments/custom_plots.R | Custom plotting functions
experiments/motivating_example.R | Code to reproduce manuscript Fig 1
experiments/stock_example.R | Code to reproduce manuscript Fig 5
run.R | A script to source() all needed files and run code that produces some of the figures into `./figs/`.
