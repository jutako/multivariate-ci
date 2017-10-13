#' Greedy heuristics for computing MultiVariate Confidence Intervals (MVCI)
#' bottom-up
#' top-down

#' Greedy methods to find MVCIs
#' TODO: contains only bottom-up, not the more frequent top-down
#'
#'@description
#' A wrapper to call different methods with different initializations
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param k [1,1] integer, Number of observations (rows) to remove
#' @param l [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#' @param method char, Method to select, allowed: {'bottomup_byrow'}
#' @return A list describing the multivariate confidence band.
#'
#' @importFrom stats median
#'
#' @export
mvci_greedy <- function(dmat, k, l, method='bottomup_byrow'){
  N <- nrow(dmat)

  if (method == 'bottomup_byrow'){
    # add median as starting point
    dmat <- rbind(dmat, apply(dmat, 2, stats::median))
    cires <- findcb_bottomup_byrow(dmat, k, l, nrow(dmat))

    # remove starting point
    cires$Z <- cires$Z[1:N,]
    cires$upmask <- cires$upmask[1:N,]
    cires$downmask <- cires$downmask[1:N,]
  } else {
    stop(sprintf('Unknown method  \'%s\'.', method))
  }

  cires
}


#' Bottom-up method that adds full rows of data
#'
#'@description
#' Find MVCI using a bottom up greedy algorithm.
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param k [1,1] integer, Number of observations (rows) to remove
#' @param l [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#' @param ciSeedInds [1,J] integer, Row indices of rows that should be used as seed. Typically the median.
#'
#' @return A list describing the multivariate confidence band.
#'
#' @export
findcb_bottomup_byrow <- function(dmat, k, l, ciSeedInds){

  N <- nrow(dmat)
  M <- ncol(dmat)

  # data structure for storing Z
  Z <- matrix(F, nrow = N, ncol = M)
  Z[ciSeedInds,] <- T #initialize CI

  # data structure for storing current CI
  mvci <- createEnvelopeMasks(dmat, Z)

  # data structure for storing candidates
  candidateRowInds <- setdiff(1:nrow(dmat), ciSeedInds)

  # Loop 1:(N-k):
  # Evaluate candidates for adding
  # Exclude l most deviant points for each drow
  # Compute marginal cost by summin the rest

  # Select drow to add and update records
  for (b in 1:(N-k)){
    costArr <- vector(mode = "double",
                      length = length(candidateRowInds))
    for (ind in 1:length(candidateRowInds)){
      itmp <- points_outside_ci(candidateRowInds[ind], dmat, mvci)
      # note: itmp is ordered into ascending cost
      costArr[ind] <- sum(itmp$cost[1:(nrow(itmp)-l)])
    }
    rowToAdd <- candidateRowInds[which.min(costArr)]

    pointsDF <- points_outside_ci(rowToAdd, dmat, mvci)
    upto <- nrow(pointsDF)-l

    #browser()
    #cat(sprintf('upto: %d\n',upto))
    if (upto > 0){
      Z[rowToAdd, pointsDF$colind[1:upto]] <- T
      mvci <- updateEnvelopeMasks(mvci, rowToAdd, pointsDF[1:upto,])
      #browser()
      #plot.data.ci(dmat, Z, mvci$upmask, mvci$downmask, type = 'inside')
    } else {
      Z[rowToAdd,] <- T
    }

    candidateRowInds <- setdiff(candidateRowInds, rowToAdd)
  }

  list(Z = Z, upmask = mvci$upmask, downmask = mvci$downmask)
}



#' Identify rows that are inside/outside confidence band based on Z
#'
#' @description
#' Identify rows that lie completely within the mvci using mvci$Z
#'
#' @param Z [N,M] logical, Output of some mvci method, mvci$Z
#' @param type char {'inside', 'outside'}, inside -> rows completely TODO
#'
#' @return [N,1] logical, A logical vector indicating rows that match 'type'
#'
#' @export
getRows <- function(Z, type){
  switch(type,
         inside = (rowSums(Z) != 0),
         outside = (rowSums(Z) == 0) )
}


#' Initial creation of the envelope mask data structure
#'
#'@description
#' Create mvci envelope masks
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param Z [N,M] logical, Output of some mvci method, mvci$Z
#'
#' @return A list with mvci envelope masks.
#'
#' @keywords internal
createEnvelopeMasks <- function(dmat, Z){
  upmask = matrix(F, nrow = nrow(Z), ncol = ncol(Z))
  downmask = matrix(F, nrow = nrow(Z), ncol = ncol(Z))

  for (m in 1:dim(dmat)[[2]]){
    tmpinds <- which(Z[,m])
    maxind <- which.max(dmat[tmpinds,m])
    upmask[tmpinds[maxind], m] <- T

    tmpinds <- which(Z[,m])
    minind <- which.min(dmat[tmpinds,m])
    downmask[tmpinds[minind], m] <- T
  }

  out <- list(upmask = upmask,
              downmask = downmask)

  out
}

# updateEnvelopeMasks0 <- function(mvci, rowInd, upColInds, downColInds){
#   mvci$upmask[, upColInds] <- F #remove existing
#   mvci$upmask[rowInd, upColInds] <- T # add new
#   mvci$up <- mvci$upmask
#
#   mvci$downmask[, downColInds] <- F #remove existing
#   mvci$downmask[rowInd, downColInds] <- T # add new
#
#   mvci
# }


#' Update mvci envelope mask
#'
#' @description
#' Update the mvci envelope mask data structure by adding a new row
#'
#' @param mvci list, mvci result list
#' @param rowInd integer, Index of the row to add
#' @param pointDF data.frame, Data point data.frame
#'
#' @return An updated mvci list
#'
#' @keywords internal
updateEnvelopeMasks <- function(mvci, rowInd, pointDF){

  pdfs <- subset(pointDF, pointDF$isover)
  if (nrow(pdfs) > 0){
    mvci$upmask[, pdfs$colind] <- F #remove existing
    mvci$upmask[rowInd, pdfs$colind] <- T # add new
  }

  pdfs <- subset(pointDF, !pointDF$isover)
  if (nrow(pdfs) > 0){
    mvci$downmask[, pdfs$colind] <- F #remove existing
    mvci$downmask[rowInd, pdfs$colind] <- T # add new
  }

  mvci
}


#' Get confidence band (envelope) based on Z matrix
#'
#'@description
#' Not used much since envelope masks are easier to use
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param Z [N,M] logical, Output of some mvci method, mvci$Z
#'
#' @return A list with two vectors indicating values of mvci for each dimension
#'
#' @export
getEnvelope <- function(dmat, Z){
  out <- list(down = vector(mode = 'double', length = dim(dmat)[[2]]),
              up = vector(mode = 'double', length = dim(dmat)[[2]]) )

  for (m in 1:dim(dmat)[[2]]){
    out$down[m] <- min(dmat[Z[,m],m])
    out$up[m] <- max(dmat[Z[,m],m])
  }

  out
}


#' Get confidence band (envelope) based on row indices
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param rowinds [] integer, TODO
#'
#' @return A list with two vectors indicating values of mvci for each dimension
#'
#' @export
getEnvelopeRow <- function(dmat, rowinds){
  list( down = apply(dmat[rowinds,, drop = F], 2, min),
        up =   apply(dmat[rowinds,, drop = F], 2, max) )
}


# points_outside_ci <- function(drow, mvci){
#
#   over_inds <- which(drow > mvci$up)
#   over_costs <- drow[over_inds] - mvci$up[over_inds]
#
#   under_inds <- which(drow < mvci$down)
#   under_costs <- mvci$down[under_inds] - drow[under_inds]
#
#   list(inds = c(over_inds, under_inds),
#        isover = c( rep(T, length(over_inds)), rep(F, length(under_inds)) ),
#        costs = c(over_costs, under_costs))
# }


#' Find points on a row that lie outside confidence band
#'
#' @description
#' Report each point along with associated costs
#'
#' @param rowind integer, A row index
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param mvci list, mvci result list
#'
#' @return A data.frame with column indices, $isover and cost
#'
#' @keywords internal
points_outside_ci <- function(rowind, dmat, mvci){

  upArr <- dmat[mvci$upmask]
  over_inds <- which(dmat[rowind,] > upArr)
  over_costs <- dmat[rowind, over_inds] - upArr[over_inds]

  downArr <- dmat[mvci$downmask]
  under_inds <- which(dmat[rowind,] < downArr)
  under_costs <- downArr[under_inds] - dmat[rowind, under_inds]

  out <- data.frame(
            colind = c(over_inds, under_inds),
            isover = c( rep(T, length(over_inds)),
                        rep(F, length(under_inds)) ),
            cost = c(over_costs, under_costs) )
  out <- out[order(out$cost),]
  out
}




#'  Greedy top-down multivariate confidence interval (mvci)
#'
#' @description
#' Manuscript algorithm 1: Greedy
#'
#' Algorithm:
#'  do k times:
#'    1. solve k = 0, l > 0 to get CB
#'    2. decide which row to remove and remove it
#'  For the N-k row remaining rows, solve k = 0, l > 0 to get CB.
#'
#'  Output:
#'    Z [N, M] logical, indicates included points/rows
#'    upmask, downmask [N, M] logical, indicates positions of CB border, one TRUE per column
#'    row.inc.idx [1:(N-K)] integer, Indices of rows that are within CB
#'    row.exc.idx [1:K] integer,  Indices of rows that are outside CB
#'    L, K Input parameters for documentation
#'
#' TODO: Tries to process insane inputs such as L>M
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param K [1,1] integer, Number of observations (rows) to remove
#' @param L [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#' @param verbose logical, Verbosity of output
#'
#' @return A list describing the multivariate confidence band.
#'
#' @export
findcb_topdown <- function(dmat, K, L, verbose = F){ #(K(M^2logM + NlogN))

  N <- nrow(dmat) #observations / time
  M <- ncol(dmat) #variables


  row.inc.idx <- 1:N # indices of rows within confidence band

  # Greedily remove rows (i.e. observations) K -times
  if (K > 0){
    for (k in 1:K){#O(K)
      if (verbose){ cat(sprintf('Searching for row %d/%d to remove ...\n', k, K)) }

      cb <- remove_points_cb(dmat[row.inc.idx,], L) #O(M^2logM)
      costs <- collect_costs(dmat[row.inc.idx,], cb$d.order) #O(M)
      maxind <- which.max(costs$cost) #O(NlogN)
      row.inc.idx <- row.inc.idx[-maxind]

    }
  }

  # Compute CB using the L-criterion for the remaining rows
  cb <- remove_points_cb(dmat[row.inc.idx,], L)

  # Create outputs
  Z <- matrix(F, nrow = N, ncol = M)
  Z[row.inc.idx,] <- cb$Z

  downmask <- matrix(F, nrow = N, ncol = M)
  downmask[row.inc.idx,] <- cb$downmask

  upmask <- matrix(F, nrow = N, ncol = M)
  upmask[row.inc.idx,] <- cb$upmask

  list(Z = Z,
       upmask = upmask,
       downmask = downmask,
       row.inc.idx = row.inc.idx,
       row.exc.idx = setdiff(1:N, row.inc.idx),
       L = L,
       K = K)
}


#' Minimize envelope for k=0, l>0 (sub-step of findcb_topdown())
#'
#' @description
#' Solve problem k = 0, l > 0 by greedily removing points starting from the
#' whole data envelope.
#' To be used as a substep of the greedy top-down approach
#' Manuscript algorithm 2: FindEnvelope(X, I, L)
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param l [1,1] integer, Max number of outlier dimensions for a data row that is counted as within MVCI
#'
#' @return A list describing the multivariate confidence band.
#'
#' @importFrom utils tail
#'
#' @keywords internal
remove_points_cb <- function(dmat, l){ #O(M^2logM)

  N <- nrow(dmat)
  M <- ncol(dmat)

  # Points inside/outside selected set
  # Z is a N-by-M logical matrix indicating inclusion/exclusion of points
  # T -> keep point, F -> discard point
  Z <- matrix(T, nrow = N, ncol = M)

  # ordering of values
  # d.order is a M element list of N element integer vectors giving the ordering
  # permutations for each column of dmat.
  # Values are row indices [1,...,N].
  # Smallest element first!
  d.order <- rep(list(vector('integer',N)), M)
  for (m in 1:M){ #O(M)
    d.order[[m]] <- order(dmat[,m])
  }

  # record of points removed
  # d.removed.count is a N element integer vector keeping track of points
  # from each row of Z
  d.removed.count <- vector('integer', N)

  # record of gains attainable from CI border
  # initalize
  d.gains <- data.frame(gain = vector('numeric', 2*M),
                        col = vector('integer', 2*M),
                        islow = vector('logical', 2*M))
  # populate
  ind <- 1
  for (m in 1:M){ #O(M)
    d.gains$gain[ind] <- dmat[d.order[[m]][2], m] - dmat[d.order[[m]][1], m]
    d.gains$col[ind] <- m
    d.gains$islow[ind] <- T

    d.gains$gain[ind+1] <-  dmat[utils::tail(d.order[[m]],1), m] -
      dmat[utils::tail(d.order[[m]],2)[1], m]
    d.gains$col[ind+1] <- m
    d.gains$islow[ind+1] <- F

    ind <- ind + 2
  }

  # Do as long as there are potential candidates left
  points.left <- T
  while(points.left){ #O(M)

    if (nrow(d.gains) == 0){
      #exhausted all points -> terminate
      points.left <- F

    } else {
      # Search for the next point to remove
      tmpidx <- which.max(d.gains$gain) #O(MlogM)
      #previous step could be made faster using some fancy sorted list

      # get respective col and row
      ccol <- d.gains$col[tmpidx]
      if (d.gains$islow[tmpidx]){
        cislow <- T
        crow <- d.order[[ ccol ]][1]
      } else {
        cislow <- F
        crow <- utils::tail(d.order[[ ccol ]], 1)
      }

      # check if point can be removed
      if (d.removed.count[crow] >= l){
        # l criterion exceeded -> cannot remove from CB
        d.gains <- d.gains[-tmpidx,] #remove from search list

      } else {
        # all good -> proceed to removal

        # mark point as removed
        Z[crow, ccol] <- F
        d.removed.count[crow] <- d.removed.count[crow] + 1

        # update gains and CI border
        # d.gains$col and d.gains$islow remain intact
        if (length(d.order[[ccol]]) > 2){
          # on rare occasions almost all samples of one column migth have been rejected
          if (cislow){
            d.order[[ccol]] <- d.order[[ccol]][-1] #pop first
            d.gains$gain[tmpidx] <- dmat[d.order[[ccol]][2], ccol] -
              dmat[d.order[[ccol]][1], ccol]


          } else {
            d.order[[ccol]] <- d.order[[ccol]][-length(d.order[[ccol]])] #pop last
            d.gains$gain[tmpidx] <- dmat[utils::tail(d.order[[ccol]],1), ccol] -
              dmat[utils::tail(d.order[[ccol]],2)[1], ccol]
          }
          #if (is.na(d.gains$gain[tmpidx])){ browser() }
        } else {
          d.gains <- d.gains[-tmpidx,] #cannot update -> remove
        }
      }

    }

  }# of while

  # Create outputs
  downmask <- matrix(F, nrow = N, ncol = M)
  upmask <- matrix(F, nrow = N, ncol = M)
  for (m in 1:M){ #O(M)
    downmask[d.order[[m]][1], m] <- T
    upmask[utils::tail(d.order[[m]],1), m] <- T
  }

  list(Z = Z, upmask = upmask, downmask = downmask, d.order = d.order)
}


#' Collect boundary row costs (sub-step of findcb_topdown())
#'
#' @description
#' Collect costs attainable from rows that lie at the envelope boundary
#' Uses d.order structure to specify points at the border
#'
#' @param dmat [N,M] numeric, Data, vector valued M-dimensional observations on rows,
#' @param d.order TODO
#'
#' @return A data.frame with costs for each row of data
#'
#' @importFrom utils tail
#'
#' @keywords internal
collect_costs <- function(dmat, d.order){ #O(M)
  N <- nrow(dmat)
  M <- ncol(dmat)
  costs <- data.frame(row = 1:N,
                      cost = rep(0,N))
  for (m in 1:M){
    #low
    costs$cost[ d.order[[m]][1] ] <-  costs$cost[d.order[[m]][1]] +
                                      abs(dmat[d.order[[m]][2], m] -
                                          dmat[d.order[[m]][1], m])

    #up
    costs$cost[ utils::tail(d.order[[m]],1) ] <-  costs$cost[utils::tail(d.order[[m]],1)] +
                                              abs(dmat[utils::tail(d.order[[m]],1), m] -
                                                  dmat[utils::tail(d.order[[m]],2)[1], m])
  }

  costs
}
