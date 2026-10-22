# QP branch-length fitting and kernel checks (Sections 4.2-4.3).
# For each network: build the reference lengths x0, solve min ||Ax - d||^2 s.t. Ax >= d, x >= 0 with a small ridge, and check the kernel formula (Corollary A.2).
library(quadprog)

this_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))])
BASE <- normalizePath(file.path(dirname(this_file), "..", "data", "networks"),
                      winslash = "/")
SEEDS <- c("s117","s120","s175","s25","s179","s187")

run_one <- function(seed) {
  d0 <- file.path(BASE, seed)
  o <- readRDS(file.path(d0, "phy.rds"))
  phy <- o$phy

  tips <- read.csv(file.path(d0, "net_extant_tips.csv"))$phy_node

  E <- phy$edge
  len.t <- phy$edge.length.time
  nE <- nrow(E)

  ret <- read.csv(file.path(d0, "net_retic_edges.csv"))
  ret <- ret[ret$from != "", ]
    x0 <- len.t
  # the tree edge entering a hybrid vertex gets the hybrid's blended distance
  # to its tree parent, from the simulation's dists2 bookkeeping
  nt0 <- read.csv(file.path(d0, 'net_name_table.csv'))
  node2sp0 <- setNames(as.integer(nt0$species), nt0$phy_node)
  dm0 <- read.csv(file.path(d0, 'dists2.csv'), row.names = 1, check.names = FALSE)
  for (k in seq_len(nrow(ret))) {
    hyb <- ret$to[k]
    j <- which(E[,2] == hyb)
    if (length(j) != 1) next
    tp <- E[j,1]
    hyb.sp <- as.character(node2sp0[as.character(hyb)])
    tp.sp <- as.character(node2sp0[as.character(tp)])
    if (hyb.sp %in% rownames(dm0) && tp.sp %in% colnames(dm0)) {
      x0[j] <- dm0[hyb.sp, tp.sp]
    } else {
      inher <- read.csv(file.path(d0, 'inher.csv'))
      g <- inher$inher[k]
      x0[j] <- g * ret$len[k] / (1 - g)
    }
  }

  tipset <- tips

  ch <- split(E[,2], E[,1])
  has_tip <- new.env()
  has <- function(v) {
    key <- as.character(v)
    if (exists(key, envir = has_tip)) return(get(key, envir = has_tip))
    if (v %in% tipset) { assign(key, TRUE, envir = has_tip); return(TRUE) }
    kids <- ch[[as.character(v)]]
    r <- FALSE
    if (!is.null(kids)) for (k in kids) if (has(k)) { r <- TRUE; break }
    assign(key, r, envir = has_tip)
    r
  }
  keep <- sapply(E[,1], has) & sapply(E[,2], has)
  E.p <- E[keep, , drop = FALSE]
  x0.p <- x0[keep]
  len.p <- len.t[keep]
  m <- nrow(E.p)

  dm <- read.csv(file.path(d0, "dists2.csv"), row.names = 1, check.names = FALSE)
  nt <- read.csv(file.path(d0, "net_name_table.csv"))
  sp2node <- setNames(nt$phy_node, nt$species)
  node2sp <- setNames(as.integer(nt$species), nt$phy_node)
  tips.sp <- as.character(node2sp[as.character(tips)])
  Dobs <- as.matrix(dm[tips.sp, tips.sp])

  pairs <- combn(seq_along(tips), 2)
  d <- Dobs[cbind(pairs[1,], pairs[2,])]
  # path-edge incidence matrix A

  root <- setdiff(E.p[,1], E.p[,2])

  par <- setNames(E.p[,1], E.p[,2])
  ch2 <- split(E.p[,2], E.p[,1])

  tree_dist <- function(ln) {

    t <- numeric()
    stack <- as.character(root)
    t[stack] <- 0
    while (length(stack)) {
      v <- stack[1]; stack <- stack[-1]
      for (w in ch2[[as.character(v)]]) {
        t[as.character(w)] <- t[as.character(v)] + ln[which(E.p[,2] == w)[1]]
        stack <- c(stack, as.character(w))
      }
    }
    t
  }
  Dtau <- tree_dist(len.p)
  Dref <- tree_dist(x0.p)

  A <- matrix(0, nrow = length(d), ncol = m)

  for (col in seq_len(m)) {
    e <- E.p[col, ]
    child <- e[2]

    sub <- c()
    stack <- child
    while (length(stack)) {
      v <- stack[1]; stack <- stack[-1]
      sub <- c(sub, v)
      for (w in ch2[[as.character(v)]]) stack <- c(stack, as.character(w))
    }
    sub <- as.integer(sub)

    ti <- tips %in% sub

    n_t <- length(tips)
    idx <- matrix(FALSE, n_t, n_t)
    idx[ti, !ti] <- TRUE
    idx[!ti, ti] <- TRUE

    A[, col] <- idx[cbind(pairs[1,], pairs[2,])]
  }
  # QP with ridge
  Q <- crossprod(A)
  eps <- 1e-10 * max(1, max(diag(Q)))
  Dmat <- Q + diag(eps, m)
  dvec <- as.numeric(crossprod(A, d))

  Amat <- t(rbind(A, diag(m)))
  bvec <- c(d, rep(0, m))
  sol <- tryCatch(solve.QP(Dmat, dvec, Amat, bvec, meq = 0),
                  error = function(e) NULL)
  if (is.null(sol)) return(list(seed = seed, error = "QP failed"))
  xhat <- sol$solution

  Ax0 <- as.numeric(A %*% x0.p)
  feas <- min(Ax0 - d)

  Dhat <- as.numeric(A %*% xhat)
  Efit <- sqrt(sum((Dhat - Dtau[cbind(pairs[1,], pairs[2,])])^2))
  Eref <- sqrt(sum((Ax0 - Dtau[cbind(pairs[1,], pairs[2,])])^2))
  # nullity check (Corollary A.2)
  rk <- qr(A)$rank
  nullity <- m - rk


  SS <- sum((Dhat - Ax0)^2)
  # hybrid-entry edges: those where x0 differs from the time-based length
  hyb.edges <- which(abs(x0.p - len.p) > 1e-9)
  pair.hyb <- (rowSums(A[, hyb.edges, drop = FALSE]) > 0)
  SS.hyb <- sum((Dhat - Ax0)[pair.hyb]^2)
  SS.oth <- SS - SS.hyb
  dev <- xhat - x0.p
  top.edges <- order(abs(dev), decreasing = TRUE)[1:5]
  top.is.hyb <- (top.edges %in% hyb.edges)
  write.csv(data.frame(from = E.p[,1], to = E.p[,2], x0 = x0.p, xhat = xhat, dev = dev),
            file.path(d0, "edge_dev_xhat_x0.csv"), row.names = FALSE)
  list(seed = seed, tips = length(tips), m = m, N = length(d),
       feas = feas, SS = SS, SS2 = sqrt(SS), SSref = sum((Ax0 - d)^2),
       SS.hyb = SS.hyb, SS.oth = SS.oth, n.hyb.pairs = sum(pair.hyb),
       n.hyb.edges = length(hyb.edges), top.dev = max(abs(dev)),
       top.edges = top.edges, top.is.hyb = top.is.hyb,
       nullity = nullity)
}

for (s in SEEDS) {
  r <- run_one(s)
  if (!is.null(r$error)) { cat(s, 'ERROR', r$error, '\n'); next }
  cat(sprintf('%s: SS=%.4f (hyb pairs %.4f = %.1f%% of SS, other %.4f) | hyb pairs %d/%d | hyb edges %d | top5 dev %.4f hyb? %s\n',
              r$seed, r$SS, r$SS.hyb, 100 * r$SS.hyb / r$SS, r$SS.oth,
              r$n.hyb.pairs, r$N, r$n.hyb.edges, r$top.dev,
              paste(r$top.is.hyb, collapse = ',')))
}
