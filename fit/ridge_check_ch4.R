# Ridge insensitivity check (Appendix B).

this_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))])
BASE <- normalizePath(file.path(dirname(this_file), "..", "data", "networks"),
                      winslash = "/")
library(quadprog)

d0 <- file.path(BASE, "s25")
o <- readRDS(file.path(d0, "phy.rds"))
phy <- o$phy
tips <- read.csv(file.path(d0, "net_extant_tips.csv"))$phy_node
E <- phy$edge
len.t <- phy$edge.length.time
ret <- read.csv(file.path(d0, "net_retic_edges.csv"))
ret <- ret[ret$from != "", ]

# reference lengths x0: time-based lengths; hybrid-entry edges get the
# hybrid's blended distance to its tree parent (as in the simulation)
x0 <- len.t
nt0 <- read.csv(file.path(d0, "net_name_table.csv"))
node2sp0 <- setNames(as.integer(nt0$species), nt0$phy_node)
dm0 <- read.csv(file.path(d0, "dists2.csv"), row.names = 1, check.names = FALSE)
for (k in seq_len(nrow(ret))) {
  hyb <- ret$to[k]
  j <- which(E[, 2] == hyb)
  if (length(j) != 1) next
  tp <- E[j, 1]
  hyb.sp <- as.character(node2sp0[as.character(hyb)])
  tp.sp <- as.character(node2sp0[as.character(tp)])
  if (hyb.sp %in% rownames(dm0) && tp.sp %in% colnames(dm0)) {
    x0[j] <- dm0[hyb.sp, tp.sp]
  } else {
    inher <- read.csv(file.path(d0, "inher.csv"))
    g <- inher$inher[k]
    x0[j] <- g * ret$len[k] / (1 - g)
  }
}

  # prune to edges whose both endpoints keep at least one extant tip
tipset <- tips
ch <- split(E[, 2], E[, 1])
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
keep <- sapply(E[, 1], has) & sapply(E[, 2], has)
E.p <- E[keep, , drop = FALSE]
x0.p <- x0[keep]
m <- nrow(E.p)


nt <- read.csv(file.path(d0, "net_name_table.csv"))
sp2node <- setNames(nt$phy_node, nt$species)
node2sp <- setNames(as.integer(nt$species), nt$phy_node)
tips.sp <- as.character(node2sp[as.character(tips)])
dm <- read.csv(file.path(d0, "dists2.csv"), row.names = 1, check.names = FALSE)
Dobs <- as.matrix(dm[tips.sp, tips.sp])
pairs <- combn(seq_along(tips), 2)
d <- Dobs[cbind(pairs[1, ], pairs[2, ])]

  # path-edge incidence matrix A
root <- setdiff(E.p[, 1], E.p[, 2])
ch2 <- split(E.p[, 2], E.p[, 1])
A <- matrix(0, nrow = length(d), ncol = m)
for (col in seq_len(m)) {
  child <- E.p[col, 2]
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
  A[, col] <- idx[cbind(pairs[1, ], pairs[2, ])]
}

solve_qp <- function(ridge.scale) {
  Q <- crossprod(A)
  eps <- ridge.scale * max(1, max(diag(Q)))
  Dmat <- Q + diag(eps, m)
  dvec <- as.numeric(crossprod(A, d))
  Amat <- t(rbind(A, diag(m)))
  bvec <- c(d, rep(0, m))
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)
  as.numeric(A %*% sol$solution)
}

scales <- c(1e-10, 1e-8, 1e-6, 1e-4)
ref <- solve_qp(1e-10)
cat("ridge.scale  ||Dhat(s) - Dhat(1e-10)||\n")
for (s in scales) {
  Ds <- solve_qp(s)
  cat(sprintf("%-12g %.3e\n", s, sqrt(sum((Ds - ref)^2))))
}
