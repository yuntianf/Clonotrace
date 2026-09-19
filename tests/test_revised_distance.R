# The revised clone distance on a small synthetic data set: runs the default clone_disance() path
# end to end and checks the properties the revision promises.
library(Clonotrace)
library(Matrix)
set.seed(1)

# 3 cell states in 5 dimensions; 12 clones, each drawing its cells from one or two states
n_cells <- 900
centres <- matrix(c(0, 0, 0, 0, 0,  6, 0, 0, 0, 0,  0, 6, 0, 0, 0), 3, byrow = TRUE)
state <- sample(3, n_cells, replace = TRUE)
emb <- centres[state, ] + matrix(rnorm(n_cells * 5), n_cells)
rownames(emb) <- paste0("cell", seq_len(n_cells))
colnames(emb) <- paste0("PC", 1:5)
clone_mix <- lapply(1:12, function(k) if (k <= 6) ((k - 1) %% 3) + 1 else c(((k - 1) %% 3) + 1, (k %% 3) + 1))
clone <- rep(NA_character_, n_cells)
for (k in 1:12) {
  pool <- which(state %in% clone_mix[[k]] & is.na(clone))
  clone[sample(pool, 25)] <- paste0("clone", k)
}
names(clone) <- rownames(emb)
clones <- paste0("clone", 1:12)
labelled <- which(!is.na(clone))
prob <- sparseMatrix(i = labelled, j = match(clone[labelled], clones), x = 1,
                     dims = c(n_cells, 12), dimnames = list(rownames(emb), clones))

out <- file.path(tempdir(), "revised_test")
dir.create(out, showWarnings = FALSE)
cell_clone <- data.frame(cell = names(clone), clone = clone)

# the default is the revised distance
res <- clone_disance(emb, prob, outpath = out, clone_labels = cell_clone, local_k = 3, verbose = FALSE)
s <- attr(res, "mean_weight")
stopifnot(all(c("group1", "group2", "dis", "dis_clonotrace") %in% colnames(res)),
          nrow(res) == choose(12, 2), all(is.finite(res$dis)), all(res$dis >= 0),
          length(s) == 1, s >= 0, s <= 1)

# revised = FALSE is the 1.1.0 behaviour, and the cached graph distance is reused unchanged
old <- clone_disance(emb, prob, outpath = out, revised = FALSE)
stopifnot(isTRUE(all.equal(as.numeric(old[, "dis"]), res$dis_clonotrace)))

# the long result matches the square result of revise_clone_distance() at the same weight
sq <- revise_clone_distance(old, emb, clone_labels = cell_clone, clone_names = clones,
                            local_k = 3, mean_weight = s, verbose = FALSE)
stopifnot(isSymmetric(unname(sq)), all(diag(sq) == 0),
          isTRUE(all.equal(sq[cbind(res$group1, res$group2)], res$dis)))

# the blend, written out: medians over all pairs, local scale = distance to the 3rd-nearest clone
D <- attr(sq, "clonotrace")
sig <- apply(D, 1, function(x) sort(x)[4])
L <- D / sqrt(outer(sig, sig))
M <- t(sapply(clones, function(k) colMeans(emb[which(clone == k), ])))
Dm <- as.matrix(dist(M))
u <- upper.tri(D)
manual <- (1 - s) * L / median(L[u]) + s * Dm / median(Dm[u])
stopifnot(isTRUE(all.equal(unname(sq[u]), unname(manual[u]))))

# the endpoints are the two distances alone
s0 <- revise_clone_distance(old, emb, clone_labels = cell_clone, clone_names = clones, local_k = 3,
                            mean_weight = 0, verbose = FALSE)
s1 <- revise_clone_distance(old, emb, clone_labels = cell_clone, clone_names = clones, local_k = 3,
                            mean_weight = 1, verbose = FALSE)
stopifnot(isTRUE(all.equal(unname(s0[u]), unname((L / median(L[u]))[u]))),
          isTRUE(all.equal(unname(s1[u]), unname((Dm / median(Dm[u]))[u]))))

# without observed labels the clone probabilities weight the cells; with one-hot probabilities
# that is the same clone mean
sp <- revise_clone_distance(old, emb, cell_clone_prob = prob, local_k = 3, mean_weight = s, verbose = FALSE)
stopifnot(isTRUE(all.equal(unname(sp), unname(sq), check.attributes = FALSE)))

# within-sample pairs: medians and the weight come from same-sample pairs only
smp <- setNames(rep(c("p1", "p2"), 6), clones)
ws <- revise_clone_distance(old, emb, clone_labels = cell_clone, clone_names = clones, sample = smp,
                            local_k = 3, verbose = FALSE)
stopifnot(attr(ws, "mean_weight") >= 0, attr(ws, "mean_weight") <= 1)

# supplied clusters replace the Leiden step
cl <- data.frame(r1 = state, row.names = rownames(emb))
wc <- revise_clone_distance(old, emb, clone_labels = cell_clone, clone_names = clones, clusters = cl,
                            local_k = 3, verbose = FALSE)
stopifnot(length(attr(wc, "picks")) == 1)

cat("revised distance: all checks passed (s =", s, ")\n")
