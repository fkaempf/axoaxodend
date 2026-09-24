# Axo-dendritic vs axo-axonic input, for any presynaptic x postsynaptic set.
#
#   source("axoaxo.R")   # from the repo root
#   res <- axoaxo(pre = "AN09B017f", post = "AVLP743m", name = "an09b017f_avlp743m")
#
# pre / post: bodyids and/or type names, mixed freely (types resolve via neuprint).
# For every postsynaptic body the arbor is split into axon and dendrite by synaptic
# flow centrality (Schneider-Mizell et al. 2016), exactly as in
# AVLP743m_connectomics/R/circuit/06_axon_dendrite.R: the edge carrying the most
# input-output paths is the cut, and the side holding most output sites is the axon.
# Each synapse from `pre` is then labelled axo-dendritic or axo-axonic by the
# domain its nearest skeleton node falls in.
#
# Writes to scenes/<name>/:
#   synapses.csv  one row per pre->post synapse, with domain
#   summary.csv   per post body x pre type: counts, % axo-axonic, segregation index
#   scene.json    Clio-NG state, one layer per neuron type and per pre type x domain; scene_url.txt holds the same state inline in a
#                 clio-ng URL, so it opens without hosting the JSON anywhere
#   <name>.webloc clickable link to that URL: double-click in Finder; opened
#                 automatically when run interactively (open = FALSE to skip)

suppressMessages({
  library(neuprintr); library(nat); library(igraph); library(nabor)
  library(dplyr); library(jsonlite)
})
# the folder this file lives in, when source()d; otherwise the working directory
AXO_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
source(file.path(AXO_DIR, "neuroglancer_scene.R"))   # make_mcns_scene

AXO_CACHE    <- file.path(AXO_DIR, "cache")
MCNS_DATASET <- "male-cns:v1.0"
MCNS_SERVER  <- Sys.getenv("NEUPRINT_HOST", "https://neuprint.janelia.org")
# bright on neuroglancer's black background
DOMAIN_COLS  <- c("axo-dendritic" = "#40C4FF", "axo-axonic" = "#FF5252")
# one colour per neuron type: warm for post types, cool for pre types
POST_PAL <- c("#ff0000", "#ff9100", "#ff4081", "#ffd600", "#d500f9", "#ff6e40")
PRE_PAL  <- c("#00e5ff", "#76ff03", "#2979ff", "#1de9b6", "#b388ff", "#c6ff00")
dir.create(AXO_CACHE, showWarnings = FALSE, recursive = TRUE)

# bodyids pass through; anything else is a type name
resolve_ids <- function(x, conn) {
  x <- as.character(x)
  is_id <- grepl("^[0-9]+$", x)
  types <- x[!is_id]
  hits <- if (length(types)) neuprint_search(
    sprintf("^(%s)$", paste(types, collapse = "|")), field = "type",
    dataset = MCNS_DATASET, conn = conn) else NULL
  miss <- setdiff(types, hits$type)
  if (length(miss)) stop("unknown type(s): ", paste(miss, collapse = ", "))
  unique(as.numeric(c(x[is_id], hits$bodyid)))
}

# bodyid -> type, "untyped" when neuprint has none
body_types <- function(ids, conn) {
  m <- neuprint_get_meta(ids, dataset = MCNS_DATASET, conn = conn)
  data.frame(bodyid = as.numeric(m$bodyid),
             type = ifelse(is.na(m$type) | m$type == "", "untyped", m$type))
}

conn_mcns <- function() {
  tryCatch(neuprint_login(server = MCNS_SERVER, dataset = MCNS_DATASET),
           error = function(e) neuprint_login(server = "https://neuprint-cns.janelia.org",
                                              dataset = MCNS_DATASET))
}

# skeleton + all synapses of one body, cached on disk
fetch_body <- function(bid, conn) {
  f <- file.path(AXO_CACHE, paste0(bid, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  n <- neuprint_read_neurons(bid, dataset = MCNS_DATASET, conn = conn)[[1]]
  s <- neuprint_get_synapses(bid, dataset = MCNS_DATASET, conn = conn)
  s <- s[, c("connector_id", "x", "y", "z", "bodyid", "partner", "prepost")]
  out <- list(neuron = n, syn = s)
  saveRDS(out, f)
  out
}

grp_entropy <- function(p) if (p <= 0 || p >= 1) 0 else -(p * log(p) + (1 - p) * log(1 - p))

# flow-centrality split, from 06_axon_dendrite.R. Returns the axon node set and the
# segregation index, or NULL when the body has too few sites to split.
split_body <- function(n, si, so, min_sites = 20) {
  if (nrow(si) < min_sites || nrow(so) < min_sites) return(NULL)
  xyz <- nat::xyzmatrix(n)
  si_node <- nabor::knn(xyz, as.matrix(si[, c("x", "y", "z")]), k = 1)$nn.idx[, 1]
  so_node <- nabor::knn(xyz, as.matrix(so[, c("x", "y", "z")]), k = 1)$nn.idx[, 1]

  gu <- as.undirected(as.ngraph(n, weights = TRUE), mode = "collapse",
                      edge.attr.comb = "first")
  root <- n$StartPoint
  bf <- igraph::bfs(gu, root = root, father = TRUE, order = TRUE, unreachable = FALSE)
  ord <- as.integer(bf$order); ord <- ord[!is.na(ord)]
  fat <- as.integer(bf$father)
  nv <- vcount(gu)

  cin  <- tabulate(si_node, nbins = nv)
  cout <- tabulate(so_node, nbins = nv)
  for (v in rev(ord)) {
    p <- fat[v]
    if (!is.na(p) && p > 0) { cin[p] <- cin[p] + cin[v]; cout[p] <- cout[p] + cout[v] }
  }
  tot_in <- nrow(si); tot_out <- nrow(so)

  cand <- setdiff(ord, root)
  flow <- cin[cand] * (tot_out - cout[cand]) + (tot_in - cin[cand]) * cout[cand]
  best <- cand[which.max(flow)]
  side <- as.integer(igraph::subcomponent(
    igraph::delete_edges(gu, igraph::get.edge.ids(gu, c(fat[best], best))),
    best, mode = "all"))
  axon_nodes <- if (cout[best] / tot_out >= 0.5) side else setdiff(seq_len(nv), side)

  dom_in  <- ifelse(si_node %in% axon_nodes, "axon", "dendrite")
  dom_out <- ifelse(so_node %in% axon_nodes, "axon", "dendrite")
  h_groups <- sum(vapply(c("axon", "dendrite"), function(dm) {
    ni <- sum(dom_in == dm); no <- sum(dom_out == dm)
    if (ni + no == 0) 0 else (ni + no) * grp_entropy(ni / (ni + no))
  }, numeric(1)))
  n_all <- tot_in + tot_out
  h_all <- n_all * grp_entropy(tot_in / n_all)

  list(dom_in = dom_in, seg_index = if (h_all > 0) 1 - h_groups / h_all else NA_real_)
}

point_layer <- function(df, name, colour) {
  dims <- list(x = list(8e-9, "m"), y = list(8e-9, "m"), z = list(8e-9, "m"))
  anns <- lapply(seq_len(nrow(df)), function(i)
    list(point = c(df$x[i], df$y[i], df$z[i]), type = "point", id = as.character(i)))
  list(type = "annotation",
       source = list(url = "local://annotations",
                     transform = list(outputDimensions = dims)),
       tool = "annotatePoint", annotations = anns, annotationColor = colour,
       shader = "void main() { setColor(defaultColor()); setPointMarkerSize(8.0); }",
       tab = "annotations", name = sprintf("%s (%d syn)", name, nrow(df)))
}

axoaxo <- function(pre, post, name, min_sites = 20, open = interactive()) {
  conn <- conn_mcns()
  pre_ids  <- resolve_ids(pre, conn)
  post_ids <- resolve_ids(post, conn)
  types    <- body_types(unique(c(pre_ids, post_ids)), conn)
  out_dir  <- file.path(AXO_DIR, "scenes", name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  per_body <- lapply(post_ids, function(b) {
    d  <- fetch_body(b, conn)
    si <- d$syn[d$syn$prepost == 1, ]
    so <- d$syn[d$syn$prepost == 0, ]
    sp <- split_body(d$neuron, si, so, min_sites)
    si$domain <- if (is.null(sp)) NA_character_ else
      ifelse(sp$dom_in == "axon", "axo-axonic", "axo-dendritic")
    si$seg_index <- if (is.null(sp)) NA_real_ else sp$seg_index
    si[si$partner %in% pre_ids, ]
  })
  syn <- bind_rows(per_body) %>%
    rename(post = bodyid, pre = partner) %>%
    left_join(types %>% rename(pre_type = type), by = c("pre" = "bodyid")) %>%
    left_join(types %>% rename(post_type = type), by = c("post" = "bodyid")) %>%
    select(connector_id, pre, pre_type, post, post_type, domain, seg_index, x, y, z)

  summ <- syn %>%
    group_by(post, post_type, pre_type) %>%
    summarise(n_syn = n(),
              n_axo_dendritic = sum(domain == "axo-dendritic"),
              n_axo_axonic    = sum(domain == "axo-axonic"),
              pct_axo_axonic  = 100 * n_axo_axonic / n_syn,
              seg_index       = first(seg_index), .groups = "drop") %>%
    arrange(post_type, post, desc(n_syn))
  unsplit <- setdiff(post_ids, syn$post[!is.na(syn$domain)])
  if (length(unsplit))
    message("no split (< ", min_sites, " input or output sites, or no synapses from pre): ",
            paste(unsplit, collapse = ", "))

  write.csv(syn,  file.path(out_dir, "synapses.csv"), row.names = FALSE)
  write.csv(summ, file.path(out_dir, "summary.csv"),  row.names = FALSE)

  # scene: one segmentation layer per neuron type (post warm, pre cool) and one
  # point layer per pre type x domain, so each type can be toggled on its own
  shown_pre <- intersect(pre_ids, syn$pre)
  ids <- c(post_ids, shown_pre)
  scene <- make_mcns_scene(ids, title = name)
  base  <- scene$layers[[2]]
  type_of <- function(b) types$type[match(as.numeric(b), types$bodyid)]
  seg_layers <- function(b, role, pal) {
    by_type <- split(as.character(b), type_of(b))
    cols <- rep_len(pal, length(by_type))
    Map(function(bb, colour) {
      l <- base
      l$segments <- as.list(bb)
      l$segmentQuery <- paste(bb, collapse = " ")
      l$segmentColors <- setNames(as.list(rep(colour, length(bb))), bb)
      l$name <- sprintf("%s %s (%d)", role, type_of(bb[1]), length(bb))
      l
    }, by_type, cols)
  }
  lab <- syn %>% filter(!is.na(domain))
  pt_layers <- list()
  for (tp in sort(unique(lab$pre_type)))
    for (dm in names(DOMAIN_COLS)) {
      d <- lab[lab$pre_type == tp & lab$domain == dm, ]
      if (nrow(d)) pt_layers <- c(pt_layers, list(
        point_layer(d, paste(tp, dm), unname(DOMAIN_COLS[dm]))))
    }
  scene$layers <- c(scene$layers[-2],
                    unname(seg_layers(post_ids, "post", POST_PAL)),
                    unname(seg_layers(shown_pre, "pre", PRE_PAL)),
                    pt_layers)
  if (nrow(syn)) {
    scene$position <- c(mean(syn$x), mean(syn$y), mean(syn$z))
    span <- max(diff(range(syn$x)), diff(range(syn$y)), diff(range(syn$z)))
    scene$projectionScale <- max(2500, round(span * 1.8))
  }
  js <- toJSON(scene, auto_unbox = TRUE, pretty = FALSE, digits = 8)
  write(js, file.path(out_dir, "scene.json"))
  url <- paste0("https://clio-ng.janelia.org/#!", URLencode(as.character(js), reserved = TRUE))
  write(url, file.path(out_dir, "scene_url.txt"))
  # the URL is ~60k characters, too long to copy around; a .webloc is a macOS
  # link file, so double-clicking it in Finder opens the scene in the browser
  link <- file.path(out_dir, paste0(name, ".webloc"))
  write(sprintf(paste0('<?xml version="1.0" encoding="UTF-8"?>\n',
                       '<plist version="1.0"><dict><key>URL</key><string>%s</string></dict></plist>'),
                gsub("&", "&amp;", url, fixed = TRUE)), link)
  if (open) system2("open", shQuote(link))

  cat(sprintf("%s: %d synapses, %d post bodies, %.1f%% axo-axonic overall\n", name,
              nrow(syn), length(post_ids), 100 * mean(syn$domain == "axo-axonic", na.rm = TRUE)))
  cat("written to", out_dir, "\n")
  cat("scene link:", link, "\n")
  invisible(list(synapses = syn, summary = summ, url = url))
}
