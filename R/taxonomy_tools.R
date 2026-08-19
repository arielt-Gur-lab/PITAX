# ============================================================
# Taxonomic interpretation helpers
# Sanger Sequence Pipeline v2.10.2
# Score-aware, rank-aware identification logic
# ============================================================

as_num_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_taxon_text <- function(x) {
  x <- ifelse(is.na(x), "", trimws(as.character(x)))
  gsub("\\s+", " ", x)
}

# ---------------------------------------------------------------------------
# NCBI taxonomy resolution
# ---------------------------------------------------------------------------

resolve_missing_taxids_rentrez <- function(hits, sleep_seconds = 0.35) {
  if (is.null(hits) || !nrow(hits)) return(hits)
  if (!"taxid" %in% names(hits)) hits$taxid <- NA_integer_
  if (!"organism" %in% names(hits)) hits$organism <- ""

  missing <- which(is.na(suppressWarnings(as.numeric(hits$taxid))) & nzchar(clean_taxon_text(hits$organism)))
  if (!length(missing)) return(hits)

  cache <- new.env(parent = emptyenv())
  for (i in missing) {
    org <- clean_taxon_text(hits$organism[i])
    if (!nzchar(org)) next

    if (exists(org, envir = cache, inherits = FALSE)) {
      hits$taxid[i] <- get(org, envir = cache, inherits = FALSE)
      next
    }

    srch <- tryCatch(
      rentrez::entrez_search(
        db = "taxonomy",
        term = paste0('"', org, '"[Scientific Name]'),
        retmax = 2
      ),
      error = function(e) NULL
    )

    tx <- NA_integer_
    if (!is.null(srch) && length(srch$ids) == 1) {
      tx <- suppressWarnings(as.integer(srch$ids[1]))
    }
    assign(org, tx, envir = cache)
    hits$taxid[i] <- tx
    Sys.sleep(sleep_seconds)
  }
  hits
}

classification_for_taxid <- function(class_list, taxid, fallback_index = NULL) {
  if (is.null(class_list) || !length(class_list)) return(NULL)
  tx <- as.character(taxid)

  if (!is.null(names(class_list))) {
    idx <- which(names(class_list) == tx)
    if (length(idx)) return(class_list[[idx[1]]])
  }

  for (i in seq_along(class_list)) {
    cl <- class_list[[i]]
    if (is.data.frame(cl) && nrow(cl) && "id" %in% names(cl)) {
      ids <- as.character(cl$id)
      if (tx %in% ids) return(cl)
    }
  }

  if (!is.null(fallback_index) && fallback_index <= length(class_list)) {
    return(class_list[[fallback_index]])
  }
  NULL
}

rank_value <- function(classification_df, rank_names) {
  if (is.null(classification_df) || !is.data.frame(classification_df) || !nrow(classification_df)) return("")
  if (!all(c("name", "rank") %in% names(classification_df))) return("")
  rr <- tolower(trimws(as.character(classification_df$rank)))
  nn <- trimws(as.character(classification_df$name))
  for (rk in tolower(rank_names)) {
    idx <- which(rr == rk)
    if (length(idx)) return(nn[idx[length(idx)]])
  }
  ""
}

# ---------------------------------------------------------------------------
# Reference-quality annotation
# ---------------------------------------------------------------------------

is_unresolved_taxon_label <- function(x, rank = c("any", "genus", "species")) {
  rank <- match.arg(rank)
  x <- clean_taxon_text(x)
  bad <- !nzchar(x) |
    grepl("(^|\\s)(sp\\.?|cf\\.?|aff\\.?)(\\s|$)", x, ignore.case = TRUE) |
    grepl("uncultured|unclassified|unidentified|unknown|environmental sample|metagenome", x, ignore.case = TRUE)
  if (rank == "species") bad <- bad | lengths(strsplit(x, "\\s+")) < 2
  bad
}

reference_quality_label <- function(accession, title = "", organism = "") {
  accession <- clean_taxon_text(accession)
  title <- clean_taxon_text(title)
  organism <- clean_taxon_text(organism)

  type_pattern <- paste(
    c("type material", "type strain", "ex-type", "ex type", "holotype", "isotype",
      "lectotype", "neotype", "epitype", "paratype", "syntype"),
    collapse = "|"
  )

  if (grepl(type_pattern, title, ignore.case = TRUE)) return("Type material")

  # RefSeq accessions contain an underscore (e.g. NR_, NG_, NC_, NM_). For
  # fungal ITS, NR_ is especially relevant. This marks RefSeq status; it does
  # not claim that every RefSeq record is type-derived.
  if (grepl("^[A-Z]{2}_[0-9]+", accession)) return("RefSeq")

  if (is_unresolved_taxon_label(organism, "any")) return("Low / unresolved annotation")
  "Standard GenBank"
}

reference_quality_rank <- function(x) {
  map <- c(
    "Low / unresolved annotation" = 0,
    "Standard GenBank" = 1,
    "RefSeq" = 2,
    "Type material" = 3
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- 0
  out
}

best_reference_quality <- function(x) {
  if (!length(x)) return("Unavailable")
  r <- reference_quality_rank(x)
  if (!length(r) || all(!is.finite(r))) return("Unavailable")
  as.character(x[which.max(r)][1])
}

# ---------------------------------------------------------------------------
# Attach NCBI taxonomy + reference annotation
# ---------------------------------------------------------------------------

enrich_hits_with_taxonomy <- function(hits) {
  if (is.null(hits) || !nrow(hits)) return(hits)

  hits <- resolve_missing_taxids_rentrez(hits)
  ids <- unique(as.character(hits$taxid[!is.na(suppressWarnings(as.numeric(hits$taxid)))]))
  ids <- ids[nzchar(ids)]

  empty_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (nm in empty_cols) if (!nm %in% names(hits)) hits[[nm]] <- ""

  if (length(ids)) {
    classes <- tryCatch(
      taxize::classification(ids, db = "ncbi"),
      error = function(e) {
        x <- list()
        attr(x, "taxonomy_error") <- conditionMessage(e)
        x
      }
    )

    for (i in seq_len(nrow(hits))) {
      tx <- suppressWarnings(as.integer(hits$taxid[i]))
      if (is.na(tx)) next
      idx_id <- match(as.character(tx), ids)
      cl <- classification_for_taxid(classes, tx, fallback_index = idx_id)
      if (is.null(cl) || !is.data.frame(cl)) next

      hits$kingdom[i] <- rank_value(cl, c("kingdom", "superkingdom"))
      hits$phylum[i]  <- rank_value(cl, "phylum")
      hits$class[i]   <- rank_value(cl, "class")
      hits$order[i]   <- rank_value(cl, "order")
      hits$family[i]  <- rank_value(cl, "family")
      hits$genus[i]   <- rank_value(cl, "genus")
      hits$species[i] <- rank_value(cl, "species")
    }
  } else {
    classes <- list()
  }

  # Taxonomy services can occasionally fail while BLAST already supplied a
  # usable scientific name. Recover genus/species conservatively from the
  # NCBI organism label so a transient lineage lookup does not erase evidence.
  if ("organism" %in% names(hits)) {
    org <- clean_taxon_text(hits$organism)
    for (i in seq_len(nrow(hits))) {
      if (!nzchar(hits$genus[i]) && !is_unresolved_taxon_label(org[i], "any")) {
        tok <- strsplit(org[i], "\\s+")[[1]]
        if (length(tok)) hits$genus[i] <- tok[1]
      }
      if (!nzchar(hits$species[i]) && !is_unresolved_taxon_label(org[i], "species")) {
        tok <- strsplit(org[i], "\\s+")[[1]]
        if (length(tok) >= 2) hits$species[i] <- paste(tok[1:2], collapse = " ")
      }
    }
  }

  if (!"record_title" %in% names(hits)) hits$record_title <- ""
  if (!"accession" %in% names(hits)) hits$accession <- ""
  if (!"organism" %in% names(hits)) hits$organism <- ""
  hits$reference_quality <- mapply(
    reference_quality_label,
    hits$accession,
    hits$record_title,
    hits$organism,
    USE.NAMES = FALSE
  )

  attr(hits, "classification_list") <- classes
  attr(hits, "classification_ids") <- ids
  attr(hits, "taxonomy_error") <- attr(classes, "taxonomy_error")
  hits
}

# ---------------------------------------------------------------------------
# Locus context
# ---------------------------------------------------------------------------

locus_interpretation_note <- function(target) {
  target <- toupper(clean_taxon_text(target))
  if (grepl("^ITS$", target)) {
    return("ITS species-level resolution is taxon-dependent. Review difficult species complexes with an additional locus and/or curated type/reference evidence; methodological details are in Help / About.")
  }
  if (grepl("LSU", target)) {
    return("LSU can provide strong higher-level placement, but species-level resolution varies among groups. See Help / About for interpretation details.")
  }
  if (grepl("TEF1|EF1|RPB2|BETA|TUBULIN|CYP51|SDHB", target)) {
    return("Species-level discriminatory power for this locus is clade-dependent. See Help / About for interpretation details.")
  }
  "Species-level resolution depends on locus, clade and reference-database quality. See Help / About for interpretation details."
}

locus_benchmark <- function(target) {
  target_up <- toupper(clean_taxon_text(target))
  if (identical(target_up, "ITS")) {
    return(list(
      name = "ITS fungal broad benchmark",
      species_identity = 99.6,
      genus_identity = 94.3,
      source = "Vu et al. (2018), broad filamentous-fungi benchmark",
      universal = FALSE
    ))
  }
  list(
    name = "No universal identity cutoff applied",
    species_identity = NA_real_,
    genus_identity = NA_real_,
    source = "Direct molecular match + close alternatives + sequence/reference evidence",
    universal = FALSE
  )
}

# ---------------------------------------------------------------------------
# Score-aware helper functions
# ---------------------------------------------------------------------------

resolved_taxon_values <- function(x, rank = c("genus", "species")) {
  rank <- match.arg(rank)
  x <- clean_taxon_text(x)
  x[is_unresolved_taxon_label(x, rank)] <- ""
  x
}

order_hits_by_score <- function(hits) {
  bits <- as_num_safe(hits$bit_score)
  ev <- as_num_safe(hits$evalue)
  id <- as_num_safe(hits$identity_percent)
  cov <- as_num_safe(hits$query_coverage_percent)
  ord <- order(-bits, ev, -cov, -id, na.last = TRUE)
  hits[ord, , drop = FALSE]
}

separation_metrics <- function(candidate_hit, competitor_hit) {
  if (is.null(candidate_hit) || !nrow(candidate_hit)) {
    return(list(delta_bit = NA_real_, relative_percent = NA_real_, strength = "Unavailable", competitor = ""))
  }
  cb <- as_num_safe(candidate_hit$bit_score[1])
  if (is.null(competitor_hit) || !nrow(competitor_hit)) {
    return(list(delta_bit = NA_real_, relative_percent = NA_real_, strength = "No resolved competitor", competitor = ""))
  }
  bb <- as_num_safe(competitor_hit$bit_score[1])
  delta <- if (is.finite(cb) && is.finite(bb)) cb - bb else NA_real_
  rel <- if (is.finite(delta) && is.finite(cb) && cb > 0) 100 * delta / cb else NA_real_

  strength <- "Weak"
  if (is.finite(rel) && (rel >= 5 || (is.finite(delta) && delta >= 50))) {
    strength <- "Strong"
  } else if (is.finite(rel) && (rel >= 2 || (is.finite(delta) && delta >= 20))) {
    strength <- "Moderate"
  }

  comp_name <- ""
  if ("species" %in% names(competitor_hit) && nzchar(clean_taxon_text(competitor_hit$species[1]))) comp_name <- clean_taxon_text(competitor_hit$species[1])
  if (!nzchar(comp_name) && "genus" %in% names(competitor_hit)) comp_name <- clean_taxon_text(competitor_hit$genus[1])
  if (!nzchar(comp_name) && "organism" %in% names(competitor_hit)) comp_name <- clean_taxon_text(competitor_hit$organism[1])

  list(delta_bit = delta, relative_percent = rel, strength = strength, competitor = comp_name)
}

sequence_evidence_from_hit <- function(hit, target = "", rank = c("species", "genus")) {
  rank <- match.arg(rank)
  if (is.null(hit) || !nrow(hit)) return(list(level = "Unavailable", identity = NA_real_, coverage = NA_real_, benchmark = ""))
  identity <- as_num_safe(hit$identity_percent[1])
  coverage <- as_num_safe(hit$query_coverage_percent[1])
  bm <- locus_benchmark(target)

  level <- "Low"
  if (is.finite(identity) && is.finite(coverage)) {
    if (rank == "species" && is.finite(bm$species_identity)) {
      if (identity >= bm$species_identity && coverage >= 90) level <- "High"
      else if (identity >= max(bm$genus_identity, 97) && coverage >= 80) level <- "Moderate"
    } else if (rank == "genus" && is.finite(bm$genus_identity)) {
      if (identity >= bm$genus_identity && coverage >= 90) level <- "High"
      else if (identity >= bm$genus_identity && coverage >= 70) level <- "Moderate"
    } else {
      if (identity >= 99 && coverage >= 90) level <- "High"
      else if (identity >= 95 && coverage >= 80) level <- "Moderate"
    }
  }

  list(level = level, identity = identity, coverage = coverage, benchmark = bm$name)
}

# ---------------------------------------------------------------------------
# Evidence-first taxonomic interpretation (v2.14)
# ---------------------------------------------------------------------------

# Molecular matching is intentionally simple and conservative:
# 1) Partial high-identity alignments should not outrank near-full-length hits.
#    If any resolved hit covers >=90% of the query, only that coverage tier is
#    considered for the leading molecular match. Otherwise >=80% is preferred;
#    if neither tier exists, all usable hits remain eligible.
# 2) Inside the comparable-coverage tier, Identity is ranked first, then query
#    coverage, then Bit score and E-value as tie-breakers.
# This is an application decision rule, not a universal taxonomic threshold.
select_best_molecular_hit <- function(hits) {
  if (is.null(hits) || !nrow(hits)) return(NULL)
  h <- hits
  cov <- as_num_safe(h$query_coverage_percent)
  if (any(is.finite(cov) & cov >= 90)) {
    h <- h[is.finite(cov) & cov >= 90, , drop = FALSE]
  } else if (any(is.finite(cov) & cov >= 80)) {
    h <- h[is.finite(cov) & cov >= 80, , drop = FALSE]
  }
  if (!nrow(h)) return(NULL)

  id <- as_num_safe(h$identity_percent)
  cov <- as_num_safe(h$query_coverage_percent)
  bits <- as_num_safe(h$bit_score)
  ev <- as_num_safe(h$evalue)
  id_ord <- ifelse(is.finite(id), id, -Inf)
  cov_ord <- ifelse(is.finite(cov), cov, -Inf)
  bit_ord <- ifelse(is.finite(bits), bits, -Inf)
  ev_ord <- ifelse(is.finite(ev), ev, Inf)
  ord <- order(-id_ord, -cov_ord, -bit_ord, ev_ord, na.last = TRUE)
  h[ord[1], , drop = FALSE]
}

build_rank_evidence_profile <- function(hits, rank = c("species", "genus"),
                                        close_identity_pp = 0.5,
                                        close_coverage_pp = 2.0) {
  rank <- match.arg(rank)
  if (is.null(hits) || !nrow(hits) || !rank %in% names(hits)) return(data.frame())

  vals <- resolved_taxon_values(hits[[rank]], rank)
  keep <- nzchar(vals)
  if (!any(keep)) return(data.frame())
  h <- hits[keep, , drop = FALSE]
  vals <- vals[keep]
  groups <- split(seq_len(nrow(h)), vals)

  rows <- lapply(names(groups), function(taxon) {
    g <- h[groups[[taxon]], , drop = FALSE]
    rep_hit <- select_best_molecular_hit(g)
    if (is.null(rep_hit) || !nrow(rep_hit)) return(NULL)
    genus_name <- if (rank == "genus") taxon else {
      gv <- resolved_taxon_values(rep_hit$genus, "genus")
      if (length(gv) && nzchar(gv[1])) gv[1] else ""
    }
    data.frame(
      level = rank,
      taxon = taxon,
      genus = genus_name,
      best_accession = if ("accession" %in% names(rep_hit)) as.character(rep_hit$accession[1]) else "",
      best_identity_percent = as_num_safe(rep_hit$identity_percent[1]),
      best_query_coverage_percent = as_num_safe(rep_hit$query_coverage_percent[1]),
      best_bit_score = as_num_safe(rep_hit$bit_score[1]),
      best_evalue = as_num_safe(rep_hit$evalue[1]),
      accession_count = nrow(g),
      reference_support = best_reference_quality(g$reference_quality),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  cov <- as_num_safe(out$best_query_coverage_percent)
  tier <- ifelse(is.finite(cov) & cov >= 90, 2L, ifelse(is.finite(cov) & cov >= 80, 1L, 0L))
  if (any(tier == 2L)) eligible <- tier == 2L else if (any(tier == 1L)) eligible <- tier >= 1L else eligible <- rep(TRUE, nrow(out))

  id <- as_num_safe(out$best_identity_percent)
  bits <- as_num_safe(out$best_bit_score)
  ev <- as_num_safe(out$best_evalue)
  ord_metric <- order(
    !eligible,
    -ifelse(is.finite(id), id, -Inf),
    -ifelse(is.finite(cov), cov, -Inf),
    -ifelse(is.finite(bits), bits, -Inf),
    ifelse(is.finite(ev), ev, Inf),
    na.last = TRUE
  )
  out <- out[ord_metric, , drop = FALSE]
  rownames(out) <- NULL

  best_id <- as_num_safe(out$best_identity_percent[1])
  best_cov <- as_num_safe(out$best_query_coverage_percent[1])
  out$delta_identity_pp <- if (is.finite(best_id)) round(best_id - as_num_safe(out$best_identity_percent), 3) else NA_real_
  out$delta_coverage_pp <- if (is.finite(best_cov)) round(best_cov - as_num_safe(out$best_query_coverage_percent), 3) else NA_real_
  out$is_best_match <- seq_len(nrow(out)) == 1L
  out$is_close_alternative <- FALSE
  if (nrow(out) > 1 && is.finite(best_id) && is.finite(best_cov)) {
    alt_id <- as_num_safe(out$best_identity_percent)
    alt_cov <- as_num_safe(out$best_query_coverage_percent)
    out$is_close_alternative[-1] <- is.finite(alt_id[-1]) & is.finite(alt_cov[-1]) &
      alt_id[-1] >= (best_id - close_identity_pp) &
      alt_cov[-1] >= (best_cov - close_coverage_pp)
  }
  out$interpretation <- ifelse(out$is_best_match, "Best molecular match",
                               ifelse(out$is_close_alternative, "Close alternative", "Other database match"))
  out$database_share_percent <- round(100 * out$accession_count / sum(out$accession_count), 1)
  out$evidence_rank <- seq_len(nrow(out))
  out
}

profile_hit <- function(hits, profile_row, rank = c("species", "genus")) {
  rank <- match.arg(rank)
  if (is.null(profile_row) || !nrow(profile_row) || is.null(hits) || !nrow(hits)) return(NULL)
  vals <- resolved_taxon_values(hits[[rank]], rank)
  g <- hits[vals == as.character(profile_row$taxon[1]), , drop = FALSE]
  select_best_molecular_hit(g)
}

build_taxonomic_consensus <- function(enriched_hits, target = "", top_n = 25L) {
  if (is.null(enriched_hits) || !nrow(enriched_hits)) {
    return(list(summary = data.frame(), hits = data.frame(), counts = data.frame(), lca = data.frame(), score_cluster = data.frame()))
  }

  hits <- order_hits_by_score(enriched_hits)
  top_n <- max(1L, min(as.integer(top_n), nrow(hits)))
  hits <- hits[seq_len(top_n), , drop = FALSE]
  hits$analysis_rank <- seq_len(nrow(hits))

  if (!"reference_quality" %in% names(hits)) {
    if (!"record_title" %in% names(hits)) hits$record_title <- ""
    hits$reference_quality <- mapply(reference_quality_label, hits$accession, hits$record_title, hits$organism, USE.NAMES = FALSE)
  }

  bits <- as_num_safe(hits$bit_score)
  best_bit <- if (any(is.finite(bits))) max(bits, na.rm = TRUE) else NA_real_
  hits$delta_bit_from_best <- if (is.finite(best_bit)) round(best_bit - bits, 3) else NA_real_
  hits$relative_bit_score_percent <- if (is.finite(best_bit) && best_bit > 0) round(100 * bits / best_bit, 2) else NA_real_

  # Close-match windows are intentionally narrow and are used to detect whether
  # a different named taxon is practically indistinguishable from the best
  # molecular match. They reduce confidence; database frequency does not vote.
  close_identity_pp <- 0.5
  close_coverage_pp <- 2.0

  species_profile <- build_rank_evidence_profile(hits, "species", close_identity_pp, close_coverage_pp)
  genus_profile <- build_rank_evidence_profile(hits, "genus", close_identity_pp, close_coverage_pp)

  best_species_row <- if (nrow(species_profile)) species_profile[1, , drop = FALSE] else NULL
  best_genus_row <- if (nrow(genus_profile)) genus_profile[1, , drop = FALSE] else NULL
  best_species_hit <- profile_hit(hits, best_species_row, "species")
  best_genus_hit <- profile_hit(hits, best_genus_row, "genus")

  best_species <- if (!is.null(best_species_row)) as.character(best_species_row$taxon[1]) else ""
  best_genus <- if (!is.null(best_genus_row)) as.character(best_genus_row$taxon[1]) else ""

  close_species <- if (nrow(species_profile)) species_profile[species_profile$is_close_alternative, , drop = FALSE] else data.frame()
  close_genera <- if (nrow(genus_profile)) genus_profile[genus_profile$is_close_alternative, , drop = FALSE] else data.frame()
  closest_species_row <- if (nrow(close_species)) close_species[1, , drop = FALSE] else NULL
  closest_genus_row <- if (nrow(close_genera)) close_genera[1, , drop = FALSE] else NULL
  closest_species_hit <- profile_hit(hits, closest_species_row, "species")
  closest_genus_hit <- profile_hit(hits, closest_genus_row, "genus")

  best_species_genus <- if (!is.null(best_species_row) && "genus" %in% names(best_species_row)) as.character(best_species_row$genus[1]) else best_genus
  close_species_genera <- if (nrow(close_species)) unique(clean_taxon_text(close_species$genus[nzchar(clean_taxon_text(close_species$genus))])) else character()
  close_cross_genus_species <- if (nzchar(best_species_genus)) close_species_genera[close_species_genera != best_species_genus] else close_species_genera

  species_seq <- sequence_evidence_from_hit(best_species_hit, target, "species")
  genus_seq <- sequence_evidence_from_hit(best_genus_hit, target, "genus")
  species_ref <- if (!is.null(best_species_row)) as.character(best_species_row$reference_support[1]) else "Unavailable"
  genus_ref <- if (!is.null(best_genus_row)) as.character(best_genus_row$reference_support[1]) else "Unavailable"

  species_sep <- separation_metrics(best_species_hit, closest_species_hit)
  genus_sep <- separation_metrics(best_genus_hit, closest_genus_hit)

  species_discrimination <- if (!nzchar(best_species)) {
    "Unavailable"
  } else if (nrow(close_species) > 0) {
    "Poor"
  } else if (species_seq$level %in% c("High", "Moderate")) {
    "Good"
  } else {
    "Not established"
  }

  genus_discrimination <- if (!nzchar(best_genus)) {
    "Unavailable"
  } else if (nrow(close_genera) > 0 || length(close_cross_genus_species) > 0) {
    "Poor"
  } else {
    "Good"
  }

  locus_discrimination <- if (identical(genus_discrimination, "Poor")) {
    "Poor at genus and species level"
  } else if (identical(species_discrimination, "Poor")) {
    "Poor at species level"
  } else if (identical(species_discrimination, "Good")) {
    "Good at species level"
  } else {
    "Insufficient species-level evidence"
  }

  locus_flag <- ""
  if (species_seq$level == "High" && identical(species_discrimination, "Poor")) {
    locus_flag <- "High-quality sequence; multiple species are nearly indistinguishable at this locus."
  }
  if (genus_seq$level == "High" && identical(genus_discrimination, "Poor")) {
    locus_flag <- "High-quality sequence; close matches extend across different genera at this locus."
  }

  best_species_count <- if (!is.null(best_species_row)) as.integer(best_species_row$accession_count[1]) else 0L
  best_genus_count <- if (!is.null(best_genus_row)) as.integer(best_genus_row$accession_count[1]) else 0L

  species_promotable <- nzchar(best_species) &&
    species_seq$level == "High" &&
    nrow(close_species) == 0 &&
    reference_quality_rank(species_ref) >= 1

  genus_supported <- nzchar(best_genus) && identical(genus_discrimination, "Good")

  recommendation <- "Unresolved"
  rec_rank <- "unresolved"
  confidence <- "Low / review"
  taxonomic_support <- "Insufficient"
  sequence_evidence <- if (nzchar(best_genus)) genus_seq$level else species_seq$level
  selected_ref <- if (nzchar(best_genus)) genus_ref else species_ref

  if (species_promotable) {
    recommendation <- best_species
    rec_rank <- "species"
    taxonomic_support <- "High"
    sequence_evidence <- species_seq$level
    selected_ref <- species_ref
    confidence <- if (best_species_count >= 2 || reference_quality_rank(species_ref) >= 2) "High" else "Moderate"
  } else if (genus_supported) {
    recommendation <- best_genus
    rec_rank <- "genus"
    taxonomic_support <- "Very high"
    sequence_evidence <- genus_seq$level
    selected_ref <- genus_ref
    if (genus_seq$level == "High" && reference_quality_rank(genus_ref) >= 1) {
      confidence <- if (best_genus_count >= 2) "High" else "Moderate"
    } else if (genus_seq$level == "Moderate") {
      confidence <- "Moderate"
    } else {
      confidence <- "Low / review"
    }
  }

  # LCA is only a fallback when close molecular matches cross genus boundaries.
  # It is calculated from the genuinely close evidence, not the weak BLAST tail.
  classes <- attr(enriched_hits, "classification_list")
  lca_df <- data.frame(name = "", rank = "", id = "", stringsAsFactors = FALSE)
  if (identical(rec_rank, "unresolved") && identical(genus_discrimination, "Poor") && !is.null(classes) && length(classes)) {
    close_taxa <- c(best_genus, if (nrow(close_genera)) as.character(close_genera$taxon) else character())
    close_taxa <- unique(close_taxa[nzchar(close_taxa)])
    genus_vals <- resolved_taxon_values(hits$genus, "genus")
    close_hit_rows <- hits[genus_vals %in% close_taxa, , drop = FALSE]
    ids <- unique(as.character(close_hit_rows$taxid[!is.na(suppressWarnings(as.numeric(close_hit_rows$taxid)))]))
    ids <- ids[nzchar(ids)]
    if (length(ids) >= 2) {
      class_subset <- lapply(seq_along(ids), function(i) {
        classification_for_taxid(classes, ids[i], fallback_index = match(ids[i], attr(enriched_hits, "classification_ids")))
      })
      names(class_subset) <- ids
      class(class_subset) <- class(classes)
      lca_try <- tryCatch(taxize::lowest_common(ids, class_list = class_subset), error = function(e) NULL)
      if (!is.null(lca_try) && is.data.frame(lca_try) && nrow(lca_try)) {
        lca_df <- data.frame(
          name = as.character(lca_try$name[1]),
          rank = as.character(lca_try$rank[1]),
          id = if ("id" %in% names(lca_try)) as.character(lca_try$id[1]) else "",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (identical(rec_rank, "unresolved") && nzchar(lca_df$name[1])) {
    recommendation <- lca_df$name[1]
    rec_rank <- paste0("LCA: ", lca_df$rank[1])
    confidence <- "Low / review"
    taxonomic_support <- "Low / review"
  }

  # Mark accession-level evidence for plotting and audit. These flags describe
  # molecular closeness, not a majority vote.
  hits$is_best_molecular_match <- FALSE
  if (!is.null(best_species_hit) && nrow(best_species_hit) && "accession" %in% names(hits)) {
    hits$is_best_molecular_match <- as.character(hits$accession) == as.character(best_species_hit$accession[1])
  }
  species_vals <- resolved_taxon_values(hits$species, "species")
  close_species_names <- if (nrow(close_species)) as.character(close_species$taxon) else character()
  hits$is_close_species_alternative <- nzchar(species_vals) & species_vals %in% close_species_names
  genus_vals <- resolved_taxon_values(hits$genus, "genus")
  close_genus_names <- if (nrow(close_genera)) as.character(close_genera$taxon) else character()
  hits$is_close_genus_alternative <- nzchar(genus_vals) & genus_vals %in% close_genus_names
  hits$molecular_evidence_status <- ifelse(hits$is_best_molecular_match, "Best molecular match",
                                           ifelse(hits$is_close_species_alternative | hits$is_close_genus_alternative,
                                                  "Close alternative", "Other database match"))

  best_id <- if (!is.null(best_species_row)) as_num_safe(best_species_row$best_identity_percent[1]) else if (!is.null(best_genus_row)) as_num_safe(best_genus_row$best_identity_percent[1]) else NA_real_
  best_cov <- if (!is.null(best_species_row)) as_num_safe(best_species_row$best_query_coverage_percent[1]) else if (!is.null(best_genus_row)) as_num_safe(best_genus_row$best_query_coverage_percent[1]) else NA_real_
  best_acc <- if (!is.null(best_species_row)) as.character(best_species_row$best_accession[1]) else if (!is.null(best_genus_row)) as.character(best_genus_row$best_accession[1]) else ""
  best_match_name <- if (nzchar(best_species)) best_species else best_genus

  alt_name <- if (!is.null(closest_species_row)) as.character(closest_species_row$taxon[1]) else ""
  alt_acc <- if (!is.null(closest_species_row)) as.character(closest_species_row$best_accession[1]) else ""
  alt_id <- if (!is.null(closest_species_row)) as_num_safe(closest_species_row$best_identity_percent[1]) else NA_real_
  alt_cov <- if (!is.null(closest_species_row)) as_num_safe(closest_species_row$best_query_coverage_percent[1]) else NA_real_
  alt_count <- if (!is.null(closest_species_row)) as.integer(closest_species_row$accession_count[1]) else 0L

  species_conclusion <- if (rec_rank == "species") {
    paste0("Resolved to best molecular match: ", best_species)
  } else if (nrow(close_species) > 0) {
    "Unresolved: one or more species have nearly indistinguishable molecular matches"
  } else if (nzchar(best_species)) {
    "Not promoted to species: sequence/reference evidence is not strong enough"
  } else {
    "No resolved species-level match"
  }

  fmt_pct <- function(x, digits = 2) if (is.finite(x)) paste0(round(x, digits), "%") else "unavailable"
  database_note <- if (nzchar(best_species)) paste0("Database representation: ", best_species, " = ", best_species_count, " accession(s)",
                                                   if (nzchar(alt_name)) paste0("; ", alt_name, " = ", alt_count, " accession(s)") else "",
                                                   ". Accession counts are context only and do not vote on the identification.") else ""

  decision_reason <- if (rec_rank == "species") {
    paste0(
      "Best molecular match: ", best_species, " (", fmt_pct(best_id), " identity; ", fmt_pct(best_cov, 1), " query coverage). ",
      "No different resolved species falls inside the close-match window (within ", close_identity_pp, " identity percentage points and ", close_coverage_pp, " coverage points). ",
      database_note
    )
  } else if (rec_rank == "genus" && nrow(close_species) > 0) {
    paste0(
      "Best molecular match: ", best_match_name, " (", fmt_pct(best_id), " identity; ", fmt_pct(best_cov, 1), " query coverage). ",
      nrow(close_species), " close species alternative(s) remain, but the close evidence stays within genus ", best_genus, ". ",
      "The locus therefore supports the genus while species-level discrimination is poor. ", database_note
    )
  } else if (rec_rank == "genus") {
    paste0(
      "Best molecular match: ", best_match_name, " (", fmt_pct(best_id), " identity; ", fmt_pct(best_cov, 1), " query coverage). ",
      "No close alternative genus was detected, but species-level evidence was not strong enough for a conservative species recommendation. ", database_note
    )
  } else if (grepl("^LCA:", rec_rank)) {
    paste0(
      "Best molecular match: ", best_match_name, ". Close molecular alternatives extend across different genera, so genus-level identification is withheld. ",
      "The lowest common taxon among the close evidence is ", recommendation, "."
    )
  } else {
    paste0(
      "Best molecular match: ", ifelse(nzchar(best_match_name), best_match_name, "no resolved taxon"), ". ",
      "Close molecular alternatives and/or insufficient sequence evidence prevent a conservative genus/species recommendation."
    )
  }

  bm <- locus_benchmark(target)

  # Species evidence profile: one row per named species. Accession counts show
  # database representation only; they are not used as majority votes.
  counts <- species_profile
  if (nrow(counts)) {
    counts$close_identity_window_pp <- close_identity_pp
    counts$close_coverage_window_pp <- close_coverage_pp
  }

  summary <- data.frame(
    algorithm_version = "evidence-first-v2.14.0",
    recommended_identification = recommendation,
    recommended_level = rec_rank,
    confidence = confidence,
    taxonomic_support = taxonomic_support,
    sequence_evidence = sequence_evidence,
    reference_support = selected_ref,
    decision_reason = decision_reason,
    hits_used = nrow(hits),
    best_molecular_match = best_match_name,
    best_match_species = best_species,
    best_match_genus = best_genus,
    best_match_accession = best_acc,
    best_match_identity_percent = ifelse(is.finite(best_id), round(best_id, 3), NA_real_),
    best_match_query_coverage_percent = ifelse(is.finite(best_cov), round(best_cov, 2), NA_real_),
    best_match_reference_support = species_ref,
    closest_alternative_species = alt_name,
    closest_alternative_accession = alt_acc,
    closest_alternative_identity_percent = ifelse(is.finite(alt_id), round(alt_id, 3), NA_real_),
    closest_alternative_query_coverage_percent = ifelse(is.finite(alt_cov), round(alt_cov, 2), NA_real_),
    closest_alternative_delta_identity_pp = if (!is.null(closest_species_row)) round(as_num_safe(closest_species_row$delta_identity_pp[1]), 3) else NA_real_,
    closest_alternative_delta_coverage_pp = if (!is.null(closest_species_row)) round(as_num_safe(closest_species_row$delta_coverage_pp[1]), 3) else NA_real_,
    close_species_count = nrow(close_species),
    close_genus_count = nrow(close_genera),
    species_level_conclusion = species_conclusion,
    species_discrimination = species_discrimination,
    genus_discrimination = genus_discrimination,
    locus_discrimination = locus_discrimination,
    locus_flag = locus_flag,
    best_species_database_accessions = best_species_count,
    closest_species_database_accessions = alt_count,
    close_identity_window_pp = close_identity_pp,
    close_coverage_window_pp = close_coverage_pp,
    # Compatibility fields retained for older exports/projects; they no longer
    # drive the decision engine.
    candidate_genus = best_genus,
    species_candidate = best_species,
    species_candidate_confidence = if (rec_rank == "species") confidence else if (nzchar(best_species)) "Low / review" else "Insufficient",
    genus_reference_support = genus_ref,
    species_reference_support = species_ref,
    genus_best_competitor = if (!is.null(closest_genus_row)) as.character(closest_genus_row$taxon[1]) else "",
    genus_delta_bit = ifelse(is.finite(genus_sep$delta_bit), round(genus_sep$delta_bit, 2), NA_real_),
    genus_separation = genus_sep$strength,
    species_best_competitor = alt_name,
    species_delta_bit = ifelse(is.finite(species_sep$delta_bit), round(species_sep$delta_bit, 2), NA_real_),
    species_separation = species_sep$strength,
    candidate_identity_percent = ifelse(is.finite(best_id), round(best_id, 3), NA_real_),
    candidate_query_coverage_percent = ifelse(is.finite(best_cov), round(best_cov, 2), NA_real_),
    its_species_identity_benchmark = ifelse(is.finite(bm$species_identity), bm$species_identity, NA_real_),
    its_genus_identity_benchmark = ifelse(is.finite(bm$genus_identity), bm$genus_identity, NA_real_),
    benchmark_note = paste0(bm$name, "; ", bm$source, "; not a universal hard cutoff"),
    genus_resolved_hits = sum(nzchar(resolved_taxon_values(hits$genus, "genus"))),
    species_resolved_hits = sum(nzchar(resolved_taxon_values(hits$species, "species"))),
    lowest_common_taxon = lca_df$name[1],
    lowest_common_rank = lca_df$rank[1],
    locus_note = locus_interpretation_note(target),
    stringsAsFactors = FALSE
  )

  list(summary = summary, hits = hits, counts = counts, lca = lca_df, score_cluster = data.frame())
}


# ---------------------------------------------------------------------------
# Export taxonomic checkpoint
# ---------------------------------------------------------------------------

make_taxonomy_checkpoint_zip <- function(file, tax_summary, tax_hits, tax_counts, blast_hits = NULL) {
  temp_dir <- tempfile("sanger_taxonomy_")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  write.csv(tax_summary, file.path(temp_dir, "taxonomic_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(tax_hits, file.path(temp_dir, "taxonomic_hits_enriched.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(tax_counts, file.path(temp_dir, "species_evidence_profile.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  if (!is.null(blast_hits) && nrow(blast_hits)) {
    write.csv(blast_hits, file.path(temp_dir, "BLAST_hits_source.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  }

  app_version <- tryCatch(get("APP_VERSION", inherits = TRUE), error = function(e) "unknown")
  note <- c(
    paste0("Application version: ", app_version),
    paste0("Exported at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "v2.14 uses evidence-first accession-level BLAST interpretation.",
    "Best molecular match is selected from comparable-coverage hits by Identity first, then query coverage; Bit score and E-value are tie-breakers.",
    "Near-full-length hits (>=90% query coverage) are preferred over partial high-identity matches when available; >=80% is the next coverage tier.",
    "A different species is treated as a close alternative when it is within 0.5 identity percentage points and no more than 2 coverage points below the best molecular match.",
    "Close alternatives reduce the taxonomic level/confidence. If close species remain within one genus, genus identification can remain strong while species is unresolved.",
    "Accession counts describe database representation only and are not majority votes.",
    "Reference quality (type-material wording / RefSeq / standard GenBank / unresolved annotation) is reported separately.",
    "For ITS, 99.6% species and 94.3% genus identity values remain broad literature context, not universal hard cutoffs.",
    "A high-quality sequence with several nearly indistinguishable species is flagged as poor species-level discrimination by the locus."
  )
  writeLines(note, file.path(temp_dir, "README_taxonomic_interpretation.txt"))

  files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  zip::zipr(file, files, root = temp_dir)
}
