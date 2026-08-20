# ============================================================
# Export / checkpoint / BLAST helpers
# ============================================================

clean_fasta_name <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

wrap_sequence <- function(seq, width = 80) {
  if (is.null(seq) || !nzchar(seq)) return("")
  starts <- seq(1, nchar(seq), by = width)
  paste(substring(seq, starts, pmin(starts + width - 1, nchar(seq))), collapse = "\n")
}

make_fasta <- function(records, include_metadata = FALSE, summary_df = NULL) {
  out <- character()
  if (!length(records)) return("")
  for (original_name in names(records)) {
    record <- records[[original_name]]
    if (is.null(record) || !nzchar(record$seq)) next
    final_name <- clean_fasta_name(record$final_name)
    header <- paste0(">", final_name)
    if (include_metadata && !is.null(summary_df)) {
      sm <- summary_df[summary_df$sample_id == original_name, , drop = FALSE]
      if (nrow(sm) == 1) {
        header <- paste0(header,
          " target=", gsub("\\s+", "_", sm$target),
          " length=", sm$trimmed_length,
          " trim=", sm$trim_start, "-", sm$trim_end,
          " reason=", sm$reason)
      }
    }
    out <- c(out, header, wrap_sequence(record$seq, 80))
  }
  paste(out, collapse = "\n")
}

write_assignment_checkpoint_zip <- function(file, rename_df, read_assignments, architecture, settings) {
  temp_dir <- tempfile("pitax_assignment_")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  write.csv(rename_df, file.path(temp_dir, "rename_map.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(read_assignments, file.path(temp_dir, "read_assignments.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  architecture_dir <- file.path(temp_dir, "project_architecture")
  dir.create(architecture_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(architecture$isolates, file.path(architecture_dir, "isolates.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(architecture$loci, file.path(architecture_dir, "loci.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(architecture$reads, file.path(architecture_dir, "reads.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  app_version <- tryCatch(get("APP_VERSION", inherits = TRUE), error = function(e) "unknown")
  settings_lines <- c(
    paste0("application_version: ", app_version),
    paste0("exported_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "checkpoint_stage: renamed_and_assigned_before_trimming",
    unlist(lapply(names(settings), function(nm) paste0(nm, ": ", settings[[nm]])))
  )
  writeLines(settings_lines, file.path(temp_dir, "run_settings.txt"))
  files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  zip::zipr(file, files, root = temp_dir)
}

write_checkpoint_zip <- function(file, stage, records, summary_df, settings, rename_df = NULL, results = NULL,
                                 read_assignments = NULL, architecture = NULL) {
  temp_dir <- tempfile(paste0("sanger_", stage, "_"))
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  target_name <- clean_fasta_name(settings$target)
  write.csv(summary_df, file.path(temp_dir, paste0(target_name, "_", stage, "_summary.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(make_fasta(records, FALSE), file.path(temp_dir, paste0(target_name, "_", stage, ".fasta")))
  if (!is.null(rename_df)) {
    write.csv(rename_df, file.path(temp_dir, "rename_map.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  }
  if (is.data.frame(read_assignments) && nrow(read_assignments)) {
    write.csv(read_assignments, file.path(temp_dir, "read_assignments.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  }
  if (is.list(architecture)) {
    architecture_dir <- file.path(temp_dir, "project_architecture")
    dir.create(architecture_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.data.frame(architecture$isolates)) write.csv(architecture$isolates, file.path(architecture_dir, "isolates.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    if (is.data.frame(architecture$loci)) write.csv(architecture$loci, file.path(architecture_dir, "loci.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    if (is.data.frame(architecture$reads)) write.csv(architecture$reads, file.path(architecture_dir, "reads.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  }
  app_version <- tryCatch(get("APP_VERSION", inherits = TRUE), error = function(e) "unknown")
  settings_lines <- c(
    paste0("application_version: ", app_version),
    paste0("exported_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("checkpoint_stage: ", stage),
    unlist(lapply(names(settings), function(nm) paste0(nm, ": ", settings[[nm]])))
  )
  writeLines(settings_lines, file.path(temp_dir, "run_settings.txt"))

  # Keep the original QC evidence with every processing checkpoint. These are
  # deliberately rendered as portable PNG files even though the chromatogram
  # itself is interactive in the Shiny UI.
  if (!is.null(results) && length(results)) {
    qc_dir <- file.path(temp_dir, "QC_plots")
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
    for (sample_name in names(results)) {
      result <- results[[sample_name]]
      if (is.null(result) || is.null(result$metrics)) next
      png_file <- file.path(qc_dir, paste0(clean_fasta_name(sample_name), "_QC.png"))
      try(write_qc_plot_png(result, settings, png_file), silent = TRUE)
    }
  }

  # Export ambiguous-peak review flags as a separate auditable QC table.
  # The default export scope is the retained/trimmed sequence, matching the
  # default review scope in the application.
  if (!is.null(results) && length(results)) {
    qc_flags <- tryCatch(collect_ambiguous_peak_flags(results, scope = "trimmed", settings = settings), error = function(e) data.frame())
    if (!ncol(qc_flags)) {
      qc_flags <- data.frame(
        Sample=character(), Position=integer(), Call=character(), Competing_channel=character(),
        Called_signal=numeric(), Competitor_signal=numeric(), Peak_ratio=numeric(),
        Competitor_percent=numeric(), Severity=character(), Flag=character(),
        stringsAsFactors=FALSE
      )
    }
    write.csv(qc_flags, file.path(temp_dir, paste0(target_name, "_QC_ambiguous_peak_flags.csv")),
              row.names = FALSE, fileEncoding = "UTF-8")
  }

  # Manual curation is exported separately from automatic QC so every edited
  # base, manual trim, review decision, undo and redo remains auditable.
  if (!is.null(results) && length(results)) {
    curation_log <- tryCatch(collect_curation_log(results), error = function(e) data.frame())
    if (!ncol(curation_log)) {
      curation_log <- data.frame(
        Sample=character(), Timestamp=character(), Transaction_ID=character(), Revision=integer(),
        Action=character(), Position=integer(), Before=character(), After=character(), Method=character(),
        Evidence=character(), Details=character(), stringsAsFactors=FALSE
      )
    }
    write.csv(curation_log, file.path(temp_dir, paste0(target_name, "_manual_curation_log.csv")),
              row.names = FALSE, fileEncoding = "UTF-8")
  }

  files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  zip::zipr(file, files, root = temp_dir)
}

# Preliminary support classification from BLAST-style metrics.
# This is deliberately labeled "support" rather than taxonomic certainty.
blast_match_support <- function(identity, coverage, evalue) {
  identity <- suppressWarnings(as.numeric(identity))
  coverage <- suppressWarnings(as.numeric(coverage))
  evalue <- suppressWarnings(as.numeric(evalue))
  if (any(is.na(c(identity, coverage, evalue)))) return("Unscored")
  if (identity >= 99 && coverage >= 95 && evalue <= 1e-50) return("Strong")
  if (identity >= 97 && coverage >= 90 && evalue <= 1e-20) return("Moderate")
  "Weak / review"
}

extract_accession_from_hit <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(NA_character_)
  x <- as.character(x)
  # Common BLAST identifiers: ref|NR_...|, gb|AB123...|, emb|..., or a bare accession.version.
  m <- regmatches(x, regexpr("(?:ref|gb|emb|dbj|tpg|tpe|tpd)\\|[A-Z0-9_.-]+\\|", x, perl=TRUE))
  if (length(m) && nzchar(m)) return(sub("^[^|]+\\|([^|]+)\\|$", "\\1", m, perl=TRUE))
  m2 <- regmatches(x, regexpr("[A-Z]{1,4}_[A-Z0-9]+(?:\\.[0-9]+)?|[A-Z]{1,2}[0-9]{5,8}(?:\\.[0-9]+)?", x, perl=TRUE))
  if (length(m2) && nzchar(m2)) return(m2)
  NA_character_
}

xml_unescape_small <- function(x) {
  if (is.null(x) || is.na(x)) return("")
  x <- gsub("&amp;", "&", x, fixed=TRUE)
  x <- gsub("&lt;", "<", x, fixed=TRUE)
  x <- gsub("&gt;", ">", x, fixed=TRUE)
  x <- gsub("&quot;", '"', x, fixed=TRUE)
  x <- gsub("&apos;", "'", x, fixed=TRUE)
  x
}

extract_xml_tag <- function(xml, tag) {
  if (is.null(xml) || !nzchar(xml)) return("")
  pat <- paste0("<", tag, ">([\\s\\S]*?)</", tag, ">")
  m <- regexec(pat, xml, perl=TRUE)
  r <- regmatches(xml, m)[[1]]
  if (length(r) < 2) "" else xml_unescape_small(r[2])
}

# Retrieve authoritative title and organism from the NCBI Nucleotide record that
# corresponds to the BLAST accession. This avoids trying to infer a taxon solely
# from the accession string or a truncated BLAST subject identifier.
fetch_ncbi_nucleotide_metadata <- function(accession) {
  accession <- ifelse(is.null(accession) || is.na(accession), "", trimws(as.character(accession)))
  if (!nzchar(accession)) return(list(title="", organism=""))

  req <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi") |>
    httr2::req_url_query(db="nuccore", id=accession, rettype="gb", retmode="xml") |>
    httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI E-utilities client")
  xml <- httr2::resp_body_string(httr2::req_perform(req))
  list(
    title = extract_xml_tag(xml, "GBSeq_definition"),
    organism = extract_xml_tag(xml, "GBSeq_organism")
  )
}


# Retrieve title / organism metadata for many BLAST accessions in batches. The
# BLAST tabular fallback contains accession and alignment statistics but no
# subject title or organism columns, so every returned hit must be enriched --
# not only the first hit.
fetch_ncbi_nucleotide_metadata_batch <- function(accessions, batch_size = 20L) {
  accessions <- unique(trimws(as.character(accessions)))
  accessions <- accessions[!is.na(accessions) & nzchar(accessions)]
  if (!length(accessions)) {
    return(data.frame(
      accession = character(), title = character(), organism = character(),
      taxid = integer(), stringsAsFactors = FALSE
    ))
  }

  groups <- split(accessions, ceiling(seq_along(accessions) / max(1L, as.integer(batch_size))))
  rows <- list()

  for (g in seq_along(groups)) {
    ids <- groups[[g]]
    req <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi") |>
      httr2::req_url_query(
        db = "nuccore",
        id = paste(ids, collapse = ","),
        rettype = "gb",
        retmode = "xml"
      ) |>
      httr2::req_user_agent("Local Sanger Sequence Pipeline / NCBI E-utilities client")

    xml_txt <- tryCatch(
      httr2::resp_body_string(httr2::req_perform(req)),
      error = function(e) ""
    )

    if (nzchar(xml_txt)) {
      doc <- tryCatch(xml2::read_xml(xml_txt), error = function(e) NULL)
      if (!is.null(doc)) {
        seq_nodes <- xml2::xml_find_all(doc, ".//*[local-name()='GBSeq']")
        if (length(seq_nodes)) {
          for (k in seq_along(seq_nodes)) {
            node <- seq_nodes[[k]]
            text1 <- function(xpath) {
              z <- xml2::xml_find_first(node, xpath)
              if (inherits(z, "xml_missing")) "" else trimws(xml2::xml_text(z))
            }
            accession_version <- text1("./*[local-name()='GBSeq_accession-version']")
            primary_accession <- text1("./*[local-name()='GBSeq_primary-accession']")
            definition <- text1("./*[local-name()='GBSeq_definition']")
            organism <- text1("./*[local-name()='GBSeq_organism']")

            tax_values <- xml2::xml_find_all(
              node,
              ".//*[local-name()='GBQualifier'][./*[local-name()='GBQualifier_name']='db_xref']/*[local-name()='GBQualifier_value']"
            )
            tax_values <- if (length(tax_values)) xml2::xml_text(tax_values) else character()
            tax_match <- tax_values[grepl("^taxon:[0-9]+$", tax_values)]
            taxid <- if (length(tax_match)) suppressWarnings(as.integer(sub("^taxon:", "", tax_match[1]))) else NA_integer_

            acc <- if (nzchar(accession_version)) accession_version else primary_accession
            rows[[length(rows) + 1]] <- data.frame(
              accession = acc,
              primary_accession = primary_accession,
              title = definition,
              organism = organism,
              taxid = taxid,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    if (g < length(groups)) Sys.sleep(0.35)
  }

  if (!length(rows)) {
    return(data.frame(
      accession = character(), primary_accession = character(), title = character(),
      organism = character(), taxid = integer(), stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

organism_from_title_fallback <- function(title) {
  if (is.null(title) || is.na(title) || !nzchar(title)) return("")
  title <- trimws(as.character(title))
  # Prefer a bracketed taxon if present.
  b <- regmatches(title, regexpr("\\[[^][]+\\]$", title, perl=TRUE))
  if (length(b) && nzchar(b)) return(sub("^\\[|\\]$", "", b))
  # Conservative fallback: first binomial-looking pair only. This is a display
  # aid, not a taxonomic assertion; EFetch metadata is preferred whenever possible.
  m <- regmatches(title, regexpr("\\b[A-Z][a-z][A-Za-z.-]+\\s+[a-z][A-Za-z.-]+\\b", title, perl=TRUE))
  if (length(m) && nzchar(m)) m else ""
}

parse_blast_csv_top_hit <- function(text) {
  if (is.null(text) || !nzchar(text)) return(NULL)
  if (grepl("Status=WAITING|Status=UNKNOWN|Status=FAILED", text, fixed = FALSE)) return(NULL)
  df <- tryCatch(read.csv(text = text, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)

  nms <- tolower(names(df))
  find_col <- function(patterns) {
    idx <- which(vapply(nms, function(x) any(grepl(paste(patterns, collapse = "|"), x)), logical(1)))
    if (length(idx)) idx[1] else NA_integer_
  }

  subj_col <- find_col(c("subject", "description", "title"))
  acc_col <- find_col(c("accession", "acc"))
  ident_col <- find_col(c("identity", "pident"))
  cover_col <- find_col(c("query cover", "qcov"))
  eval_col <- find_col(c("e.?value", "expect"))

  row <- df[1, , drop = FALSE]
  getv <- function(col) if (is.na(col)) NA else row[[col]][1]
  top_hit <- as.character(getv(subj_col))
  accession <- as.character(getv(acc_col))
  if (is.na(accession) || !nzchar(accession)) accession <- extract_accession_from_hit(top_hit)
  ident <- suppressWarnings(as.numeric(gsub("%", "", getv(ident_col))))
  cover <- suppressWarnings(as.numeric(gsub("%", "", getv(cover_col))))
  ev <- suppressWarnings(as.numeric(getv(eval_col)))

  data.frame(
    top_hit = top_hit,
    organism = organism_from_title_fallback(top_hit),
    record_title = top_hit,
    accession = accession,
    identity_percent = ident,
    query_coverage_percent = cover,
    evalue = ev,
    match_support = blast_match_support(ident, cover, ev),
    stringsAsFactors = FALSE
  )
}



# ============================================================
# NCBI BLAST XML2 parsing
# ============================================================

xml_local_text <- function(node, tag) {
  if (inherits(node, "xml_missing") || length(node) == 0) return("")
  # NCBI XML2 serializations may use either concise element names (e.g.
  # accession) or schema-style names (e.g. HitDescr_accession). Accept both.
  xpath <- paste0(
    ".//*[local-name()='", tag, "' or substring-after(local-name(), '_')='", tag, "']"
  )
  found <- xml2::xml_find_first(node, xpath)
  if (inherits(found, "xml_missing")) "" else trimws(xml2::xml_text(found))
}

covered_query_bases <- function(starts, ends) {
  starts <- suppressWarnings(as.integer(starts))
  ends <- suppressWarnings(as.integer(ends))
  good <- is.finite(starts) & is.finite(ends)
  starts <- starts[good]; ends <- ends[good]
  if (!length(starts)) return(0L)
  lo <- pmin(starts, ends); hi <- pmax(starts, ends)
  ord <- order(lo, hi)
  lo <- lo[ord]; hi <- hi[ord]
  total <- 0L
  cur_lo <- lo[1]; cur_hi <- hi[1]
  if (length(lo) > 1) {
    for (i in 2:length(lo)) {
      if (lo[i] <= cur_hi + 1L) {
        cur_hi <- max(cur_hi, hi[i])
      } else {
        total <- total + (cur_hi - cur_lo + 1L)
        cur_lo <- lo[i]; cur_hi <- hi[i]
      }
    }
  }
  total + (cur_hi - cur_lo + 1L)
}

# Parse the supported NCBI BLAST XML2 result format. Unlike the tabular CSV
# response, XML2 carries subject titles, accessions and (for current NCBI
# reports) scientific-name/taxonomy fields, so the application does not have to
# guess column positions or infer a species from an accession number.
parse_blast_xml2_hits <- function(text, query_length = NA_integer_, max_hits = Inf) {
  if (is.null(text) || !nzchar(text)) return(data.frame())
  if (grepl("Status=WAITING|Status=UNKNOWN|Status=FAILED", text, perl = TRUE)) return(data.frame())

  doc <- tryCatch(xml2::read_xml(text), error = function(e) NULL)
  if (is.null(doc)) return(data.frame())

  if (!is.finite(query_length) || query_length <= 0) {
    qlen_text <- xml_local_text(doc, "query-len")
    query_length <- suppressWarnings(as.integer(qlen_text))
  }

  hits <- xml2::xml_find_all(doc, ".//*[local-name()='Hit' or substring-after(local-name(), '_')='Hit']")
  if (!length(hits)) return(data.frame())
  if (is.finite(max_hits)) hits <- hits[seq_len(min(length(hits), as.integer(max_hits)))]

  rows <- vector("list", length(hits))
  for (i in seq_along(hits)) {
    hit <- hits[[i]]
    descr <- xml2::xml_find_first(hit, ".//*[local-name()='HitDescr' or substring-after(local-name(), '_')='HitDescr']")

    accession <- xml_local_text(descr, "accession")
    hit_id <- xml_local_text(descr, "id")
    title <- xml_local_text(descr, "title")
    organism <- xml_local_text(descr, "sciname")
    taxid <- suppressWarnings(as.integer(xml_local_text(descr, "taxid")))
    if (!nzchar(accession)) accession <- extract_accession_from_hit(hit_id)
    if (!nzchar(title)) title <- hit_id
    if (!nzchar(organism)) organism <- organism_from_title_fallback(title)

    hsps <- xml2::xml_find_all(hit, ".//*[local-name()='Hsp' or substring-after(local-name(), '_')='Hsp']")
    if (!length(hsps)) next

    hsp_df <- do.call(rbind, lapply(seq_along(hsps), function(j) {
      h <- hsps[[j]]
      data.frame(
        identity = suppressWarnings(as.numeric(xml_local_text(h, "identity"))),
        align_len = suppressWarnings(as.numeric(xml_local_text(h, "align-len"))),
        q_from = suppressWarnings(as.integer(xml_local_text(h, "query-from"))),
        q_to = suppressWarnings(as.integer(xml_local_text(h, "query-to"))),
        evalue = suppressWarnings(as.numeric(xml_local_text(h, "evalue"))),
        bit_score = suppressWarnings(as.numeric(xml_local_text(h, "bit-score"))),
        stringsAsFactors = FALSE
      )
    }))

    hsp_df <- hsp_df[is.finite(hsp_df$align_len) & hsp_df$align_len > 0, , drop = FALSE]
    if (!nrow(hsp_df)) next
    # Subject-level score and representative identity are intentionally
    # separated. BLAST ranking is driven by the strongest HSP (best E-value /
    # Bit score), while Identity should describe the broadest individual HSP
    # rather than a tiny perfect local segment.
    hsp_df$q_span <- abs(hsp_df$q_to - hsp_df$q_from) + 1
    score_order <- order(hsp_df$evalue, -hsp_df$bit_score, na.last = TRUE)
    score_best <- hsp_df[score_order[1], , drop = FALSE]
    representative_order <- order(-hsp_df$q_span, hsp_df$evalue, -hsp_df$bit_score, na.last = TRUE)
    representative <- hsp_df[representative_order[1], , drop = FALSE]

    identity_pct <- 100 * representative$identity / representative$align_len
    covered <- covered_query_bases(hsp_df$q_from, hsp_df$q_to)
    coverage_pct <- if (is.finite(query_length) && query_length > 0) 100 * covered / query_length else NA_real_

    rows[[i]] <- data.frame(
      rank = i,
      organism = organism,
      record_title = title,
      accession = accession,
      taxid = taxid,
      identity_percent = round(identity_pct, 3),
      query_coverage_percent = if (is.finite(coverage_pct)) round(min(100, coverage_pct), 2) else NA_real_,
      evalue = score_best$evalue,
      bit_score = score_best$bit_score,
      alignment_length = as.integer(representative$align_len),
      hsp_count = nrow(hsp_df),
      match_support = blast_match_support(identity_pct, coverage_pct, score_best$evalue),
      stringsAsFactors = FALSE
    )
  }

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Fallback for the exact 12-column CSV returned by NCBI's tabular CSV view.
# NCBI does not include a header in this response; the fields correspond to
# qseqid,sacc,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore.
#
# IMPORTANT (v2.9): the tabular response can contain multiple HSP rows for the
# same accession. Those are local alignment segments belonging to ONE BLAST hit,
# not independent database hits. We therefore aggregate by accession before the
# rows enter taxonomy/confidence calculations.
parse_blast_csv_hits_fallback <- function(text, query_length = NA_integer_) {
  if (is.null(text) || !nzchar(text)) return(data.frame())
  lines <- strsplit(text, "\\r?\\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !grepl("^#|^Status=", lines)]
  if (!length(lines)) return(data.frame())

  df <- tryCatch(utils::read.csv(
    text = paste(lines, collapse = "\n"), header = FALSE,
    stringsAsFactors = FALSE, check.names = FALSE
  ), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 12 || !nrow(df)) return(data.frame())
  df <- df[, 1:12, drop = FALSE]
  names(df) <- c("qseqid","accession","identity_percent","alignment_length","mismatch","gapopen",
                 "qstart","qend","sstart","send","evalue","bit_score")

  df$accession <- as.character(df$accession)
  df$identity_percent <- suppressWarnings(as.numeric(df$identity_percent))
  df$alignment_length <- suppressWarnings(as.integer(df$alignment_length))
  df$qstart <- suppressWarnings(as.integer(df$qstart))
  df$qend <- suppressWarnings(as.integer(df$qend))
  df$evalue <- suppressWarnings(as.numeric(df$evalue))
  df$bit_score <- suppressWarnings(as.numeric(df$bit_score))

  accessions <- unique(df$accession[nzchar(df$accession)])
  if (!length(accessions)) return(data.frame())

  rows <- lapply(accessions, function(acc) {
    g <- df[df$accession == acc, , drop = FALSE]

    # Subject-level score comes from the strongest HSP. Representative
    # Identity comes from the HSP spanning the largest portion of the query.
    score_ord <- order(g$evalue, -g$bit_score, -g$alignment_length, na.last = TRUE)
    score_best <- g[score_ord[1], , drop = FALSE]
    q_span <- abs(g$qend - g$qstart) + 1
    rep_ord <- order(-q_span, g$evalue, -g$bit_score, na.last = TRUE)
    representative <- g[rep_ord[1], , drop = FALSE]

    covered <- covered_query_bases(g$qstart, g$qend)
    coverage_pct <- if (is.finite(query_length) && query_length > 0) {
      min(100, 100 * covered / query_length)
    } else {
      NA_real_
    }

    data.frame(
      rank = NA_integer_,
      organism = "",
      record_title = "",
      accession = acc,
      taxid = NA_integer_,
      identity_percent = round(representative$identity_percent[1], 3),
      query_coverage_percent = if (is.finite(coverage_pct)) round(coverage_pct, 2) else NA_real_,
      evalue = score_best$evalue[1],
      bit_score = score_best$bit_score[1],
      alignment_length = representative$alignment_length[1],
      hsp_count = nrow(g),
      match_support = blast_match_support(representative$identity_percent[1], coverage_pct, score_best$evalue[1]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  ord <- order(out$evalue, -out$bit_score, -out$identity_percent, -out$query_coverage_percent, na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

# Normalize a stored BLAST hit table so each RID/accession contributes at most
# one row to downstream taxonomy. This also repairs v2.8 project data in which
# CSV HSP rows may already have been stored as if they were separate hits.
normalize_blast_hits_unique_accession <- function(hits) {
  if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) return(data.frame())
  if (!"accession" %in% names(hits)) return(hits)

  # Preserve rows without an accession separately; they cannot be safely merged.
  acc <- as.character(hits$accession)
  with_acc <- hits[nzchar(acc) & !is.na(acc), , drop = FALSE]
  without_acc <- hits[!nzchar(acc) | is.na(acc), , drop = FALSE]
  if (!nrow(with_acc)) return(hits)

  key_cols <- intersect(c("original_name", "rid", "accession"), names(with_acc))
  if (!length(key_cols)) key_cols <- "accession"
  keys <- do.call(paste, c(lapply(with_acc[key_cols], as.character), sep = "\r"))
  groups <- split(seq_len(nrow(with_acc)), keys)

  first_nonempty <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x)) x[1] else ""
  }

  rows <- lapply(groups, function(idx) {
    g <- with_acc[idx, , drop = FALSE]

    cov <- if ("query_coverage_percent" %in% names(g)) suppressWarnings(as.numeric(g$query_coverage_percent)) else rep(NA_real_, nrow(g))
    ev <- if ("evalue" %in% names(g)) suppressWarnings(as.numeric(g$evalue)) else rep(NA_real_, nrow(g))
    bs <- if ("bit_score" %in% names(g)) suppressWarnings(as.numeric(g$bit_score)) else rep(NA_real_, nrow(g))
    al <- if ("alignment_length" %in% names(g)) suppressWarnings(as.numeric(g$alignment_length)) else rep(NA_real_, nrow(g))

    # For legacy HSP rows, prefer the row with the widest query coverage first;
    # then strongest E-value/bit score. This recovers the biologically relevant
    # full-length HSP in the common case seen in v2.8.
    ord <- order(-cov, ev, -bs, -al, na.last = TRUE)
    rep_row <- g[ord[1], , drop = FALSE]

    # Preserve representative Identity/Coverage from the widest row, but use
    # the strongest subject-level score across all stored HSP rows.
    if ("bit_score" %in% names(rep_row) && any(is.finite(bs))) rep_row$bit_score <- max(bs[is.finite(bs)], na.rm = TRUE)
    if ("evalue" %in% names(rep_row) && any(is.finite(ev))) rep_row$evalue <- min(ev[is.finite(ev)], na.rm = TRUE)

    if ("organism" %in% names(rep_row)) rep_row$organism <- first_nonempty(g$organism)
    if ("record_title" %in% names(rep_row)) rep_row$record_title <- first_nonempty(g$record_title)
    if ("taxid" %in% names(rep_row)) {
      tx <- suppressWarnings(as.integer(g$taxid))
      tx <- tx[is.finite(tx)]
      rep_row$taxid <- if (length(tx)) tx[1] else NA_integer_
    }
    if ("hsp_count" %in% names(rep_row)) {
      hc <- suppressWarnings(as.integer(g$hsp_count))
      hc[!is.finite(hc)] <- 1L
      rep_row$hsp_count <- sum(hc)
    }
    if ("match_support" %in% names(rep_row) && all(c("identity_percent","query_coverage_percent","evalue") %in% names(rep_row))) {
      rep_row$match_support <- blast_match_support(
        suppressWarnings(as.numeric(rep_row$identity_percent[1])),
        suppressWarnings(as.numeric(rep_row$query_coverage_percent[1])),
        suppressWarnings(as.numeric(rep_row$evalue[1]))
      )
    }
    rep_row
  })

  out <- do.call(rbind, rows)
  if (nrow(without_acc)) out <- rbind(out, without_acc)

  # Re-rank inside each RID (or globally if RID is unavailable).
  rerank_group <- function(g) {
    ev <- if ("evalue" %in% names(g)) suppressWarnings(as.numeric(g$evalue)) else rep(NA_real_, nrow(g))
    bs <- if ("bit_score" %in% names(g)) suppressWarnings(as.numeric(g$bit_score)) else rep(NA_real_, nrow(g))
    id <- if ("identity_percent" %in% names(g)) suppressWarnings(as.numeric(g$identity_percent)) else rep(NA_real_, nrow(g))
    cv <- if ("query_coverage_percent" %in% names(g)) suppressWarnings(as.numeric(g$query_coverage_percent)) else rep(NA_real_, nrow(g))
    ord <- order(ev, -bs, -id, -cv, na.last = TRUE)
    g <- g[ord, , drop = FALSE]
    if ("rank" %in% names(g)) g$rank <- seq_len(nrow(g))
    g
  }

  if ("rid" %in% names(out)) {
    parts <- split(seq_len(nrow(out)), as.character(out$rid))
    out <- do.call(rbind, lapply(parts, function(idx) rerank_group(out[idx, , drop = FALSE])))
  } else {
    out <- rerank_group(out)
  }
  rownames(out) <- NULL
  out
}
