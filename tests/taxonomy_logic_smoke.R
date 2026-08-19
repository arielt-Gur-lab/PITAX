# Smoke tests for the v2.14 evidence-first taxonomy engine.
# Run with: Rscript tests/taxonomy_logic_smoke.R

get_this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }
  normalizePath(getwd(), mustWork = TRUE)
}

test_dir <- get_this_script_dir()
app_dir <- normalizePath(file.path(test_dir, ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "taxonomy_tools.R"))

make_hits <- function(species, genus, bits, identity, coverage, accession_prefix="ACC") {
  n <- length(bits)
  data.frame(
    rank=seq_len(n), organism=species,
    record_title=paste(species, "sequence"),
    accession=paste0(accession_prefix, seq_len(n), ".1"),
    taxid=seq_len(n)+1000L,
    identity_percent=identity,
    query_coverage_percent=coverage,
    evalue=rep(0, n), bit_score=bits,
    alignment_length=rep(500L,n), hsp_count=rep(1L,n), match_support=rep("Strong",n),
    kingdom=rep("Fungi",n), phylum=rep("Ascomycota",n), class=rep("Sordariomycetes",n),
    order=rep("Hypocreales",n), family=rep("Exampleaceae",n), genus=genus, species=species,
    reference_quality=rep("Standard GenBank",n), stringsAsFactors=FALSE
  )
}

# 1. Strong best species match but ITS sequence evidence is only Moderate:
# retain the genus and keep the best molecular species visible.
h1 <- make_hits(
  c(rep("Granulobasidium vellereum",15), "Cunninghammyces umbonatus"),
  c(rep("Granulobasidium",15), "Cunninghammyces"),
  c(1114,1110,1096,1089,1085,1082,1079,1076,1060,1015,1014,1014,1011,980,969,750),
  c(rep(98.2,15),86.9), rep(98,16), "GRA")
r1 <- build_taxonomic_consensus(h1, target="ITS", top_n=16)$summary
stopifnot(r1$recommended_identification == "Granulobasidium")
stopifnot(r1$recommended_level == "genus")
stopifnot(r1$best_match_species == "Granulobasidium vellereum")

# 2. High-quality full-length match with no close named species alternative.
h2 <- make_hits(rep("Pleurotus pulmonarius",10), rep("Pleurotus",10),
                seq(900,873,length.out=10), rep(99.8,10), rep(100,10), "PLE")
r2 <- build_taxonomic_consensus(h2, target="ITS", top_n=10)$summary
stopifnot(r2$recommended_identification == "Pleurotus pulmonarius")
stopifnot(r2$recommended_level == "species")
stopifnot(r2$close_species_count == 0)

# 3. Nearly indistinguishable species from the same genus: species unresolved,
# genus retained.
h3 <- make_hits(
  c("Fusarium alpha","Fusarium beta","Fusarium alpha","Fusarium beta","Fusarium alpha"),
  rep("Fusarium",5), c(900,899,898,897,896),
  c(99.8,99.7,99.6,99.6,99.5), rep(100,5), "FUS")
r3 <- build_taxonomic_consensus(h3, target="ITS", top_n=5)$summary
stopifnot(r3$recommended_identification == "Fusarium")
stopifnot(r3$recommended_level == "genus")
stopifnot(r3$species_discrimination == "Poor")

# 4. High agreement does not manufacture sequence quality: 93% ITS can still
# support genus with Low/review sequence evidence.
h4 <- make_hits(rep("Pleurotus pulmonarius",8), rep("Pleurotus",8),
                seq(572,565,length.out=8), rep(93.0,8), rep(100,8), "P93")
r4 <- build_taxonomic_consensus(h4, target="ITS", top_n=8)$summary
stopifnot(r4$recommended_identification == "Pleurotus")
stopifnot(r4$recommended_level == "genus")
stopifnot(r4$sequence_evidence == "Low")

# 5. A close match in a different genus prevents a genus-level call.
h5 <- make_hits(
  c(rep("Pleurotus pulmonarius",5), "Lentinus sajor-caju"),
  c(rep("Pleurotus",5), "Lentinus"),
  c(1107,1105,1104,1103,1102,1101),
  c(rep(99.84,5),99.60), rep(100,6), "PLX")
r5 <- build_taxonomic_consensus(h5, target="ITS", top_n=6)$summary
stopifnot(r5$recommended_level == "unresolved" || grepl("^LCA:", r5$recommended_level))
stopifnot(r5$genus_discrimination == "Poor")

# 6. Database abundance is context, not a vote. The best molecular match is
# F. ipomoeae (100/99.2), although F. chlamydosporum has many more accessions.
sp <- c(rep("Fusarium chlamydosporum",23), rep("Fusarium incarnatum",4),
        rep("Fusarium oxysporum",4), rep("Fusarium equiseti",2), rep("Fusarium ipomoeae",2))
id <- c(rep(100.0,23), rep(99.7,4), rep(99.5,4), rep(99.4,2), rep(100.0,2))
cov <- c(rep(99.0,23), rep(99.1,4), rep(99.0,4), rep(98.9,2), rep(99.2,2))
bits <- seq(1107, 1107-length(sp)+1)
h6 <- make_hits(sp, rep("Fusarium",length(sp)), bits, id, cov, "FCX")
r6all <- build_taxonomic_consensus(h6, target="ITS", top_n=nrow(h6))
r6 <- r6all$summary
p6 <- r6all$counts
stopifnot(r6$best_match_species == "Fusarium ipomoeae")
stopifnot(r6$closest_alternative_species == "Fusarium chlamydosporum")
stopifnot(r6$recommended_identification == "Fusarium")
stopifnot(r6$recommended_level == "genus")
stopifnot(r6$species_discrimination == "Poor")
stopifnot(r6$locus_discrimination == "Poor at species level")
stopifnot(p6$accession_count[p6$taxon == "Fusarium chlamydosporum"] == 23)
stopifnot(p6$accession_count[p6$taxon == "Fusarium ipomoeae"] == 2)

# 7. A short 100%-identity fragment must not outrank a near-full-length 99.8%
# hit when a near-full coverage tier exists.
h7 <- make_hits(c("Fusarium shortmatch","Fusarium fullmatch"), rep("Fusarium",2),
                c(600,900), c(100,99.8), c(55,100), "COV")
r7 <- build_taxonomic_consensus(h7, target="ITS", top_n=2)$summary
stopifnot(r7$best_match_species == "Fusarium fullmatch")


# 8. Within the same species, a marginal Identity advantage must not select a
# noticeably shorter representative when a near-full hit is available.
h8 <- make_hits(rep("Pleurotus pulmonarius",2), rep("Pleurotus",2),
                c(1138,1193), c(99.84,99.80), c(94.9,100.0), "PCM")
r8 <- build_taxonomic_consensus(h8, target="ITS", top_n=2)$summary
stopifnot(r8$best_match_accession == "PCM2.1")
stopifnot(r8$best_match_query_coverage_percent == 100)

cat("v2.14.2 taxonomy smoke tests passed.\n")
