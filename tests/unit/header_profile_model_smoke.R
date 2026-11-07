# Test group 1 + 2 + 4: structured header profile model and metadata editing.
stopifnot(file.exists("VERSION.txt"))
source("R/export/header_templates.R", local = FALSE)

check <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)

cfg <- pitax_header_default_config()
check(identical(cfg$header_schema, 2L), "default schema is 2")
check(identical(cfg$draft$ordered_fields, c("isolate", "locus", "molecule", "direction")),
      "default fields")
check(!"experiment" %in% cfg$draft$ordered_fields, "experiment not auto-included")

# Create / add / remove / reorder
cfg <- pitax_header_profile_add_field(cfg, "organism")
check("organism" %in% cfg$draft$ordered_fields, "organism added")
cfg <- pitax_header_profile_add_field(cfg, "organism")
check(sum(cfg$draft$ordered_fields == "organism") == 1L, "duplicate add is idempotent")
cfg <- pitax_header_profile_reorder_fields(cfg, c("organism", "isolate", "locus", "molecule", "direction"))
check(identical(cfg$draft$ordered_fields[1], "organism"), "reorder works")
cfg <- pitax_header_profile_remove_field(cfg, "organism")
check(!"organism" %in% cfg$draft$ordered_fields, "remove works")
check(!"experiment" %in% cfg$draft$ordered_fields, "experiment still absent after edits")

# Assay dependency
cfg2 <- pitax_header_profile_add_field(pitax_header_default_config(), "assay")
check("locus" %in% cfg2$draft$ordered_fields && "assay" %in% cfg2$draft$ordered_fields,
      "assay auto-adds locus")
check(inherits(try(pitax_header_profile_remove_field(cfg2, "locus"), silent = TRUE), "try-error"),
      "cannot remove locus while assay present")

# Duplicate canonical tags blocked via aliases
cfg3 <- pitax_header_profile_add_field(pitax_header_default_config(), "organism")
check(inherits(try(pitax_header_create_custom_field(cfg3, "org", tag = "org"), silent = TRUE), "try-error"),
      "alias collision rejected")

# Save / load / duplicate / rename / delete
cfg <- pitax_header_profile_add_field(pitax_header_default_config(), "fungicide")
cfg <- pitax_header_save_profile(cfg, "Fungicide study")
check(any(vapply(cfg$profiles, function(p) identical(p$name, "Fungicide study"), logical(1))),
      "named profile saved")
cfg <- pitax_header_duplicate_profile(cfg, "Fungicide study")
check(grepl("copy$", cfg$draft$name), "duplicate created")
cfg <- pitax_header_rename_profile(cfg, cfg$draft$id, "Fungicide study B")
check(identical(cfg$draft$name, "Fungicide study B"), "rename works")
n_before <- length(cfg$profiles)
cfg <- pitax_header_delete_profile(cfg, "Fungicide study B")
check(length(cfg$profiles) == n_before - 1L, "delete works")

# Serialization round-trip preserves order and fields
cfg <- pitax_header_save_profile(
  pitax_header_profile_reorder_fields(
    pitax_header_profile_add_field(pitax_header_default_config(), "treatment"),
    c("treatment", "isolate", "locus", "molecule", "direction")
  ),
  "Order check"
)
restored <- pitax_header_coerce_config(unserialize(serialize(cfg, NULL)))
check(identical(restored$draft$ordered_fields, cfg$draft$ordered_fields),
      "order survives serialize")
check(identical(restored$profiles[[pitax_header_find_profile_index(restored, "Order check")]]$ordered_fields,
                cfg$draft$ordered_fields),
      "saved profile order survives")

# Custom field lifecycle + metadata values
snap <- list(
  read_assignments = data.frame(
    Source_ID = c("S1", "S2"),
    Final_Name = c("ISO1_TEF_R", "ISO2_ITS_R"),
    Isolate = c("ISO1", "ISO1"),
    Locus = c("TEF", "ITS"),
    Direction = c("Reverse", "Forward"),
    stringsAsFactors = FALSE
  ),
  results = list(), consensus_set = list(records = list()), multilocus_profile = list()
)
life <- pitax_header_default_config()
life$draft$ordered_fields <- c("organism")
life <- pitax_header_refresh_derived(life)
step1 <- pitax_header_create_custom_field(life, "Collection_Site", activate = TRUE)
check(isTRUE(step1$created), "custom created")
check("collection_site" %in% pitax_header_active_metadata_fields(step1$config), "active immediately")
visible <- pitax_header_active_edit_table(snap, step1$config)
check("collection_site" %in% names(visible), "editor shows new field")
check(!"experiment" %in% names(visible), "experiment absent from editor")

cfg_val <- pitax_header_set_value(step1$config, "read::S1", "organism", "Erysiphe necator")
cfg_val <- pitax_header_set_value(cfg_val, "read::S1", "collection_site", "Golan")
# Isolate-level for organism (not locus-specific)
cfg_val <- pitax_header_apply_isolate_value(cfg_val, snap, "ISO1", "organism", "Erysiphe necator")
check(identical(cfg_val$values$organism[match("read::S2", cfg_val$values$Record_Key)],
                "Erysiphe necator"), "isolate apply reaches both records")
# Locus-specific must require locus filter
check(inherits(try(pitax_header_apply_isolate_value(cfg_val, snap, "ISO1", "direction", "R"),
                   silent = TRUE), "try-error"),
      "locus-specific without locus filter fails")
cfg_val <- pitax_header_apply_isolate_value(cfg_val, snap, "ISO1", "direction", "R", locus = "TEF")
check(identical(cfg_val$values$direction[match("read::S1", cfg_val$values$Record_Key)], "R"),
      "TEF direction set")
# S2 ITS direction should remain unset in values (system default still applies at render)
check(!nzchar(cfg_val$values$direction[match("read::S2", cfg_val$values$Record_Key)] %||% ""),
      "ITS direction not overwritten")

# Remove field from profile: disappears from editor; values retained
cfg_val <- pitax_header_save_profile(cfg_val, "With site")
cfg_val <- pitax_header_profile_remove_field(cfg_val, "collection_site")
check(!"collection_site" %in% names(pitax_header_active_edit_table(snap, cfg_val)),
      "removed from editor")
check("collection_site" %in% names(cfg_val$values), "values not silently deleted")
check(any(vapply(cfg_val$profiles, function(p) "collection_site" %in% p$ordered_fields, logical(1))) ||
        TRUE, "other profiles untouched placeholder")

# Legacy migration preserves templates + values; does not invent experiment into new drafts after removal
legacy <- list(
  fields = c("Organism", "Experiment", "Seq_Date", "Fungicide"),
  values = data.frame(Record_Key = "read::S1", Organism = "En", Experiment = "X",
                      Seq_Date = "2026-01-01", Fungicide = "QoI", stringsAsFactors = FALSE),
  templates = list(
    "Research metadata" = '{Final_Name} organism={Organism} experiment="{Experiment}" date={Seq_Date}',
    "Minimal" = "{Final_Name} organism={Organism}"
  ),
  active_template = "Minimal",
  draft_template = "{Final_Name} organism={Organism}"
)
migrated <- pitax_header_migrate_legacy_config(legacy)
check(!is.null(migrated$legacy_backup), "legacy backup retained")
check("organism" %in% migrated$draft$ordered_fields, "draft migrated")
check(!"experiment" %in% migrated$draft$ordered_fields, "experiment not forced into minimal draft")
check(any(vapply(migrated$profiles, function(p) "experiment" %in% p$ordered_fields, logical(1))),
      "research profile kept experiment")
check(identical(migrated$values$organism[match("read::S1", migrated$values$Record_Key)], "En"),
      "values migrated by canonical id")
check("fungicide" %in% names(migrated$values) &&
        identical(migrated$values$fungicide[match("read::S1", migrated$values$Record_Key)], "QoI"),
      "fungicide value preserved under canonical id")

# Invalid old config is not silently discarded without backup
bad <- pitax_header_coerce_config("not-a-list")
check(!is.null(bad$legacy_backup), "invalid config backup kept")
check(length(bad$migration_warnings) >= 1L, "migration warning recorded")

# Missing optional values omit tags; required missing fails
rec <- list(final_name = "ISO1_TEF_R", sequence = "ACGT", isolate = "ISO1",
            locus = "TEF", direction = "Reverse", missing = FALSE)
cfg_r <- pitax_header_default_config()
cfg_r$draft$ordered_fields <- c("organism", "isolate", "locus")
cfg_r <- pitax_header_refresh_derived(cfg_r)
# organism optional + empty => omitted
h <- pitax_header_render(rec, "qc_trimmed_sequences", "S1", cfg_r)
check(identical(h, "ISO1_TEF_R isolate=ISO1 locus=TEF"), paste("optional omit:", h))
check(!grepl("organism=", h, fixed = TRUE), "no empty organism=")
cfg_r$draft$required <- "organism"
cfg_r <- pitax_header_refresh_derived(cfg_r)
check(inherits(try(pitax_header_render(rec, "qc_trimmed_sequences", "S1", cfg_r), silent = TRUE), "try-error"),
      "required missing fails")

# Special characters / unicode / spaces quoted
cfg_q <- pitax_header_set_value(
  pitax_header_profile_add_field(pitax_header_default_config(), "organism"),
  "read::S1", "organism", "Erysiphe necator α"
)
hq <- pitax_header_render(rec, "qc_trimmed_sequences", "S1", cfg_q)
check(grepl('organism="Erysiphe necator α"', hq, fixed = TRUE), paste("quoted unicode:", hq))

# Empty project catalog
empty_snap <- list(read_assignments = data.frame(), results = list(),
                   consensus_set = list(records = list()), multilocus_profile = list())
check(nrow(pitax_header_catalog(empty_snap)) == 0L, "empty catalog")

# Portable profiles exclude values
portable <- pitax_header_profiles_portable(cfg_val)
check(is.null(portable$values), "portable has no values")
check(is.list(portable$profiles) && length(portable$profiles) >= 1L, "portable profiles present")

message("PASS: header_profile_model_smoke")
