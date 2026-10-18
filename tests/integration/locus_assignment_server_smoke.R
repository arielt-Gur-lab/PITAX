# Exercise the real server observers, project round-trip and result invalidation.
app_env <- new.env(parent = globalenv())
sys.source("app.R", envir = app_env)
check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

publish_editor <- function(session, input, assay_id, revision) {
  values <- shiny::reactiveValuesToList(input)
  values$assay_editor_snapshot <- NULL
  session$setInputs(assay_editor_snapshot = list(
    assay_id = assay_id, revision = revision, values = values
  ))
}

saved_custom_project <- NULL
shiny::testServer(app_env$server, {
  session$flushReact()
  session$setInputs(
    target = "ITS", assay_name = "ITS", sequencing_primer = "Forward",
    expected_amplicon_len = 650L, absolute_max_base_index = 680L
  )
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  fixture <- normalizePath("tests/fixtures/Stage3_synthetic_pair/TESTPAIR001_forward.ab1", mustWork = TRUE)
  uploads <- data.frame(name = c("BAR.1.ab1", "BAR.2.ab1"), size = rep(file.info(fixture)$size, 2),
                        type = "application/octet-stream", datapath = rep(fixture, 2), stringsAsFactors = FALSE)
  session$setInputs(ab1_files = uploads)
  session$flushReact()
  check(nrow(rv$read_assignments) == 2L, "Upload did not create assignment rows.")

  session$setInputs(locus_add = 1)
  session$setInputs(locus_new_id = "sdhc", locus_display_name = "Subunit C", locus_save = 1)
  check("SDHC" %in% rv$locus_vocabulary$Locus_ID, "Add-locus observer did not register the entry.")
  check(rv$assay_profiles$Locus_ID[1] == "SDHC", "New locus was not applied to the active assay.")
  # Mock sessions do not apply updateSelectInput messages back to their inputs.
  session$setInputs(target = "SDHC")
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  check(all(rv$read_assignments$Locus == "SDHC"), "Assigned reads did not inherit the new locus.")

  # First demonstrate why searching an immutable barcode in blank isolate fields cannot match.
  session$setInputs(batch_text_source = "isolate", batch_isolate_find = "BAR.", batch_isolate_replace = "FB",
                    batch_isolate_prefix = "", batch_isolate_suffix = "", batch_assay = "", batch_direction = "",
                    apply_assignment_batch = 1)
  check(all(rv$read_assignments$Isolate == ""), "Blank isolates were inferred from source IDs.")
  session$setInputs(batch_text_source = "source", apply_assignment_batch = 2)
  check(identical(rv$read_assignments$Isolate, c("FB1", "FB2")), "Barcode-to-isolate observer failed.")
  check(identical(rv$read_assignments$Source_ID, c("BAR.1", "BAR.2")), "Source IDs were changed.")
  check(identical(rv$read_assignments$Final_Name, c("FB1_SDHC_F", "FB2_SDHC_F")), "Generated names did not update.")

  session$setInputs(window = 25L, min_peak_ratio = 3, min_relative_signal = 0.2,
                    min_len_before_collapse = 350L, bad_run_windows = 12L, min_usable_len = 100L,
                    consensus_min_overlap = 100L, consensus_min_identity = 90,
                    consensus_quality_delta = 10L, consensus_strong_quality = 30L,
                    run_trimming = 1)
  check(length(rv$results) == 2L, "Actual AB1 processing failed for a custom locus.")
  check(all(vapply(rv$results, function(r) identical(r$processing_settings$target, "SDHC"), logical(1))),
        "Actual processing remapped the custom locus.")
  check(isTRUE(build_analysis_sequences(notify = FALSE)), "Analysis sequence build failed for a custom locus.")
  check(length(rv$consensus_set$records) == 2L, "Custom locus analysis records were lost.")

  # Collect still-pending table edits and only replace the checked row.
  args <- list(batch_text_source = "isolate", batch_isolate_find = "DRAFT", batch_isolate_replace = "DONE",
               apply_assignment_batch = 3)
  args[[assignment_input_id("assign_isolate", rv$read_assignments$Read_ID[1])]] <- "DRAFT1"
  args[[assignment_input_id("assign_isolate", rv$read_assignments$Read_ID[2])]] <- "DRAFT2"
  args[[assignment_input_id("assign_select", rv$read_assignments$Read_ID[1])]] <- TRUE
  do.call(session$setInputs, args)
  check(identical(rv$read_assignments$Isolate, c("DONE1", "DRAFT2")), "Pending edits or checkbox scope were lost.")
  check(length(rv$consensus_set$records) == 0L, "Renaming left old analysis sequences active.")

  # Renaming a display label preserves identity and survives save/load, including unused entries.
  session$setInputs(locus_edit = 1)
  session$setInputs(locus_display_name = "Succinate C", locus_save = 2)
  check(rv$assay_profiles$Locus_Display_Name[1] == "Succinate C", "Edited display label did not propagate.")
  rv$locus_vocabulary <- app_env$pitax_save_locus(rv$locus_vocabulary, "SDHD", "Subunit D")
  saved <- make_project_bundle()
  save_path <- tempfile(fileext = ".sangerproject")
  saveRDS(saved, save_path)
  saved <- readRDS(save_path); unlink(save_path)
  saved_custom_project <<- saved
  reset_pipeline_state()
  check(!"SDHC" %in% rv$locus_vocabulary$Locus_ID, "New-run reset leaked the prior catalog.")
  loaded <- apply_loaded_project_object(saved, notify = FALSE)
  check(isTRUE(loaded$ok), paste("Custom-locus project failed to load:", loaded$message))
  check(all(c("SDHC", "SDHD") %in% rv$locus_vocabulary$Locus_ID), "Project load lost custom or unused loci.")
  check(rv$assay_profiles$Locus_ID[1] == "SDHC" && rv$assay_profiles$Locus_Display_Name[1] == "Succinate C",
        "Project load remapped locus identity or display name.")
  snap <- app_env$pitax_export_snapshot_from_rv(rv, input)
  check("SDHD" %in% snap$locus_vocabulary$Locus_ID, "Export snapshot lost the project catalog.")

  invalid <- saved
  invalid$state$locus_vocabulary <- NULL
  rejected <- apply_loaded_project_object(invalid, notify = FALSE)
  check(!isTRUE(rejected$ok), "An unregistered custom profile was silently accepted on load.")
  check("SDHC" %in% rv$locus_vocabulary$Locus_ID, "Failed load mutated the current catalog.")
  legacy <- saved
  legacy$state$locus_vocabulary <- NULL
  legacy$state$assay_profiles <- app_env$assay_default_profiles()
  legacy$state$read_assignments <- app_env$stage2_empty_assignments()
  legacy$state$results <- list()
  legacy$state$summary <- NULL
  legacy$state$consensus_set <- app_env$stage3_empty_consensus_set()
  legacy$state$settings <- app_env$assay_resolve_read_settings(legacy$state$assay_profiles)
  check(isTRUE(apply_loaded_project_object(legacy, notify = FALSE)$ok), "Old schema-6 projects no longer load.")
})
# A loaded project has no browser file-upload input. Dropdown edits must still
# synchronize its saved read assignments.
shiny::testServer(app_env$server, {
  session$flushReact()
  check(isTRUE(apply_loaded_project_object(saved_custom_project, notify = FALSE)$ok), "Fresh-session restore failed.")
  session$flushReact()
  check(is.null(input$ab1_files), "This regression requires a project with no uploaded files.")
  session$setInputs(target = "SDHD", assay_name = "Custom assay", sequencing_primer = "Forward",
                    expected_amplicon_len = 650L, absolute_max_base_index = 680L)
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  check(all(rv$read_assignments$Locus == "SDHD"), "Locus dropdown changes did not update a loaded project.")
  check(all(grepl("_SDHD_F$", rv$read_assignments$Final_Name)), "Loaded-project final names did not follow the locus.")
})
# Naming mode belongs to the assay, including when the active editor changes.
saved_naming_project <- NULL
shiny::testServer(app_env$server, {
  session$flushReact()
  check(isTRUE(apply_loaded_project_object(saved_custom_project, notify = FALSE)$ok), "Naming fixture restore failed.")
  session$flushReact()
  session$setInputs(target = "SDHC", assay_name = "Alpha Assay", assay_name_source = "assay",
                    sequencing_primer = "Forward", expected_amplicon_len = 650L, absolute_max_base_index = 680L)
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  check(all(grepl("_Alpha-Assay_F$", rv$read_assignments$Final_Name)), "Naming toggle did not update existing reads.")
  check(all(rv$read_assignments$Locus == "SDHC"), "Naming toggle changed the biological locus.")
  first_id <- selected_assay_id()
  session$setInputs(assay_add_profile = 1)
  second_id <- selected_assay_id()
  session$setInputs(assay_editor_select = second_id)
  session$setInputs(target = "LSU", assay_name = "Beta Assay", assay_name_source = "locus")
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  rv$read_assignments$Assay_ID[2] <- second_id
  sync_assignment_state()
  check(rv$assay_profiles$Name_Source[match(first_id, rv$assay_profiles$Assay_ID)] == "assay", "Naming choice leaked from the second profile into the first.")
  check(identical(rv$read_assignments$Name_Label, c("Alpha Assay", "LSU")), "Per-profile naming choices were not applied independently.")
  session$setInputs(consensus_min_overlap = 100L, consensus_min_identity = 90,
                    consensus_quality_delta = 10L, consensus_strong_quality = 30L)
  check(isTRUE(build_analysis_sequences(notify = FALSE)), "Analysis setup failed for naming invalidation test.")
  session$setInputs(assay_name_source = "assay")
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  check(rv$read_assignments$Name_Label[2] == "Beta Assay", "Second profile toggle was ignored.")
  check(length(rv$consensus_set$records) == 0L, "Naming change left outdated analysis names active.")
  unchanged <- rv$assay_profiles
  session$setInputs(assay_editor_snapshot = list(assay_id = first_id, revision = assay_editor_revision(),
                    values = list(assay_name = "Stale wrong profile", target = "ITS", assay_name_source = "locus",
                                  expected_amplicon_len = 650L, absolute_max_base_index = 680L)))
  check(identical(rv$assay_profiles, unchanged), "An editor event from the previous profile overwrote the active profile.")
  session$setInputs(assay_editor_snapshot = list(assay_id = second_id, revision = assay_editor_revision() - 1L,
                    values = list(assay_name = "Stale revision", target = "ITS", assay_name_source = "locus",
                                  expected_amplicon_len = 650L, absolute_max_base_index = 680L)))
  check(identical(rv$assay_profiles, unchanged), "An old editor revision overwrote the current profile.")
  saved_naming_project <<- unserialize(serialize(make_project_bundle(), NULL))
})
shiny::testServer(app_env$server, {
  session$flushReact()
  loaded <- apply_loaded_project_object(saved_naming_project, notify = FALSE)
  check(isTRUE(loaded$ok), paste("Naming project restore failed:", loaded$message))
  check(selected_assay_id() == saved_naming_project$ui_state$assay_editor_select, "Restored editor selected the wrong assay profile.")
  session$flushReact()
  check(all(rv$assay_profiles$Name_Source == "assay"), "Project round-trip lost naming choices.")
  session$setInputs(target = "LSU", assay_name = "Beta Assay", assay_name_source = "locus",
                    sequencing_primer = "Forward", expected_amplicon_len = 650L, absolute_max_base_index = 680L)
  publish_editor(session, input, selected_assay_id(), assay_editor_revision())
  check(identical(rv$read_assignments$Name_Label, c("Alpha Assay", "LSU")), "Restored profile toggle changed another assay.")
})
cat("Assay locus editor and assignment server integration tests passed.\n")

