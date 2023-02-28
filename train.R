train_per_drug <- function(drug, gex, mut, drug_response, pas, ic50,
                           neaMat, aml_pas, target, all_targets,
                           rnaseq_geneMap, method = "svmRadial", ntree = NULL) {
  train_drug_response <- drug_response[drug_response$inhibitor %in% drug, ]

  if (ic50) {
    y <- train_drug_response$ic50
  } else {
    y <- train_drug_response$auc
  }
  names(y) <- train_drug_response$lab_id

  x <- gex[, colnames(gex) %in% train_drug_response$lab_id]
  y <- y[colnames(x)]

  n_genes              <- 100
  variance_genes       <- get_variance_genes(x, n_genes)
  pearson_genes        <- get_pearson_genes(x, y, 2 * n_genes)
  spearman_genes       <- get_spearman_genes(x, y, 2 * n_genes)
  target_genes         <- get_target_genes(drug, target, rnaseq_geneMap)
  target_genes_UD      <- get_target_genes_UD(x, y, drug, target, 
                                              all_targets, rnaseq_geneMap)
  AML_genes            <- get_AML_genes(rnaseq_geneMap)
  mutation_genes       <- get_mutation_genes(mut, rnaseq_geneMap)
  liter_mutation_genes <- get_literature_mutation_genes(rnaseq_geneMap)
  flt3_genes           <- get_flt3_genes(rnaseq_geneMap)

  # 3, 5, 9
  subtype_genes <- NULL
  if (drug %in% rownames(neaMat)) {
    pos <- which(neaMat[drug, ] >= 3)
    temp <- NULL
    if (length(pos) > 0) {
      for (p in pos) temp <- c(temp, get_subtype_genes(p, rnaseq_geneMap))
    }
    subtype_genes <- temp
  }

  tp53_genes  <- get_tp53_genes(rnaseq_geneMap)
  sig86probes <- get_sig86probes(rnaseq_geneMap)

  pick <- unique(c(
    variance_genes, pearson_genes, spearman_genes, mutation_genes,
    target_genes_UD, target_genes, tp53_genes, AML_genes,
    subtype_genes, sig86probes
  ))

  pick           <- intersect(pick, rownames(gex))
  x              <- x[pick, ]
  selected_genes <- pick
  x2             <- x

  # PAS_gex
  pas_drugs <- names(pas$gex$final)
  if (drug %in% pas_drugs) {
    p0 <- pas$gex$final[[drug]]

    # consider druggable CSP
    pick_gdsc <- aml_pas$drug_name %in% drug
    if (sum(pick_gdsc) > 0) {
      pas1         <- aml_pas[pick_gdsc, ]
      pathway_gdsc <- intersect(rownames(p0), pas1$pathway)
      pick1_gdsc   <- which(pathway_gdsc %in% rownames(p0))
      if (length(pick1_gdsc) > 0) {
        p0 <- p0[pick1_gdsc, ]
      }
    }
    rownames(p0) <- paste("all", rownames(p0), sep = "_")

    # PAS_mut
    p3           <- pas$mut$up[[drug]]
    rownames(p3) <- paste("onupMut", rownames(p3), sep = "_")
    p4           <- pas$mut$down[[drug]]
    rownames(p4) <- paste("downMut", rownames(p4), sep = "_")

    pas_gex <- rbind(p3, p4) # standard
    if (sum(pick_gdsc) > 0) {
      pas_gex <- rbind(p0, p3, p4)
    }

    pas_gex      <- pas_gex[, colnames(x)] # keep only ones matched with x

    pick1        <- rowSums(pas_gex != 0) > 0
    pas_gex      <- pas_gex[pick1, drop = FALSE, ]
    selected_pas <- rownames(pas_gex)

    x2           <- rbind(x2, pas_gex)
  }

  final_genes <- rownames(x2)
  x2 <- t(x2)

  model <- train_model(x2, y, method = method, cost = 10, ntree = ntree)

  return(list(model, selected_genes, selected_pas, final_genes))
}

train_ensemble <- function(drug, gex, drug_response, pas, ic50,
                           models, selected_genes, selected_pas, final_genes,
                           method_1 = "svmRadial", method_2 = "svmRadial") {
  cat(drug, " ")
  train_drug_response <- drug_response[drug_response$inhibitor %in% drug, ]

  if (ic50) {
    y <- train_drug_response$ic50
  } else {
    y <- train_drug_response$auc
  }
  names(y) <- train_drug_response$lab_id

  x <- gex[, colnames(gex) %in% train_drug_response$lab_id]
  y <- y[colnames(x)]

  preds <- NULL
  all_drugs <- as.character(unique(drug_response$inhibitor))
  pas_drugs <- names(pas$gex$final)

  for (inhibitor in all_drugs) {
    x1 <- gex[selected_genes[[inhibitor]], , drop = FALSE]

    # PAS_gex
    if (inhibitor %in% pas_drugs) {
      p0           <- pas$gex$final[[inhibitor]]
      rownames(p0) <- paste("all", rownames(p0), sep = "_")
      p1           <- pas$gex$up[[inhibitor]]
      rownames(p1) <- paste("onup", rownames(p1), sep = "_")
      p2           <- pas$gex$down[[inhibitor]]
      rownames(p2) <- paste("down", rownames(p2), sep = "_")

      # PAS_mut
      p3           <- pas$mut$up[[inhibitor]]
      rownames(p3) <- paste("onupMut", rownames(p3), sep = "_")
      p4           <- pas$mut$down[[inhibitor]]
      rownames(p4) <- paste("downMut", rownames(p4), sep = "_")
      p5           <- pas$mut$final[[inhibitor]]
      rownames(p5) <- paste("allMut", rownames(p5), sep = "_")

      pas_gex      <- rbind(p0, p1, p2, p3, p4, p5)
      pas_gex      <- pas_gex[, colnames(x1)] # keep only ones matched with x

      pas_gex      <- pas_gex[selected_pas[[inhibitor]], , drop = FALSE]
      x1           <- rbind(x1, pas_gex)
    }

    x1    <- t(x1)
    x1    <- x1[, final_genes[[inhibitor]]]
    # if (!is.null(extra_features)) {
    #   extra <- extra_features[rownames(x1), ]
    #   x1 <- cbind(x1, extra)
    # }
    model <- models[[inhibitor]]

    pred <- predict_model(model, x1, method = method_1)

    preds <- rbind(preds, pred)
  }
  
  rownames(preds) <- all_drugs
  pred_auc        <- preds[, names(y)]

  pred_auc <- t(pred_auc)


  model <- train_model(pred_auc, y, method_2,
                       minlambda2 = 0.1, maxlambda2 = 1e6, cost = 1)

  return(model)
}

train <- function(data, pas, ic50 = TRUE, 
                  method_1 = "svmRadial", method_2 = "svmRadial", n_cores = 2) {
  load("data/aml_pas.RData")
  load("data/NEA_drugs_subtypes.RData")
  load("data/Target_UpDownStream_BeatAML.RData")
  load("data/rnaseq_geneMap.RData")

  all_targets     <- Map(union, tgt_down_all, tgt_up_all)

  gex             <- data$gex
  drug_response   <- data$drug_response
  mut             <- data$mut

  all_drugs       <- as.character(unique(drug_response$inhibitor))


  match_id        <- match(colnames(gex), colnames(mut))
  missing_samples <- which(is.na(match_id))
  if (length(missing_samples) > 0) {
    ex_dat           <- matrix(0,
                               ncol = length(missing_samples),
                               nrow = nrow(mut))
    colnames(ex_dat) <- colnames(gex)[missing_samples]
    mut              <- cbind(mut, ex_dat)
    match_id         <- match(colnames(gex), colnames(mut))
    mut              <- mut[, match_id]
  }

  pick         <- rowSums(mut > 0)
  mut          <- mut[pick > 2, ] # more than 1%
  mut[mut > 1] <- 1

  # normalise gene expression
  gex  <- gex - min(gex) # shift the coordinate to zero
  pick <- apply(gex, 1, function(r) sum(r < 1e-6))
  pick <- pick > 0.95 * ncol(gex)
  gex  <- gex[!pick, ]
  gex  <- normScaling(gex)

  ################################################################
  ### building train models

  trained_model <- mclapply(all_drugs,
    function(drug, gex, mut, drug_response, pas, ic50,
             neaMat, aml_pas, target, all_targets, 
             rnaseq_geneMap, method, ntree) {
      train_per_drug(
        drug, gex, mut, drug_response, pas, ic50,
        neaMat, aml_pas, target, all_targets, rnaseq_geneMap, method, ntree
      )
    },
    gex = gex, mut = mut, drug_response = drug_response,
    pas = pas, ic50 = ic50, neaMat = neaMat,
    aml_pas = aml_pas, target = target,
    all_targets = all_targets,
    rnaseq_geneMap = rnaseq_geneMap,
    method = method_1,
    ntree = NULL,
    mc.preschedule = FALSE, mc.cores = n_cores
  )

  models         <- lapply(trained_model, function(l) l[[1]])
  selected_genes <- lapply(trained_model, function(l) l[[2]])
  selected_pas   <- lapply(trained_model, function(l) l[[3]])
  selected_final <- lapply(trained_model, function(l) l[[4]])


  names(models)       <- names(selected_genes) <- all_drugs
  names(selected_pas) <- names(selected_final) <- all_drugs

  models2 <- mclapply(all_drugs,
    function(drug, gex, drug_response, pas, ic50,
             models, selected_genes,
             selected_pas, final_genes, method_1, method_2) {
      train_ensemble(
        drug, gex, drug_response, pas,
        ic50, models, selected_genes,
        selected_pas, final_genes, method_1, method_2
      )
    },
    gex = gex, drug_response = drug_response, pas = pas,
    ic50 = ic50, models = models,
    selected_genes = selected_genes,
    selected_pas = selected_pas,
    final_genes = selected_final,
    method_1 = method_1, method_2 = method_2,
    mc.preschedule = FALSE, mc.cores = n_cores
  )

  names(models2) <- all_drugs

  return(list(
    "models" = models, "models2" = models2,
    "selected_genes" = selected_genes,
    "selected_pas" = selected_pas,
    "selected_final" = selected_final
  ))
}
