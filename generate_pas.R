fn.gsum = function(gene, xdat){
  ## start with expression-sum of a set of genes, 
  pick = gene %in% rownames(xdat)
  gene = gene[pick]
  if (length(gene)>1) gsum = apply(xdat[gene,],2, sum)
  if (length(gene)==1) gsum= xdat[gene,]
  if (length(gene)==0) gsum = rep(0, ncol(xdat))
  return(gsum)
}

get_pas = function(drug, exp, up_gene_set, down_gene_set){
  down.gex = onup.gex = pas.gex = NULL
  if (drug %in% names(up_gene_set)){
    ingene = up_gene_set[[drug]]
    insum.gex = t(sapply(ingene, fn.gsum, exp))
    insum.gex[is.na(insum.gex)] = 0
    onup.gex = matrix(unlist(insum.gex), ncol=ncol(insum.gex))

    ingene = down_gene_set[[drug]]
    insum.gex = t(sapply(ingene, fn.gsum, exp))
    insum.gex[is.na(insum.gex)] = 0
    down.gex = matrix(unlist(insum.gex), ncol=ncol(insum.gex))

    pas.gex = onup.gex - down.gex
    colnames(pas.gex) = colnames(exp)
    colnames(onup.gex) = colnames(exp)
    colnames(down.gex) = colnames(exp)
    rownames(pas.gex) = names(ingene)
    rownames(onup.gex) = names(ingene)
    rownames(down.gex) = names(ingene)

  }

  return(list("final_pas" = pas.gex, "up_pas" = onup.gex,
              "down_pas" = down.gex))
}

generate_pas = function(data, gene_set, n_cores = 2){
  up_gene_set = gene_set$up
  down_gene_set = gene_set$down

  load("data/all_drugs.RData")
  load("data/rnaseq_geneMap.RData")


  #pas for gene expression
  pas.res = mclapply(all_drugs, 
		                 function(drug, exp, up, down) get_pas(drug, exp, up, down),
		                 exp = data$gex,
                     up = up_gene_set,
                     down = down_gene_set,
                     mc.preschedule = FALSE, mc.cores = n_cores)

  pas.gex.all = lapply(pas.res, function(l) l[[1]])
  pas.gex.onup = lapply(pas.res, function(l) l[[2]])
  pas.gex.down = lapply(pas.res, function(l) l[[3]])
  names(pas.gex.down) = names(pas.gex.onup) = names(pas.gex.all) = all_drugs

  pas_gex = list("final" = pas.gex.all, 
                 "up" = pas.gex.onup,
                 "down" = pas.gex.down)

  #pas for mutation
  drimut = data$mut 


  drimut[is.na(drimut)] = 0  ## set NA to 0
  sumdri = rowSums(drimut>0)
  pick  = sumdri > 2
  drimut = drimut[pick,]
  pick=rnaseq_geneMap$Symbol %in% rownames(drimut)
  ens_names=rnaseq_geneMap[pick,]
  drimut2=drimut[ens_names$Symbol,]
  rownames(drimut2)=ens_names$Gene
  drimut=drimut2
  drimut=drimut[rownames(drimut) %in% rownames(data$gex),]
  drimut=as.matrix(drimut)

  mut_exp = data$gex[rownames(drimut), ] * drimut
  pas.res = mclapply(all_drugs, 
		                 function(drug, exp, up, down) get_pas(drug, exp, up, down),
		                 exp = mut_exp,
                     up = up_gene_set,
                     down = down_gene_set,
                     mc.preschedule = FALSE, mc.cores = n_cores)
  pas.mut.all = lapply(pas.res, function(l) l[[1]])
  pas.mut.onup = lapply(pas.res, function(l) l[[2]])
  pas.mut.down = lapply(pas.res, function(l) l[[3]])
  names(pas.mut.down) = names(pas.mut.onup) = names(pas.mut.all) = all_drugs
  pas_mut = list("final" = pas.mut.all,
                 "up" = pas.mut.onup,
                 "down" = pas.mut.down)

  return(list("gex" = pas_gex, "mut" = pas_mut))
}
