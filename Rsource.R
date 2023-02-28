################################################################################
# pre-processing

normScaling = function(dt) {
  #1st: perform median centering:
  y <- apply(dt, 2, function(x){x-median(x)})
  #2nd: unit variance normalization (here *median*=0, SD=1)
  scale(y, center=FALSE)
}

################################################################################

train_model = function(x, y, method="randomForest",ntree=NULL,fold=NULL,minlambda1=1, maxlambda1=10,minlambda2=0.1, maxlambda2=1000,kernel="radial", cost=10){
  if (method=="randomForest"){
    if (is.null(ntree)) model =randomForest(x, y) else model =randomForest(x, y,ntree=ntree)
  }
  
  if (method=="svmLinear"){
    model = svm(x, y,kernel = "linear", cost=cost)
  }

  if (method=="svmRadial"){
    model = svm(x, y,kernel = "radial", cost=cost)
  }
  
  if(method == "xgboost"){
    model = xgboost(x, y, eta=0.05, nrounds = 200)
  }
  
  if(method == "elastic"){
    model = cv.glmnet(x, y, alpha = 0.1)
  }

  if (method=="lassoL1"){
    if (is.null(fold)) reg1 = optL1(y, penalized= x,trace=FALSE) else reg1 = optL1(y, penalized= x,fold=fold,trace=FALSE)
    model = reg1$fullfit
  }
  
  if (method=="lassoL2"){
    if (is.null(fold)) reg2 = optL2(y, penalized= x,trace=FALSE) else reg2 = optL2(y, penalized= x,fold=fold,trace=FALSE)
    model = reg2$fullfit
  }

  return(model) 
}


predict_model = function(model,x, method="randomForest"){
  if (method=="randomForest"){
    pred=predict(model,x)
  }
  
  if (method=="svmLinear"){
    pred=predict(model,x)
  }

  if (method=="svmRadial"){
    pred=predict(model,x)
  }
  
  if(method == "xgboost"){
    pred=predict(model,x)
  }
  
  if(method == "elastic"){
    pred=predict(model,x)
    pred=pred[,1]
  }

  if (method=="lassoL1"){
    if (nrow(x)>1) pred=predict(model,x)[,1] else pred=predict(model,x)[1]
  }
  
  if (method=="lassoL2"){
    if (nrow(x)>1) pred=predict(model,x)[,1] else pred=predict(model,x)[1]
  }

  return(pred) 
}

################################################################################
# feature generation

get_variance_genes = function(x, n_genes=50,rmCor=TRUE){
  gene_vars = apply(x,1,var)
  gene_vars=sort(gene_vars,decreasing = TRUE)
  pick=order(gene_vars,decreasing = TRUE)[1:n_genes]
  if (rmCor){
    x1=x[pick,]
    x2=cor(t(x1),method="spearman")
    x3=abs(x2)
    x3=ifelse(x3==1,0,x3)
    x4=x3 > 0.9
    x5=apply(x4,1,which)    
    if (length(x5)>0){
      rmID=NULL
      for (i in 1:length(x5)) if  (i>min(x5[[i]])) rmID=c(rmID,i)
      if (length(rmID)>0) pick=pick[-rmID]
    }
  }

  var_genes=rownames(x)[pick]
  return(var_genes)
}

get_mutation_genes = function(mut, rnaseq_geneMap){
  mut_genes=rownames(mut)
  mut_genes_ens=rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% mut_genes]
  return(mut_genes_ens)
}

get_pearson_genes = function(x, y, n_genes=50,rmCor=TRUE){
  correlation =apply(x,1,function(z) cor(z,y,method="pearson"))
  pick=order(abs(correlation),decreasing = TRUE)[1:n_genes]
  if (rmCor){
    x1=x[pick,]
    x2=cor(t(x1),method="spearman")
    x3=abs(x2)
    x3=ifelse(x3==1,0,x3)
    x4=x3 > 0.9
    x5=apply(x4,1,which)
    if (length(x5)>0){
      rmID=NULL
      for (i in 1:length(x5)) if (i>min(x5[[i]])) rmID=c(rmID,i)
      if (length(rmID)>0) pick=pick[-rmID]
    }
  }

  pearson_genes=rownames(x)[pick]
  return(pearson_genes)
}

get_spearman_genes = function(x, y, n_genes=50,rmCor=TRUE){
  correlation =apply(x,1,function(z) cor(z,y,method="spearman"))
  #correlation=sort(abs(correlation),decreasing=TRUE)
  pick=order(abs(correlation),decreasing = TRUE)[1:n_genes]
  if (rmCor){
    x1=x[pick,]
    x2=cor(t(x1),method="spearman")
    x3=abs(x2)
    x3=ifelse(x3==1,0,x3)
    x4=x3 > 0.9
    x5=apply(x4,1,which)
    if (length(x5)>0){
      rmID=NULL
      for (i in 1:length(x5)) if  (i>min(x5[[i]])) rmID=c(rmID,i)
      if (length(rmID)>0) pick=pick[-rmID]
    }
  }
  spearman_genes=rownames(x)[pick]
  return(spearman_genes)
}

get_target_genes = function(drug, target, rnaseq_geneMap){
  tar_genes = target[[drug]]
  tar_genes = unique(tar_genes[tar_genes %in% rnaseq_geneMap$Symbol])
  tar_genes = rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% tar_genes]
  return(tar_genes)
}

get_target_genes_UD = function(x, y, drug, target, all_targets, rnaseq_geneMap, n_genes = 200){
  tar_genes = target[[drug]]
  tar_genes = unique(tar_genes[tar_genes %in% rnaseq_geneMap$Symbol])
  tar_genes = rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% tar_genes]
  
  tar_genes_UD = all_targets[[drug]]
  tar_genes_UD = unique(tar_genes_UD[tar_genes_UD %in% rnaseq_geneMap$Symbol])
  tar_genes_UD = rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% tar_genes_UD]

  tar_genes_UD = unique(c(tar_genes_UD,tar_genes))
  if (length(tar_genes_UD) > n_genes){
    x=x[rownames(x) %in% tar_genes_UD,]
    tar_genes_UD=get_spearman_genes(x,y,n_genes)
  } 
  return(tar_genes_UD)
}

#all relevant AML genes: PMID:27276561
get_AML_genes = function(rnaseq_geneMap){
  AML_genes=c("RUNX1T1","RUNX1","PML","RARA","MYH11","CBFB","MLL","NPM1","IDH2",
             "CEBPA","ASXL1","SRSF2","TP53","FLT3","DEK","NUP214","GATA2","MECOM",
             "NRAS","DNMT3A","TET2","PTPN11","SRSF2","STAG2","KIT","WT1","ETV6",
             "PHF6","SF3B1","RAD21","MYC","ZRSR2")
  AML_genes=rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% AML_genes]
  return(AML_genes)
}

#all relevant AML genes - PMID:18716133
get_sig86probes = function(rnaseq_geneMap){
  sig86probes=c("SOCS2","FHL1","NPDC1","TM4SF1","HOPX","C9orf58","BCAT1","RAB34","MYO5C","KIAA0125","FAM92A1","ABI2","DOCK1","GUCY1A3","SHANK3","GPR56","DAPK1","CD109","SCHIP1","PRTFDC1","MIRN155","ATP8B2","RUNX1","TMEM163","C10orf128","COL6A1","MARVELD1","COL24A1","MAST4","TCF4","MAP1A","ARMCX1","FAM30A","WBP5","ACP6","TSC22D1","SYNJ2","LIMS3","LOC440995","GPSM1","MRC1","BCL11A","RPL35A","LAPTM4B","GOLGA8A","MSI2","SPARC","NGFRAP1","ZBTB8","HIST1H2AD","IL23A","PHGDH","MXRA7","RAB13","SCD","NPL","KIAA0922","HBG2","LGALS3","TESC","SLC25A37")
  sig86probes=rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% sig86probes]
  return(sig86probes)
}


#FLT3 signature genes - PMID:18309032 
get_flt3_genes = function(rnaseq_geneMap){
  flt3sig = c('CFHR1','CFH','SOCS2','STON2','LAPTM4B','DPPA4','SCHIP1',
              'PBX3','MYO1B','DAPK1','FAM38B','PDE4B','IL2RA','GOLGA8A',
              'TRH','FLJ37658','HOXB5','CS','HOXB2','GNG2','MYO1B')
  flt3sig_genes=rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% flt3sig]
  return(flt3sig_genes)
}

#TP53 signature genes - PMID:21248301
get_tp53_genes = function(rnaseq_geneMap){
  tp53sig = c("TP53","FLJ12505","GBP5","MYBL2","DHRS2","BRRN1","CHAD","SCGB3A1","DACH1","CDCA8","NPY1R","CYBRD1","DNALI1","LAF4","NY-BR-1","MYBL","CACNG4","CY-BRD1","LRP2","TFF1","TFF3","STC2","AGR2")
  tp53sig_genes=rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% tp53sig]
  return(tp53sig_genes)
}

get_subtype_genes = function(index = 0, rnaseq_geneMap){
  pml =c("PTPRG","HRASLS5","GABRE","NKX3-2","LOXL4","IL17RE",
         "SIX3","MEG3","EBF3","SLC24A3","MPO","ADAMTS2","SURF4",
         "STXBP1","UGT3A2")  #PML-RARA
  runx =c("RUNX1T1","POTEA","ARC","C1QTNF5","SIPA1L2","POU4F1",
          "LINC00958","PNMT","WFDC1","BAIAP3","MLLT4","PTGIR",
          "PALM","GRK5","ELAVL4")  #RUNX1
  cbfb = c("RPS6KL1","MYH11","GNG2","XPNPEP2","PRCD","CACNB4",
           "ZNF683","GRASP","TRIM71","ZFHX3","RPS6KA2","PAPSS2",
           "MTSS1","MGC12916","SULF2") #CBFB
  npm1 =c("HOXA-AS3","HOXB-AS3","HOXB4","HOXA3","HOXB3","PBX3",
          "HOXA5","HOXA6","HOXA4","HOXB6","HOXB2","HOXB5","HOXA9",
          "SDSL","HOXA7") #NPM1
  mll =c("SDHB","COX5B","PPP1CC","COX6A1","GCHFR","NDUFA6","COX5A",
         "SEC11C","EBP","LGALS1","HINT2","LOC442132","COX7B","COX6B1",
         "NDUFA13" ) #MLL
  ceb=c("DLC1","APBA1","PTX4","CETP","CRYGD","CLEC2L","LOC102724279",
        "HPGDS","SHD","UGT2B28","ZNF385C","UMODL1","KLHL33","SVOPL",
        "GPA33") #CEBPA_Biallelic
  idh2=c("ADRA2C","ZIM2","PLSCR4","GAL3ST1","ZC2HC1A","HIVEP2","NKX2-3",
         "TPTEP1","FXYD5","CPED1","NUDT10","CASC10","MFAP2") #IDH2_R172
  spli=c("L3MBTL4","PTK2","SETBP1","TRPS1","PLCL2","H1F0","PLXNC1",
         "MAP3K14","NCOA2","AKAP7","SLC37A3","HLA-DRA","SLC41A1",
         "MZT2A","TCF4") #splice
  inv3=c("LRRC16A","SLC9A9","LY6G6E","ATP8B2","APOL3","VEGFC","GBP1",
         "CASS4","SH3BP4","SELP","GP9","ASAP2","AFAP1L2","DPF1","PF4") #inv3
  p53c3=c("SH3BP5","COL5A1","CHIT1","COL18A1","SYNE1","GADD45A",
          "SPINT2","ITPR3","DGAT1","SVIP","CA3","ZIK1") #p53C
  subtype_genes = list(pml, runx, cbfb, npm1, mll, ceb, idh2, spli, inv3, p53c3)
 
  if(index == 0){
    subtype_genes = unlist(subtype_genes, recursive = FALSE)
  } else{
    subtype_genes = subtype_genes[[index]]
  } 
  
  subtype_genes = rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% subtype_genes]  
  return(subtype_genes)
}

get_literature_mutation_genes = function(rnaseq_geneMap){
  #mut.gene
  mutAMLGenes=c("RUNX1T1","RUNX1","PML","RARA","MYH11","CBFB","MLL","NPM1","IDH2","CEBPA","ASXL1","SRSF2","TP53","FLT3","DEK","NUP214","GATA2","MECOM","NRAS","DNMT3A","TET2","PTPN11","SRSF2","STAG2","KIT","WT1","ETV6","PHF6","SF3B1","RAD21","MYC","ZRSR2")
  mutGenes=mutAMLGenes
  mutGenes=c("TP53","ENAH","U2AF1","DSCAM","PHF6","AGTPBP1","SRSF2")
  mutGenes=unique(c(mutGenes,mutAMLGenes))
  mut_genes = rnaseq_geneMap$Gene[rnaseq_geneMap$Symbol %in% mutGenes] 
  return(mut_genes)
}
