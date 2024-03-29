# THE edgeR package
run_edgeRQLF_withCont <- function(L) {
  library(limma)
  library(edgeR)
  message("edgeRQLF")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList( L$count, genes = data.frame( ENTREZID = rownames ( L$count ) ) )
    dge <- calcNormFactors(dge)
    design <- model.matrix(~0 + L$Treat)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    for (pp in colnames(L$cont.matrix)) {
      qlf <- glmLRT(fit, contrast = L$cont.matrix[,pp])
      tt <- topTags(qlf, n = nrow(dge$counts))
    }
  })
  write.csv(tt, 'tt.csv')
  plotBCV(dge)
  plotQLDisp(fit)
  hist(tt$table$PValue, 50)
  hist(tt$table$FDR, 50)
  limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
  plotSmear(qlf)

  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$table$PValue,
                       padj = tt$table$FDR,
                       row.names = rownames(tt$table)))
}

scEISA <- function(exon_list, intron_list) {
  message("scEISA")
  session_info <- sessionInfo()
  library(dplyr)
  # the exon_list and intron list should contain condition and count matrix
  condit = c(exon_list$condt, intron_list$condt) #combine the conditions
  treatment <- unique(L$condt)
  {if (length(treatment) !=2)
    + stop ('the function only takes 2 treatments')
    }
  len = length(conditions)/4
  #len for condition 1 must equal condition 2
  conditions = c(
    rep('cont1', len),
    c(rep('cont2', len)),
    rep('cont1', len),
    c(rep('cont2', len)))
  cnt <- cbind(Ex=as.data.frame(cntEx[genes.sel,]), In=cntIn[genes.sel,])
  genes_no_overlap <- rownames(cntEx)
  genes.sel <-genes_no_overlap
  cntEx = exon_list$count
  cntIn = intron_list$count
  cnt = cbind(Ex=as.data.frame(cntEx[genes.sel,]), In=cntIn[genes.sel,])
  dge = DGEList(counts = cnt, genes = data.frame(ENTREZID=rownames(cnt)))
  dge = calcNormFactors(dge)
  libsize_sel = data.frame(lib=colnames(cnt)) %>%
    mutate(factorRegion=factor(rep(c("ex","in"),each=ncol(cntEx)), levels=c("in", "ex"))) %>%
    mutate(factorCondition = factor(conditions, levels=unique(conditions))) %>% #changed
    mutate(genotype = paste(factorRegion,"_", factorCondition,sep=""))
  Treat <- factor(libsize_sel$genotype)
  design <- model.matrix(~0+Treat)
  colnames(design) <- levels(Treat)
  cont.matrix = cont.matrix <- makeContrasts(
    in_cont1_cont2 = in_cont1 - in_cont2,
    ex_cont1_cont2 = ex_cont1 - ex_cont2,
    exonMin_intron = (ex_cont1 - in_cont1) - (ex_cont2- in_cont2),
    levels=design)  # build contrast matrix
  L <- list(count = cnt, condt = conditions, cont.matrix = cont.matrix, Treat = Treat)
  run_edgeRQLF_withCont(L)
}

