# ---------------------------------------------------------------------------------------
# Script: LineageMCMC.R
# Purpose: Per-lineage BEACON (Bayesian correlation) to quantify GED/PED within lineages
# Inputs:
#   - DepMap sample info: DepMap_data/sample_info_22Q2.csv
#   - mRNA: CCLE_expression_22Q2.csv.gz | Protein: mmc2.xlsx
#   - CRISPR dependency: DepMap_data/CRISPR_gene_effect_22Q2.csv.gz
# Outputs:
#   - One Excel per lineage: Table.<DATA>.dependency.Bayesian.lineage.<LINEAGE>.xlsx
# Reproducibility:
#   - Fixed R and JAGS seeds; posterior mean of rho is exported as column 'rho'
# Download links:
#   - sample_info.csv file     @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=35020903
#   - CCLE_expression.csv file (further gzipped)     @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34989919
#   - Supplementary data (Table S2: normalized protein expressions) from Nusinow et al. paper (doi.org/10.1016/j.cell.2019.12.023)     @ https://www.cell.com/cms/10.1016/j.cell.2019.12.023/attachment/3709dedc-3a01-4e1d-ab4c-82597295c5d2
#   - CRISPR_gene_effect.csv file (further gzipped)     @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34990036
# Expected runtime (8 vCPU/16 GB RAM; n.iter=500):
#   - ~5–15 min per lineage (varies with n cell lines); full run ~2–4 hrs
# ---------------------------------------------------------------------------------------
# ml R jags/4.3.0
# R
rm(list = ls(all.names = TRUE))
options(java.parameters = "-Xmx8000m")
library(openxlsx)
library(rjags)
set.seed(12345)  # fixed R seed

# ---------------- Parameters ----------------
n.adapt  = 200   # JAGS adaptation steps
n.update = 200   # burn-in before sampling
n.iter   = 1000   # MCMC iterations per chain

reproduce.results = TRUE
recalculate.FDR   = TRUE
out               = TRUE
continue.from     = NULL   # e.g. resume index; usually NULL

data      = 'mRNA'         # 'mRNA' | 'Protein'
panel     = ''             # e.g. '.druggable', 'DrugBank.Antibody', etc.
cell.type = 'All'          # 'All' | 'Primary' | 'Metastasis'

num.top.genes = 10
thr.FDR       = .05
# --------------------------------------------

# Output dir (portable)
dir_path <- file.path('..','out',
  paste0('jags.nadapt', n.adapt, '.update', n.update, '.mcmc', n.iter, '.simulation_SD_22Q2'))
if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

# Output path helper
out_lineage_xlsx <- function(data, lineage, panel, cell.type, continue.from=NULL) {
  file.path(
    dir_path,
    paste0('Table.', data, '.dependency.Bayesian.lineage.', lineage,
           panel, sub('^\\.All$','', paste0('.', cell.type)), continue.from, '.xlsx')
  )
}

# ---------------- Inputs (relative) ----------------
depmap_info_path <- file.path('..','..','..','..','Huang_lab_data','DepMap_data','sample_info_22Q2.csv') #| sample_info.csv file @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=35020903
mrna_path        <- file.path('..','..','..','..','Huang_lab_data','DepMap_data','CCLE_expression_22Q2.csv.gz') #log2(TPM+1) #| CCLE_expression.csv file (further gzipped) @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34989919
protein_xlsx     <- file.path('..','..','..','..','Huang_lab_data','QuantProtCCLE_Nusinow_Cell2020','mmc2.xlsx') #| Supplementary data (Table S2: normalized protein expressions) from Nusinow et al. paper (doi.org/10.1016/j.cell.2019.12.023) @ https://www.cell.com/cms/10.1016/j.cell.2019.12.023/attachment/3709dedc-3a01-4e1d-ab4c-82597295c5d2 
crispr_path      <- file.path('..','..','..','..','Huang_lab_data','DepMap_data','CRISPR_gene_effect_22Q2.csv.gz') #| CRISPR_gene_effect.csv file (further gzipped) @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34990036
# ---------------------------------------------------

# Run log
.run_start <- Sys.time()
.run_log <- list(
  script = "LineageMCMC.R",
  start_time = as.character(.run_start),
  params = list(n.adapt=n.adapt, n.update=n.update, n.iter=n.iter,
                data=data, panel=panel, cell.type=cell.type,
                thr.FDR=thr.FDR)
)

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#   intsect
# 
#   intersect two sets of which at least one contains only unique elements
#
intsect = function(foo, bar, map.to = 2) {
  comm = vector()
  ifoo = vector()
  ibar = vector()
  for (f in seq_along(foo)) {
    check = bar %in% foo[f]
    if (any(check)) {
      loc = which(check)
      comm = c(comm, rep(foo[f], length(loc)))
      ifoo = c(ifoo, rep(f, length(loc)))
      ibar = c(ibar, loc)
    }
    
  }
  if (map.to == 1) {
    s = sort(ifoo, decreasing = F, index.return = T, na.last = T)
  } else if (map.to == 2) {
    s = sort(ibar, decreasing = F, index.return = T, na.last = T)
  }
  comm = comm[s$ix]
  ifoo = ifoo[s$ix]
  ibar = ibar[s$ix]
  return(list(comm, ifoo, ibar))
}
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#   regexprTab
#
#   a table of the presence of the regular expressions (row) in the data (col)
#
regexprTab = function(expressions = c('foo','bar','^f','n$'), data = c('food','barn')){
  expressions = unique(expressions)
  contains = data.frame(matrix(F, length(expressions), length(data)))
  rownames(contains) = expressions; colnames(contains) = data
  for (e in expressions){
    contains[e,] = regexpr(e,data) > 0
  }
  return(contains)
}
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
# LOAD DATA
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

sam.dep = read.csv(depmap_info_path)
table(sam.dep$primary_or_metastasis)
if (cell.type!='All'){
  sam.dep = sam.dep[sam.dep$primary_or_metastasis==cell.type,]
}

dep.dat.map = sam.dep[,c(1,3)]
dep.dat.map$Tissue_Type = as.character(sam.dep$CCLE_Name)
sam.dep$Tissue_Type = as.character(sam.dep$CCLE_Name)
i=1
for (i in seq_len(nrow(dep.dat.map))){
  # dep.dat.map$Tissue_Type[i] = sub(paste0('^',dep.dat.map$stripped_cell_line_name[i], '\\_'),'',dep.dat.map$Tissue_Type[i])
  # sam.dep$Tissue_Type[i] =     sub(paste0('^',sam.dep$stripped_cell_line_name[i], '\\_'),'',sam.dep$Tissue_Type[i])
  # Derive tissue lineage label from CCLE name (strip leading ACH code)
  dep.dat.map$Tissue_Type[i] = sub('^[a-zA-Z0-9]*\\_','',dep.dat.map$Tissue_Type[i])
  sam.dep$Tissue_Type[i] = sub('^[a-zA-Z0-9]*\\_','',sam.dep$Tissue_Type[i])
}
# dep.dat.map0 = read.csv('~/Downloads/df_cell_line_mapping_20Q1.csv')


if (data == 'Protein') {
  # data exp
  exp.dat = openxlsx::read.xlsx(protein_xlsx, sheet = 'Normalized Protein Expression')
  exp.dat = cbind(Gene_Symbol = exp.dat$Gene_Symbol, exp.dat[,-(1:50)])
  col.drp = which(colSums(regexprTab(expressions = c('^Protein', '^Description', '^Group', '^Uniprot', 'Peptides$', '^Column'), data = colnames(exp.dat))) > 0)
  exp.dat = exp.dat[,-col.drp]
  # Filter missing value genes
  # exp.dat = exp.dat[-which(rowMeans(is.na(exp.dat)) > .9),]
  tab.gen = table(exp.dat$Gene_Symbol)
  # Protein duplicates: keep the replicate with fewest missing values per gene
  row.drp = c()
  for (gene in names(tab.gen)[which(tab.gen>1)]){
    idx = which(exp.dat$Gene_Symbol %in% gene) 
    row.drp = c(row.drp, setdiff(idx,idx[order(rowSums(is.na(exp.dat[idx,])))[1]]))
  }
  exp.dat = exp.dat[-row.drp,]
  # Filter missing value cells 
  # exp.dat = exp.dat[, -which(colMeans(is.na(exp.dat)) > .9)]
  # Cell-line duplicate TenPx columns: keep the one with fewest missing values
  tab.cel = table(gsub('\\_TenPx.*$','',colnames(exp.dat)))
  col.drp = c()
  for (cel in names(tab.cel)[which(tab.cel>1)]){
    idx = which(gsub('\\_TenPx.*$','',colnames(exp.dat)) %in% cel) 
    col.drp = c(col.drp, setdiff(idx,idx[order(colSums(is.na(exp.dat[,idx])))[1]]))
  }
  exp.dat = exp.dat[,-col.drp]
  #  
  gen.dat = as.character(exp.dat$Gene_Symbol)
  gen.dat.original = gen.dat
  cel.dat = gsub('\\_.*','',colnames(exp.dat)[-1])
  exp.dat = exp.dat[,-1]
  colnames(exp.dat) = gsub('\\_','\\.',gsub('\\_TenPx.*$','',colnames(exp.dat)))
} else if (data == 'mRNA') {
  exp.dat = read.csv(mrna_path) # log2(TPM+1)
  cel.dat = exp.dat$X
  exp.dat = t(exp.dat)[-1,]
  gen.dat = gsub('\\.\\..*$','',rownames(exp.dat))
  gen.dat.original = gen.dat
  colnames(exp.dat) = as.character(gsub('\\_','\\.',apply(dep.dat.map[match(cel.dat, dep.dat.map$DepMap_ID),-1], 1, function(x) paste(x, collapse = '.'))))
  cel.dat = apply(as.matrix(colnames(exp.dat)), 1, function(x) {strsplit(x, split = '\\.')[[1]][1]})
  # fix the ACH code not mapping to cell lines: discard
  idrop = which(cel.dat=='NA')
  if (length(idrop)>0){
    cel.dat = cel.dat[-idrop]
    exp.dat = exp.dat[,-idrop]
  }
}

# Dependency data
# ccl.dep = read.csv('~/Downloads/Achilles_gene_effect_20Q1.csv.gz')
ccl.dep = read.csv(crispr_path)
#' =================CRISPR_gene_effect_22Q2 (CERES)==========================================
#' CRISPR knockout screens published by Broad’s Achilles and Sanger’s SCORE projects.
#' Negative scores imply "cell growth inhibition and/or death" following gene knockout. 
#' The scores were normalized to ensure that nonessential genes have a median score of 0, 
#' while commonly identified essential genes have a median score of -1.
#' =================CRISPR_gene_effect_22Q2 (CERES)==========================================

map = match(ccl.dep$DepMap_ID, dep.dat.map$DepMap_ID) #***
cel.dep = as.character((dep.dat.map$stripped_cell_line_name[map[!is.na(map)]]))
# cel.dep = as.character((dep.dat.map$cell_line_name[map[!is.na(map)]]))
# typ.dep = as.character(sam.dep$primary_or_metastasis[map[!is.na(map)]])
cel.tis.dep = as.character(unique( 
  apply(dep.dat.map[map[!is.na(map)],-1],1,function(x){gsub('\\_','\\.',paste(x,collapse='.'))}) ))
gen.dep = gsub('\\.\\..*$','',colnames(ccl.dep)[-1])

# if (cell.type!='All'){
#   pic.cel = which(typ.dep==cell.type)
#   ccl.dep = ccl.dep[pic.cel,]
#   cel.dep = cel.dep[pic.cel]
#   cel.tis.dep = cel.tis.dep[pic.cel]
# }

if (panel=='.druggable'){
  # filter for druggable genes
  dru.gen = readRDS('~/Downloads/box/drug.genes.RDS')
  ccl.dep = ccl.dep[, which(gen.dep %in% dru.gen)]
  gen.dep = gen.dep[which(gen.dep %in% dru.gen)]
}  

# frame data: cell.line | tissue type
ccl.dep.nam = setNames(data.frame(t(sapply(cel.tis.dep,       function(x) {c(strsplit(x, '\\.')[[1]][1], paste(strsplit(x, '\\.')[[1]][-1], collapse = '.'))}))), c('Cell','Tissue'))
exp.dat.nam = setNames(data.frame(t(sapply(colnames(exp.dat), function(x) {c(strsplit(x, '\\.')[[1]][1], paste(strsplit(x, '\\.')[[1]][-1], collapse = '.'))}))), c('Cell','Tissue'))

# Filter lineages with <=5 cell lines (insufficient power/stability)
ccl.dep.drp = which(ccl.dep.nam$Tissue %in% names(table(ccl.dep.nam$Tissue)[table(as.character(ccl.dep.nam$Tissue)) <= 5]))
exp.dat.drp = which(exp.dat.nam$Tissue %in% names(table(exp.dat.nam$Tissue)[table(as.character(exp.dat.nam$Tissue)) <= 5]))
# filter
if (length(ccl.dep.drp)>0){
  cel.dep = cel.dep[-ccl.dep.drp]
  cel.tis.dep = cel.tis.dep[-ccl.dep.drp]
  ccl.dep = ccl.dep[-ccl.dep.drp,]
  ccl.dep.nam = ccl.dep.nam[-ccl.dep.drp,]
}
if (length(exp.dat.drp)){
  cel.dat = cel.dat[-exp.dat.drp]
  exp.dat = exp.dat[,-exp.dat.drp]
  exp.dat.nam = exp.dat.nam[-exp.dat.drp,]
}

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
# ALIGN CELL LINES — critical for pairing expression and dependency correctly
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
ali = intsect(cel.dep, cel.dat)
#
cel.dep = cel.dep[ali[[2]]]
cel.tis.dep = cel.tis.dep[ali[[2]]]
ccl.dep = ccl.dep[ali[[2]],]
ccl.dep.nam = ccl.dep.nam[ali[[2]],]
#
cel.dat = cel.dat[ali[[3]]]
exp.dat = exp.dat[,ali[[3]]]
exp.dat.nam = exp.dat.nam[ali[[3]],]

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
# 
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
if (reproduce.results) {
  # the gene's expression is the dependent variable
  dep.var = apply(exp.dat, 2, function(x){as.numeric(x)})
  # the gene's dependency is the independent variable
  indep.var = apply(ccl.dep[,-1], 1, function(x){as.numeric(x)})
} else {
  rm('ccl.dep', 'exp.dat')
}
model_string = '
  model {
    for(i in 1:n) {
      x[i,1:2] ~ dmnorm(mu[], prec[ , ])
    }
 
    # Constructing the covariance matrix and the corresponding precision matrix
    prec[1:2,1:2] <- inverse(cov[,])
    cov[1,1] <- sigma[1] * sigma[1]
    cov[1,2] <- sigma[1] * sigma[2] * rho
    cov[2,1] <- sigma[2] * sigma[1] * rho
    cov[2,2] <- sigma[2] * sigma[2]

    # Uninformative priors on all parameters
    mu[1] ~ dnorm(0, 1)
    mu[2] ~ dnorm(0, 1)
    #
    sigma[1] ~ dunif(0, 100) 
    sigma[2] ~ dunif(0, 100)
    #
    rho ~ dunif(-1, 1)

    # Generate random draws from the estimated bivariate normal distribution
    x_rand ~ dmnorm(mu[], prec[ , ])
  }
'
# For each lineage
# mcm.lin.all = c()
# mcm.lin.all.top = c()
lineage='SOFT.TISSUE'
for (lineage in sort(unique(as.character(ccl.dep.nam$Tissue)))) {#***
# for (lineage in sort(unique(as.character(ccl.dep.nam$Tissue)))[-(1:13)]) {#*** 
  message('\n[BEACON] Processing ', data, ' — ', lineage, ' ...')
  .t0 <- Sys.time()
  lin.loc = which(lineage == ccl.dep.nam$Tissue)
  opath   = out_lineage_xlsx(data, lineage, panel, cell.type, continue.from)
  if (file.exists(opath) & !reproduce.results){
    mcm.lin = xlsx::read.xlsx(opath, sheetIndex = 1)
  } else {
    mcm.lin = c()
    i=1
    for (i in seq_len(nrow(indep.var))){
      gen.loc = which(gen.dat == gen.dep[i])
      if (length(gen.loc) > 0){
        r = tryCatch({
          gg = gen.loc[1]
          for (gg in gen.loc){
            exp.dep.pair = cbind(dep.var[gg,lin.loc], indep.var[i,lin.loc])
            exp.dep.pair = exp.dep.pair[rowSums(is.na(exp.dep.pair)) == 0,]
            rho.pri = cor(exp.dep.pair, method = 'pearson')[2]
            
            dat.lis = list(x = exp.dep.pair, n = nrow(exp.dep.pair))
            ini.lis = list(mu = colMeans(exp.dep.pair, na.rm=T), rho = rho.pri, sigma = apply(exp.dep.pair, 2, function(x){sd(x, na.rm = T)}))

            # Create initialization list with RNG specifications for each chain
            ini.lis_with_seeds <- list(
              c(ini.lis,
                list(.RNG.name = "base::Mersenne-Twister",
                     .RNG.seed = 1)),
              c(ini.lis,
                list(.RNG.name = "base::Mersenne-Twister",
                     .RNG.seed = 2))
            )
            
            jag.mod = jags.model(textConnection(model_string), data = dat.lis, inits = ini.lis_with_seeds, n.adapt = n.adapt, n.chains = 3, quiet = T)
            update(jag.mod, n.update, progress.bar = 'none')
            mcm.sam = coda.samples(jag.mod, c('mu', 'rho', 'sigma', 'x_rand'), n.iter = n.iter)
            
            res0 = summary(mcm.sam)
            res = as.data.frame(res0$statistics)
            daf = res['rho', ]
            
            rownames(daf) = gen.dat.original[gg]
            daf$Gene = gen.dat[gg]
            daf$Lineage = lineage
            daf$Gene.Transcript = gen.dat.original[gg]
            mcm.lin = rbind(mcm.lin, daf)
          }
        }, error = function(e) e, finally = {})
      }
    }
    # Rename for clarity: posterior mean of rho
    colnames(mcm.lin)[colnames(mcm.lin) == 'Mean'] <- 'rho'
    mcm.lin$z = mcm.lin$rho/mcm.lin$SD
    # mcm.lin$P.Value = exp(-0.717*mcm.lin$z-0.416*(mcm.lin$z)**2)
    # mcm.lin$P.Value = pchisq(mcm.lin$z**2, df=1, lower=F) 
    mcm.lin$P.Value = pnorm(abs(mcm.lin$z), mean = 0, sd = 1, lower.tail = F)  
    if (recalculate.FDR) {
      mcm.lin$adj.P.Val = p.adjust(mcm.lin$P.Value, method = 'BH')
    }
    mcm.lin = mcm.lin[order(-log10(mcm.lin$adj.P.Val+1e-323)*mcm.lin$rho, decreasing = F), ]
    # WRITE TABLE
    if (out) xlsx::write.xlsx(mcm.lin, opath, row.names = FALSE)
  }
  # mcm.lin.top = head(mcm.lin[mcm.lin$adj.P.Val < thr.FDR,], num.top.genes) #***
  # mcm.lin.all = rbind(mcm.lin.all, mcm.lin)
  # mcm.lin.all.top = rbind(mcm.lin.all.top, mcm.lin.top)
  .t1 <- Sys.time()
  message(sprintf('[BEACON] Finished %s — %s in %.1f min',
                data, lineage, as.numeric(difftime(.t1, .t0, units='mins'))))
  rm('mcm.lin')
}
.run_end <- Sys.time()
.run_log$end_time <- as.character(.run_end)
.run_log$elapsed_mins <- round(as.numeric(difftime(.run_end, .run_start, units="mins")), 1)

# Write run log + session info
writeLines(c(
  paste0("Start: ", .run_log$start_time),
  paste0("End:   ", .run_log$end_time),
  paste0("Elapsed (min): ", .run_log$elapsed_mins),
  "",
  "SessionInfo():"
), file.path(dir_path, "RUNLOG_lineage.txt"))
sink(file.path(dir_path, "RUNLOG_lineage.txt"), append = TRUE); print(sessionInfo()); sink()

