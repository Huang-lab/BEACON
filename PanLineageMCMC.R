# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# Script: PanLineageMCMC.R
# Purpose: Pan-lineage BEACON (Bayesian correlation) to quantify GED/PED across all cell lines
# Inputs:
#   - DepMap sample info: DepMap_data/sample_info_22Q2.csv
#   - mRNA: CCLE_expression_22Q2.csv.gz  OR  Protein: mmc2.xlsx (Nusinow et al.)
#   - CRISPR dependency: DepMap_data/CRISPR_gene_effect_22Q2.csv.gz
# Outputs:
#   - Excel table with posterior summaries (per gene): Table.<DATA>.dependency.Bayesian.pancancer.xlsx
#   - Log file, reporting the session info, including runtime, and the details of the computation environment and package versions used. 
# Download links:
#   - sample_info.csv file @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=35020903
#   - CCLE_expression.csv file (further gzipped) @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34989919
#   - Supplementary data (Table S2: normalized protein expressions) from Nusinow et al. paper (doi.org/10.1016/j.cell.2019.12.023) @ https://www.cell.com/cms/10.1016/j.cell.2019.12.023/attachment/3709dedc-3a01-4e1d-ab4c-82597295c5d2 
#   - CRISPR_gene_effect.csv file (further gzipped) @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34990036
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
rm(list = ls(all.names = TRUE))
options(java.parameters = '-Xmx8000m')
library(openxlsx)
library(rjags)
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
set.seed(1234) # fixed seed for reproducibility
# Start run clock and capture session info later
.run_start <- Sys.time()
# Log capture: will be written to the output directory at the end
.run_log <- list(
  script = "PanLineageMCMC.R",
  start_time = as.character(.run_start),
  params = list()
)

# ---------------- Parameters (edit as needed) ----------------
n.adapt  = 200   # JAGS adaptation steps (short keeps runtime low; increase for stability)
n.update = 200   # burn-in updates before sampling (increase for stability)
n.iter   = 1000   # MCMC iterations per chain (increase -> tighter posteriors, longer runtime)

reproduce.results = TRUE   # if TRUE, uses current matrices as loaded
recalculate.FDR   = TRUE   # recompute BH FDR
out               = TRUE   # write Excel outputs
continue.from     = NULL   # resume (1..100 = percent chunks) or NULL
# continue.from = 7

data       = 'mRNA'        # 'mRNA' | 'Protein'
panel      = ''            # e.g. '.druggable'
cell.type  = 'All'         # 'All' | 'Primary' | 'Metastasis'
thr.FDR    = .05
# -------------------------------------------------------------

# Save parameters to run log
.run_log$params <- list(n.adapt=n.adapt, n.update=n.update, n.iter=n.iter,
                        data=data, panel=panel, cell.type=cell.type,
                        thr.FDR=thr.FDR)

# Output directory (relative)
dir.path <- file.path('..','out',
  paste0('jags.nadapt', n.adapt, '.update', n.update, '.mcmc', n.iter, '.simulation_SD_22Q2'))
if (!dir.exists(dir.path)) dir.create(dir.path, recursive = TRUE)

# Helper to build output filename
out_xlsx <- function(data, panel, continue.from=NULL) {
  file.path(dir.path, paste0('Table.', data, '.dependency.Bayesian.pancancer', panel, continue.from, '.xlsx'))
}

# ---------------- Inputs (relative) ----------------
depmap_info_path = file.path('..','..','..','..','Huang_lab_data','DepMap_data','sample_info_22Q2.csv') #| sample_info.csv file @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=35020903
mrna_path        = file.path('..','..','..','..','Huang_lab_data','DepMap_data','CCLE_expression_22Q2.csv.gz')  #| CCLE_expression.csv file (further gzipped) @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34989919
protein_xlsx     = file.path('..','..','..','..','Huang_lab_data','QuantProtCCLE_Nusinow_Cell2020','mmc2.xlsx') #| Supplementary data (Table S2: normalized protein expressions) from Nusinow et al. paper (doi.org/10.1016/j.cell.2019.12.023) @ https://www.cell.com/cms/10.1016/j.cell.2019.12.023/attachment/3709dedc-3a01-4e1d-ab4c-82597295c5d2 
crispr_path      = file.path('..','..','..','..','Huang_lab_data','DepMap_data','CRISPR_gene_effect_22Q2.csv.gz') #| CRISPR_gene_effect.csv file (further gzipped) @ https://figshare.com/articles/dataset/DepMap_22Q2_Public/19700056/2?file=34990036
# ---------------------------------------------------

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

# sam.dep = read.csv('~/Downloads/sample_info_20Q1.csv')
sam.dep = read.csv(depmap_info_path) 

table(sam.dep$primary_or_metastasis)
if (cell.type!='All'){
  sam.dep = sam.dep[sam.dep$primary_or_metastasis==cell.type,]
}

# Map DepMap IDs to CCLE names and extract tissue lineage from cell line labels
dep.dat.map = sam.dep[,c(1,3)]
dep.dat.map$Tissue_Type = as.character(sam.dep$CCLE_Name)
# sam.dep$Tissue_Type = as.character(sam.dep$CCLE_Name)
i=1
# Derive tissue lineage label from CCLE name (strip leading ACH code)
for (i in seq_len(nrow(dep.dat.map))){
  dep.dat.map$Tissue_Type[i] = sub('^[a-zA-Z0-9]*\\_','',dep.dat.map$Tissue_Type[i])
}

if (data == 'Protein') {
  # data exp
  exp.dat = openxlsx::read.xlsx(protein_xlsx, sheet='Normalized Protein Expression') 
  exp.dat = cbind(Gene_Symbol = exp.dat$Gene_Symbol, exp.dat[,-(1:50)])
  col.drp = which(colSums(regexprTab(expressions = c('^Protein', '^Description', '^Group', '^Uniprot', 'Peptides$', '^Column'), data = colnames(exp.dat))) > 0)
  exp.dat = exp.dat[,-col.drp]
  # Filter missing value genes
  # exp.dat = exp.dat[-which(rowMeans(is.na(exp.dat)) > .9),]
  tab.gen = table(exp.dat$Gene_Symbol)
  row.drp = c()
  for (gene in names(tab.gen)[which(tab.gen>1)]){
    idx = which(exp.dat$Gene_Symbol %in% gene) 
    row.drp = c(row.drp, setdiff(idx,idx[order(rowSums(is.na(exp.dat[idx,])))[1]]))
  }
  exp.dat = exp.dat[-row.drp,]
  # Filter missing value cells 
  # exp.dat = exp.dat[, -which(colMeans(is.na(exp.dat)) > .9)]
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
  # exp.dat = read.csv('~/Downloads/CCLE_expression_20Q1.csv.gz')
  exp.dat = read.csv(mrna_path) 
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

map = match(ccl.dep$DepMap_ID, dep.dat.map$DepMap_ID) #***
cel.dep = as.character((dep.dat.map$stripped_cell_line_name[map[!is.na(map)]]))
# cel.dep = as.character((dep.dat.map$cell_line_name[map[!is.na(map)]]))
# typ.dep = as.character(sam.dep$primary_or_metastasis[map[!is.na(map)]])
cel.tis.dep = as.character(unique( 
  apply(dep.dat.map[dep.dat.map$DepMap_ID %in% ccl.dep$DepMap_ID,-1],1,function(x){gsub('\\_','\\.',paste(x,collapse='.'))}) ))
  # apply(dep.dat.map[map[!is.na(map)],-1],1,function(x){gsub('\\_','\\.',paste(x,collapse='.'))}) )) #***
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

# Filter out lineages with <=5 cell lines (insufficient power for correlation model)
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
# ALIGN CELL LINES
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# Align expression and dependency matrices to the same cell line order
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
# BAYESIAN - PANCANCER
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
if (reproduce.results) {
  # the gene's expression is the dependent variable
  dep.var = apply(exp.dat, 2, function(x){as.numeric(x)})
  rownames(dep.var) = gen.dat
  # the gene's dependency is the independent variable
  indep.var = apply(ccl.dep[,-1], 1, function(x){as.numeric(x)})
  rownames(indep.var) = gen.dep
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

res.all = c()
genes.query = intersect(gen.dat, gen.dep); i = 0; L = length(genes.query) #**
if (!is.null(continue.from)){
  if (continue.from > 99 | continue.from < 0){
    print('Warning: please set "continue.from" to an integer between 0 and 99.')
  } else if (continue.from > 0){
    i = ceiling(seq(from = L/100, to = L, by = L/100))[continue.from]
    genes.query = genes.query[-(1:i)]; 
  }
}
gene = 'SOX10'
for (gene in genes.query) { #*** [1:100]
  i = i + 1 #**

  r = tryCatch({
    x = cbind(dep.var[which(gen.dat %in% gene),], indep.var[which(gen.dep %in% gene),])
    x = x[rowSums(is.na(x)) == 0,]
    rho.pri = cor(x)[2]
    
    dat.lis = list(x = x, n = nrow(x))
    ini.lis = list(mu = colMeans(x, na.rm=T), rho = rho.pri, sigma = apply(x, 2, function(x){sd(x, na.rm = T)}))

    # Create initialization list with RNG specifications for each chain
    ini.lis_with_seeds <- list(
      c(ini.lis, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1)),
      c(ini.lis, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2)),
      c(ini.lis, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 3))
    )
    
    jag.mod = jags.model(textConnection(model_string), data = dat.lis, inits = ini.lis_with_seeds, n.adapt = n.adapt, n.chains = 3, quiet = T)
    update(jag.mod, n.update, progress.bar = 'none')
    mcm.sam = coda.samples(jag.mod, c('mu', 'rho', 'sigma', 'x_rand'), n.iter = n.iter)
    
    res0 = summary(mcm.sam)
    res = as.data.frame(res0$statistics)
    res.out = res['rho', ]
    
    rownames(res.out) = gene
    res.out$Gene = gene
    colnames(res.out)[colnames(res.out) == 'Mean'] <- 'rho'  # rename for clarity

    # res.out$z = res.out$rho/res.out$`Time-series SE`
    res.out$z = res.out$rho/res.out$SD
    res.out$P.Value = pnorm(abs(res.out$z), mean = 0, sd = 1, lower.tail = F)  
    res.out
    
    
    
    res.all = rbind(res.all, res.out)
    # res.out$SD > 20*res.out$`Time-series SE` # SE is trivial 
    #' MCMC samples auto-correlate
    #' 'Time-series SE' considers auto-correlation
    #' Reduce it by using larger 'n.iter'
  }, error = function(e) e, finally = {})

  per.don = which(ceiling(seq(from = L/100, to = L, by = L/100)) %in% i) #**
  if (length(per.don) > 0) {
    message(per.don, '% done.')
    # WRITE TABLE
    if (out){
      xlsx::write.xlsx(res.all, out_xlsx(data, panel, continue.from))
    }
  } #**
  
}
# WRITE TABLE
if (out){
  xlsx::write.xlsx(res.all, out_xlsx(data, panel, continue.from))

}
                                                          
res.all$z = res.all$rho/res.all$SD
# res.all$P.Value = exp(-0.717*res.all$z-0.416*(res.all$z)**2)
# res.all$P.Value = pchisq(res.all$z**2, df=1, lower=F)  
res.all$P.Value = pnorm(abs(res.all$z), mean = 0, sd = 1, lower.tail = F)  
if (recalculate.FDR) {
  res.all$adj.P.Val = p.adjust(res.all$P.Value, method = 'BH')
}
res.all = res.all[order(-log10(res.all$adj.P.Val+1e-323)*res.all$rho, decreasing = F), ]

# WRITE TABLE
if (out){
  xlsx::write.xlsx(res.all, out_xlsx(data, panel, continue.from))
}
                                                          
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
.run_end <- Sys.time()
.run_log$end_time <- as.character(.run_end)
.run_log$elapsed_mins <- round(as.numeric(difftime(.run_end, .run_start, units="mins")), 1)
# Write session info + run log
if (exists("dir.path")) {
  writeLines(c(
    paste0("Start: ", .run_log$start_time),
    paste0("End:   ", .run_log$end_time),
    paste0("Elapsed (min): ", .run_log$elapsed_mins),
    "",
    "SessionInfo():"
  ), file.path(dir.path, "RUNLOG_pancancer.txt"))
  sink(file.path(dir.path, "RUNLOG_pancancer.txt"), append = TRUE); print(sessionInfo()); sink()
}
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

