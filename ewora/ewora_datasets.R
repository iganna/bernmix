# source("https://bioconductor.org/biocLite.R")
# biocLite("KEGGdzPathwaysGEO")
# biocLite("PADOG")
# biocLite("annotate")
# biocLite("hgu133a.db")
# biocLite("hgu133plus2.db")
# biocLite("KEGGESET")


library("KEGGdzPathwaysGEO")
library("PADOG")
library(limma)
library("annotate")
library("hgu133a.db")
library("hgu133plus2.db")
library('KEGG.db')
require(stringr)




platforms = c(hgu133a.db=hgu133a.db, hgu133plus2.db=hgu133plus2.db)

mysets=data(package="KEGGdzPathwaysGEO")$results[,"Item"]
mysets

dataset = mysets[3]
# ===============================================

list_de_probes = list()
list_de_genes = list()
id_list = 1
for (dataset in mysets)
{
data(list=dataset,package="KEGGdzPathwaysGEO")
data = get(dataset)
experiment = experimentData(data) # exp
expression_mx = exprs(data) # esetm
disease_status = factor(pData(data)$Group) # G
block_status = factor(pData(data)$Block) # block
experim_design = notes(experiment)$design # Paired - Not Paired
annotation = paste(data@annotation,".db",sep="") # annotation of probes
target_kegg_id= notes(experiment)$targetGeneSets # Target KEGG pathway


if (experim_design == "Paired") {
  model_design <- model.matrix(~0 + disease_status + block_status)
  colnames(model_design) <- c(str_sub(colnames(model_design)[c(1,2)],start = -1), # c and d columns
                              colnames(model_design)[c(-1,-2)]) # factors of blocks
} else {
  model_design <- model.matrix(~0 + disease_status)
  colnames(model_design) <- levels(disease_status)
}



fit <- lmFit(expression_mx, model_design)
contrast_mx <- makeContrasts(contrasts = "d-c", levels = model_design)
fit2 <- contrasts.fit(fit, contrast_mx)
fit2 <- eBayes(fit2)
probe_signif <- topTable(fit2, 1, number=Inf, adjust.method='fdr')


# PADOG: get rid of duplicates probesets per ENTREZ ID by keeping the probeset
# with smallest p-value (computed using limma)
# EWORA: here we use the same assumption but for GENE SYMBOL
gene_annotation = select(platforms[[annotation]], rownames(probe_signif), c("SYMBOL"))
# remove rows with NA
gene_annotation = gene_annotation[complete.cases(gene_annotation), ]
gene_unique = unique(gene_annotation[,'SYMBOL'])
gene_signif = matrix(, nrow = 0, ncol = dim(probe_signif)[2])
for(gene in gene_unique)
{
  probes = unique(gene_annotation[gene_annotation[,'SYMBOL'] == gene, 'PROBEID'])
  row_to_add = probe_signif[probes,]
  # if these probes has non-zero intersection
  if (length(probes) != 1)
    row_to_add = row_to_add[row_to_add[,'P.Value'] == min(row_to_add[,'P.Value']),]
  rownames(row_to_add) <- gene
  gene_signif = rbind(gene_signif, row_to_add)
}

# Citation from Tarca et al. 2013 ; "A Comparison of Gene Set Analysis Methods 
# in Terms of Sensitivity, Prioritization and Specificity", page 9.
# The selection of DE genes for ORA was based on a moderated t-test p-value. 
# The following strategy was used for gene selection for ORA: 
# 1) use all genes with FDR [34] adjusted p-values<0.1 if more than 200; 
# else go to next option; 
# 2) use all genes with nominal p-values<0.05 and fold change.1.5 if more than 200; 
# else go to next option; 
# 3) Use top 1% of genes ranked by p-values.

# de_genes = rownames(gene_signif[gene_signif[,'adj.P.Val'] < 0.1,])

# FOR PROBES
de_probes = rownames(probe_signif[probe_signif[,'adj.P.Val'] < 0.01,])   # 1 option
if (length(de_probes) < 200) 
  de_probes = rownames(probe_signif[(probe_signif[,'P.Value'] < 0.05) & (abs(probe_signif[,'logFC']) > 1.5),])  # 2 option
if (length(de_probes) < 200)
  de_probes = rownames(probe_signif[1:200,])                            # 3 option

# FOR GENES
de_genes = rownames(gene_signif[gene_signif[,'adj.P.Val'] < 0.01,])   # 1 option
if (length(de_genes) < 200) 
  de_genes = rownames(gene_signif[(gene_signif[,'P.Value'] < 0.05) & (abs(gene_signif[,'logFC']) > 1.5),])  # 2 option
if (length(de_genes) < 200)
  de_genes = rownames(gene_signif[1:200,])                            # 3 option

list_de_probes[[id_list]]=probe_signif
list_de_genes[[id_list]]=gene_signif
id_list = id_list + 1

}

save(list_de_probes, list_de_genes, file = "/Users/anna/OneDrive/polytech/stud_projects/enrichment/Article_bio/de_lists1.RData")

n_de_genes = c()
de_genes_extracted = list()
id_list = 1
for(gene_signif in list_de_genes)
{

  de_genes = rownames(gene_signif[gene_signif[,'adj.P.Val'] < 0.1 & (abs(gene_signif[,'logFC']) > 1),])   # 1 option
  # de_genes = rownames(gene_signif[gene_signif[,'adj.P.Val'] < 0.1,])   # 1 option
  if (length(de_genes) < 200) 
    de_genes = rownames(gene_signif[(gene_signif[,'P.Value'] < 0.05) & (abs(gene_signif[,'logFC']) > 1.5),])  # 2 option
  if (length(de_genes) < 200)
    de_genes = rownames(gene_signif[1:200,])                            # 3 option
  
  de_genes = de_genes[1:min(3000, length(de_genes))]
  n_de_genes = c(n_de_genes, length(de_genes))
  de_genes_extracted[[id_list]] = de_genes
  id_list = id_list + 1
}

write(c(de_genes_extracted[[2]],
        de_genes_extracted[[5]]), file = "/Users/anna/OneDrive/polytech/stud_projects/enrichment/Article_bio/de.txt")





