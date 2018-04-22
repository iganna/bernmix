library("KEGG.db")
library("KEGGgraph")
library('KEGGREST')
library(stringr)


# This function translates the list of hsa annotated genes to official gene symbols
hsa2symbols <- function(hsa)
{
  g_symbols = c()
  for (i in seq(1, length(hsa), 10))
  {
    genes_info <- keggGet(hsa[i:min(i+9, length(hsa))])
    for (i_gene in 1:length(genes_info))
    {
      s = genes_info[[i_gene]]$NAME
      if (is.null(s))
      {
        g_symbols = c(g_symbols, NA)
        next
      }
      gene = strsplit(s, ', ')[[1]]
      # write gene[1]
      g_symbols = c(g_symbols, gene[1])
    }
  }
  return(g_symbols)
}


# Get the list of pathways to analyse
# organism = 'hsa'
# pathways_all = as.list(KEGGPATHID2EXTID)
# pathways_org = pathways_all[grep(organism, names(pathways_all))]
# pathways_ids = sub(paste("^", organism, sep = ""), "", 
#                    names(pathways_org))
# pathways_name = unlist(as.list(KEGGPATHID2NAME)[pathways_ids])


# Get the list of pathways to analyse by 
pathways_ids = unique(keggLink("pathway", "hsa"))
pathways_ids = sapply(pathways_ids, function(s) str_sub(s, start = -5))
pathways_name = c()


# Read each pathway and write in new format
path_pref = '/Users/anna/OneDrive/polytech/stud_projects/enrichment/scripts/bernmix/'
file_kegg_formatted = paste0(c(path_pref, "data/kegg_formatted.txt"), "", collapse="")
file.create(file_kegg_formatted)
pathway_cnt = 0
n_genes = c()
for (i in 120:length(pathways_ids))
{

  file_kgml = paste0(c(path_pref, "data/kgml/hsa", pathways_ids[i],  ".xml"), "", collapse="")
  if (!file.exists(file_kgml))
    next
  
  # Get pathway structure
  path_struct <- parseKGML2Graph(file_kgml, expandGenes=FALSE)
  path_info <- getPathwayInfo(parseKGML(file_kgml))
  pathways_name[i] <- slot(path_info, 'title')
  

  # Get Nodes and Edges
  nodes_info <- getKEGGnodeData(path_struct)
  edges_info <- getKEGGedgeData(path_struct)

  n_genes[i] <- length(nodes_info)
  if (length(nodes_info) > 400)  # ANNA
    next
  pathway_cnt = pathway_cnt + 1

  # Save PATHWAY
  # write pathway name
  write(sprintf('PATHWAY\t"%s"', pathways_name[i]), file=file_kegg_formatted, append=TRUE)

  # Get genes from nodes
  node_genes = list()
  for (i_node in 1:length(nodes_info))
  {
    node = nodes_info[[i_node]]
    enteryID = slot(node, 'entryID')
    genes = slot(node, 'name')
    node_genes[[enteryID]] <- genes
  }

  # Get new annotation for genes
  all_hsa_genes = c()
  for (g in node_genes)
    all_hsa_genes = c(all_hsa_genes, g)
  all_hsa_genes = unique(all_hsa_genes)
  all_symbol_genes = setNames(hsa2symbols(all_hsa_genes), all_hsa_genes)

  for (i_gene in 1:length(node_genes))
  {
    node_genes[[i_gene]] = all_symbol_genes[node_genes[[i_gene]]]
    node_genes[[i_gene]] <- node_genes[[i_gene]][!is.na(node_genes[[i_gene]])]
  }
  all_symbol_genes <- unique(all_symbol_genes[!is.na(all_symbol_genes)])

  # Get adjacency matrix og graph
  adj_mx = matrix(FALSE, length(node_genes), length(node_genes),
           dimnames= list(names(node_genes), names(node_genes)))
  if (length(edges_info) > 0)
    for (i_edge in 1:length(edges_info))
    {
      edge = edges_info[[i_edge]]
      # white edge id
      entry1 = slot(edge, 'entry1ID')
      entry2 = slot(edge, 'entry2ID')
      adj_mx[entry1, entry2] <- adj_mx[entry2, entry1] <- TRUE
    }

  # Collapse nodes with common genes
  for (i_node1 in names(node_genes))
    for (i_node2 in names(node_genes))
    {
      if (i_node1 == i_node2)
        next

      common_genes = intersect(node_genes[[i_node1]], node_genes[[i_node2]])
      if (length(common_genes) > 0)
      {
        if (length(node_genes[[i_node1]]) >= length(node_genes[[i_node2]]))
        {
          i_big <- i_node1
          i_small <- i_node2
        }
        else
        {
          i_big <- i_node2
          i_small <- i_node1
        }
        adj_mx[i_big,] = adj_mx[i_big,] | adj_mx[i_small,]
        adj_mx[,i_big] = adj_mx[,i_big] | adj_mx[,i_small]
        node_genes[[i_small]] = setdiff(node_genes[[i_small]], common_genes)
      }
    }

  # Remove empty nodes
  node_to_remain = setNames(rep(TRUE,  length(node_genes)), names(node_genes))
  for(i_gene in names(node_genes))
  {
    if (length(node_genes[[i_gene]]) == 0)
      node_to_remain[i_gene] = FALSE
  }
  node_genes = node_genes[node_to_remain]
  if (sum(node_to_remain) == 1)
    adj_mx = matrix(adj_mx[node_to_remain,node_to_remain],
                     dimnames = list(names(node_genes), names(node_genes)))
  else
    adj_mx = adj_mx[node_to_remain,node_to_remain]

  # Save GENES
  write(paste0(c('GENES', length(all_symbol_genes), all_symbol_genes), sep = '\t', collapse = ''),
        file=file_kegg_formatted, append=TRUE)

  # Save NODES
  # write numer of nodes
  write(sprintf('NODES\t%i', length(node_genes)), file=file_kegg_formatted, append=TRUE)
  # write node id
  for (i_node in names(node_genes))
    write(paste0(c(i_node, sum(adj_mx[i_node,]), node_genes[[i_node]]), sep='\t', collapse = ''),
          file=file_kegg_formatted, append=TRUE)

}



