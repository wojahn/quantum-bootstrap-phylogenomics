# =============================================================================
# Quantum Bootstrap Phylogenomics Pipeline
# J. M. A. Wojahn, P. E. Ellestad, and S. Buerki
# Reproducible phylogenomic inference with quantum-inspired bootstrap analysis
# DOI: 10.5281/zenodo.20534276
# =============================================================================
 
 
# =============================================================================
# SECTION 1: SOFTWARE INSTALLATION
# =============================================================================
 
# 1.1 Install CRAN-hosted core phylogenetics and infrastructure tools
cran_packages <- c("ape", "phangorn", "treeio", "TreeDist", "reticulate", "remotes")
new_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
 
if (length(new_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(new_cran, collapse = ", "))
  install.packages(new_cran, repos = "https://cloud.r-project.org")
}
 
# 1.2 Install Bioconductor packages required for sequence alignment
bioc_packages <- c("DECIPHER", "Biostrings")
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
 
if (length(new_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(new_bioc, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
}
 
# 1.3 Install GitHub-hosted custom tools
if (!requireNamespace("PrideBar", quietly = TRUE)) {
  message("Installing Custom Development Package: PrideBar from GitHub")
  remotes::install_github("wojahn/PrideBar", upgrade = "never")
}
 
# External software (run once in terminal, not from R):
#   brew install iqtree2
#   git clone https://github.com/Nicofero/PhyloTree_ReConst.git
#   cd PhyloTree_ReConst
#   cp qaoa_functions.py qa_phylo_tree.ipynb qa_hybrid_trees.py qa_functions.py qa_big_trees.py ..
#   pip install numpy pandas networkx pyqubo dwave-ocean-sdk
#   dwave setup --auth
#   conda install -c bioconda seaview
 
# Load libraries
library(ape)
library(phangorn)
library(TreeDist)
library(reticulate)
library(DECIPHER)
library(Biostrings)
library(PrideBar)
 
 
# =============================================================================
# SECTION 2: REPRODUCIBILITY
# =============================================================================
 
set.seed(35753)
 
 
# =============================================================================
# SECTION 3: TAXON SELECTION AND SEQUENCE ALIGNMENT
# =============================================================================
 
# Load master plastome file from Wojahn et al. (2023)
plastomes <- ape::read.FASTA("AlignedOrderedGoodPlastomes_First.fasta")
 
# Extract the strategic 19-taxon cohort (16 ingroup, 2 Benstonea, 1 Freycinetia)
plastomes <- plastomes[c(6, 8, 9, 44, 30, 51, 59, 13, 119, 63, 84, 17, 106, 68, 35, 114, 113, 43, 31)]
 
names(plastomes) <- c(
  "Freycinetia_baueriana", "Benstonea_copelandii",
  "Benstonea_herbacea", "Pandanus_furcatus", "Pandanus_calcis",
  "Pandanus_irregularis", "Pandanus_lacuum", "Pandanus_amaryllifolius",
  "Pandanus_vitiensis", "Pandanus_luzonensis", "Pandanus_polyglossus",
  "Pandanus_aquaticus", "Pandanus_tectorius", "Pandanus_maximus",
  "Pandanus_concretus", "Pandanus_utilis", "Pandanus_tsingycola",
  "Pandanus_eydouxia", "Pandanus_callmanderianus"
)
 
ape::write.FASTA(plastomes, "plastomes.fasta")
 
# Re-align using DECIPHER profile-to-profile approach
plastomez <- readDNAStringSet("plastomes.fasta")
plastomez <- DECIPHER::RemoveGaps(plastomez)
plastomez_a <- DECIPHER::AlignSeqs(plastomez, processors = 10)
plastomez_a <- DECIPHER::AdjustAlignment(plastomez_a)
Biostrings::writeXStringSet(plastomez_a, "Pandanus_a.fasta")
 
 
# =============================================================================
# SECTION 4: CLASSICAL MAXIMUM LIKELIHOOD REFERENCE TREE
# =============================================================================
 
iqtree_cmd <- "iqtree -s Pandanus_a.fasta -m MFP -bb 1000 -nt AUTO -redo -o Freycinetia_baueriana -seed 35753 -safe"
system(iqtree_cmd)
 
# Parse and root the ML reference tree
classical_ref <- ape::read.tree("Pandanus_a.fasta.treefile")
classical_ref <- ape::root(classical_ref, "Freycinetia_baueriana", resolve.root = TRUE)
classical_ref <- ape::drop.tip(classical_ref, "Freycinetia_baueriana")
classical_ref <- ape::ladderize(classical_ref)
 
pdf("IQTREE_tree.pdf", width = 15, height = 20)
plot(classical_ref, show.node.label = FALSE, cex = 1.0, root.edge = TRUE)
nodelabels(classical_ref$node.label, adj = c(1.2, -1.5),
           frame = "none", cex = 0.8, font = 2, col = "darkblue")
dev.off()
 
 
# =============================================================================
# SECTION 5: CLASSICAL FREQUENTIST BOOTSTRAP RESAMPLING
# =============================================================================
 
# Read alignment in phyDat format
final_alignment <- phangorn::read.phyDat("Pandanus_a.fasta", format = "fasta")
 
# Generate 100 non-parametric bootstrap replicates
boot_alignments <- phangorn::bootstrap.phyDat(final_alignment, FUN = identity, bs = 100)
dir.create("bootstrap_inputs", showWarnings = FALSE)
 
# Read in raw reference for downstream use
raw_ref <- ape::read.tree("Pandanus_a.fasta.treefile")
raw_ref$tip.label <- chartr("-", "_", raw_ref$tip.label)
raw_ref$tip.label <- trimws(raw_ref$tip.label)
raw_ref_unrooted  <- ape::unroot(raw_ref)
max_possible_splits <- length(raw_ref_unrooted$tip.label) - 3
 
# Build and export cophenetic distance matrices for each replicate
PrideBar::SetPrideBar(1, 100, 2)
for (i in 1:100) {
  PrideBar::PrideBar()
  rep_data <- boot_alignments[[i]]
 
  tmp_dm   <- phangorn::dist.ml(rep_data, model = "JC69")
  tmp_tree <- ape::nj(tmp_dm)
  tmp_fit  <- phangorn::pml(tmp_tree, data = rep_data)
 
  suppressMessages(
    tmp_opt <- phangorn::optim.pml(
      tmp_fit, model = "GTR", k = 4,
      optEdge = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE,
      rearrangement = "NNI", control = phangorn::pml.control(trace = 0)
    )
  )
 
  rep_matrix <- as.matrix(cophenetic(tmp_opt$tree))
  file_name  <- paste0("bootstrap_inputs/matrix_", sprintf("%03d", i), ".csv")
  write.csv(rep_matrix, file = file_name, row.names = TRUE)
}
PrideBar::FadePrideBar()
 
 
# =============================================================================
# SECTION 6: QUANTUM / QUANTUM-INSPIRED BOOTSTRAP TREE RECONSTRUCTION
# =============================================================================
 
# --- 6.1 Physical Quantum Annealing (D-Wave QPU) ---
# Use this block if you have a D-Wave Leap API key and QPU access.
 
dwave_qpu_script <- c(
  "import os, sys, glob",
  "import pandas as pd",
  "import numpy as np",
  "import qa_functions as qaf",
  "from dwave.system import DWaveSampler, EmbeddingComposite",
  "",
  "matrix_dir   = 'bootstrap_inputs' if os.path.exists('bootstrap_inputs') else '../bootstrap_inputs'",
  "matrix_files = sorted(glob.glob(f'{matrix_dir}/matrix_*.csv'))",
  "",
  "if not matrix_files:",
  "    print(f\"Error: Cannot find matrix files in '{matrix_dir}'. Execution halted.\")",
  "    sys.exit(1)",
  "",
  "print('Connecting to D-Wave Live QPU Hardware via Leap...')",
  "try:",
  "    qpu_sampler = DWaveSampler(solver={'topology__type': 'pegasus'})",
  "    qpu_backend = EmbeddingComposite(qpu_sampler)",
  "    print(f\" Connected successfully to: {qpu_sampler.properties['id']}\")",
  "except Exception as e:",
  "    print(f'Hardware Connection Error: {e}')",
  "    sys.exit(1)",
  "",
  "bootstrap_topologies = []",
  "for i, file_path in enumerate(matrix_files, 1):",
  "    df     = pd.read_csv(file_path, index_col=0)",
  "    matrix = df.to_numpy()",
  "    taxa_labels = list(df.index)",
  "    max_dist = np.max(matrix) if np.max(matrix) > 0 else 1.0",
  "    similarity_matrix = max_dist - matrix",
  "    np.fill_diagonal(similarity_matrix, 0)",
  "    int_tags = list(range(len(taxa_labels)))",
  "    native_tree_object = qaf.sa_phylo_tree(",
  "        similarity_matrix, tags=int_tags, sampler=qpu_backend,",
  "        num_reads=1000, annealing_time=20)",
  "    newick_str = native_tree_object.to_newick(labels=taxa_labels)",
  "    bootstrap_topologies.append(newick_str)",
  "    if i % 5 == 0:",
  "        print(f'  -> Processed {i}/100 profiles on the QPU.')",
  "",
  "with open('quantum_bootstrap_trees.tre', 'w') as f:",
  "    for tree in bootstrap_topologies:",
  "        f.write(f'{tree}\\n')",
  "print('100 QPU-annealed topologies saved to quantum_bootstrap_trees.tre.')"
)
 
writeLines(dwave_qpu_script, "run_quantum_bootstrap_QPU.py")
 
# --- 6.2 Simulated Quantum Annealing (Classical Simulator) ---
# Use this block on standard hardware without QPU access.
 
dwave_sim_script <- c(
  "import os, sys, glob",
  "import pandas as pd",
  "import numpy as np",
  "import qa_functions as qaf",
  "import neal",
  "",
  "matrix_dir   = 'bootstrap_inputs' if os.path.exists('bootstrap_inputs') else '../bootstrap_inputs'",
  "matrix_files = sorted(glob.glob(f'{matrix_dir}/matrix_*.csv'))",
  "",
  "if not matrix_files:",
  "    print(f\"Error: Cannot find matrix files in '{matrix_dir}'. Execution halted.\")",
  "    sys.exit(1)",
  "",
  "print('Initializing classical Simulated Annealing simulator backend via neal...')",
  "try:",
  "    sim_backend = neal.SimulatedAnnealingSampler()",
  "    print(' Classical simulation environment established successfully.')",
  "except Exception as e:",
  "    print(f'Backend Initialization Error: {e}')",
  "    sys.exit(1)",
  "",
  "bootstrap_topologies = []",
  "for i, file_path in enumerate(matrix_files, 1):",
  "    df     = pd.read_csv(file_path, index_col=0)",
  "    matrix = df.to_numpy()",
  "    taxa_labels = list(df.index)",
  "    max_dist = np.max(matrix) if np.max(matrix) > 0 else 1.0",
  "    similarity_matrix = max_dist - matrix",
  "    np.fill_diagonal(similarity_matrix, 0)",
  "    int_tags = list(range(len(taxa_labels)))",
  "    native_tree_object = qaf.sa_phylo_tree(",
  "        similarity_matrix, tags=int_tags, sampler=sim_backend, num_reads=1000)",
  "    newick_str = native_tree_object.to_newick(labels=taxa_labels)",
  "    bootstrap_topologies.append(newick_str)",
  "    if i % 5 == 0:",
  "        print(f'  -> Processed {i}/100 profiles using local simulated annealing.')",
  "",
  "with open('quantum_inspired_bootstrap_trees.tre', 'w') as f:",
  "    for tree in bootstrap_topologies:",
  "        f.write(f'{tree}\\n')",
  "print('100 simulated annealing topologies saved to quantum_inspired_bootstrap_trees.tre.')"
)
 
writeLines(dwave_sim_script, "run_quantum_bootstrap.py")
 
# Run the simulator via reticulate (comment out to run QPU version instead)
reticulate::use_condaenv("base", required = TRUE)
reticulate::py_run_file("run_quantum_bootstrap.py")
 
 
# =============================================================================
# SECTION 7: BOOTSTRAP TREE INGESTION AND RESOLUTION AUDIT
# =============================================================================
 
target_tree_file  <- "Pandanus_a.fasta.treefile"
quantum_boot_file <- "quantum_bootstrap_trees.tre"
 
raw_ref <- ape::read.tree(target_tree_file)
raw_ref$tip.label <- chartr("-", "_", raw_ref$tip.label)
raw_ref$tip.label <- trimws(raw_ref$tip.label)
ref_labels        <- sort(raw_ref$tip.label)
raw_ref_unrooted  <- ape::unroot(raw_ref)
max_possible_splits <- length(raw_ref_unrooted$tip.label) - 3
 
boot_lines <- readLines(quantum_boot_file)
boot_lines <- boot_lines[boot_lines != ""]
 
clean_boot_list <- lapply(seq_along(boot_lines), function(i) {
  tr <- tryCatch(ape::read.tree(text = boot_lines[i]), error = function(e) NULL)
  if (!is.null(tr)) {
    tr$tip.label <- chartr("-", "_", tr$tip.label)
    tr$tip.label <- trimws(tr$tip.label)
    rogue_taxa   <- tr$tip.label[!tr$tip.label %in% ref_labels]
    if (length(rogue_taxa) > 0) tr <- ape::drop.tip(tr, rogue_taxa)
  }
  return(tr)
})
clean_boot_list <- clean_boot_list[!sapply(clean_boot_list, is.null)]
class(clean_boot_list) <- "multiPhylo"
 
resolution_percentages <- sapply(seq_along(clean_boot_list), function(i) {
  boot_tree    <- ape::unroot(clean_boot_list[[i]])
  actual_splits <- min(boot_tree$Nnode, max_possible_splits)
  return((actual_splits / max_possible_splits) * 100)
})
 
message("\n=== QUANTUM BOOTSTRAP RESOLUTION AUDIT ===")
print(summary(resolution_percentages))
 
 
# =============================================================================
# SECTION 8: QUANTUM BOOTSTRAP SUPPORT MAPPING
# =============================================================================
 
supported_tree <- phangorn::plotBS(raw_ref, clean_boot_list, type = "none", p = 0)
supported_tree <- ape::root(supported_tree, "Freycinetia_baueriana", resolve.root = TRUE)
supported_tree <- ape::drop.tip(supported_tree, "Freycinetia_baueriana")
supported_tree <- ape::ladderize(supported_tree)
qbs_values     <- supported_tree$node.label
 
pdf("quantum_bootstrap_hybrid_tree.pdf", width = 12, height = 15)
par(mar = c(5, 8, 4, 2))
ape::plot.phylo(supported_tree, cex = 1.0, font = 3, edge.width = 2, label.offset = 0.0001)
ape::nodelabels(
  text = paste0(as.character(as.numeric(qbs_values) * 100),
                ifelse(qbs_values == "", "", "%")),
  adj = c(1, -0.5), frame = "none", cex = 0.8, font = 2, col = "darkblue"
)
dev.off()
 
 
# =============================================================================
# SECTION 9: DENSITREE CLOUDOGRAM VISUALIZATION
# =============================================================================
 
rooted_boot_list <- lapply(seq_along(clean_boot_list), function(i) {
  tr            <- clean_boot_list[[i]]
  tr_rooted     <- ape::root(tr, outgroup = "Freycinetia_baueriana", resolve.root = TRUE)
  tr_bifurcated <- ape::multi2di(tr_rooted)
  tr_bifurcated$edge.length <- rep(1, nrow(tr_bifurcated$edge))
  ape::ladderize(ape::read.tree(text = ape::write.tree(tr_bifurcated)))
})
class(rooted_boot_list) <- "multiPhylo"
 
classical_ref_rooted <- phangorn::plotBS(raw_ref, clean_boot_list, type = "none", p = 0)
classical_ref_rooted <- ape::root(classical_ref_rooted, outgroup = "Freycinetia_baueriana", resolve.root = TRUE)
classical_ref_rooted <- ape::multi2di(classical_ref_rooted)
classical_ref_rooted$edge.length <- rep(1, nrow(classical_ref_rooted$edge))
 
pdf("Quantum_Bootstrap_Densitree.pdf", width = 12, height = 8)
par(mar = c(3, 3, 5, 10))
phangorn::densiTree(
  rooted_boot_list,
  type      = "cladogram",
  alpha     = 0.15,
  col       = "royalblue",
  consensus = classical_ref_rooted,
  direction = "rightwards",
  scaleX    = TRUE,
  cex       = 0.9,
  font      = 3,
  main      = "Quantum Bootstrap Densitree\n(Topological Variation Cloud vs. Classical Baseline)"
)
dev.off()
 
 
# =============================================================================
# SECTION 10: ROBINSON-FOULDS DISTANCE PROFILE
# =============================================================================
 
rf_distances <- sapply(seq_along(clean_boot_list), function(i) {
  TreeDist::RobinsonFoulds(supported_tree, clean_boot_list[[i]])
})
 
message("\n=== ROBINSON-FOULDS DISTANCE AUDIT ===")
print(summary(rf_distances))
 
pdf("Quantum_Bootstrap_RF_Histogram.pdf", width = 8, height = 6)
max_rf      <- max(rf_distances, na.rm = TRUE)
if (max_rf == 0) max_rf <- 2
hist_breaks <- seq(-0.5, max_rf + 0.5, by = 1)
 
hist(rf_distances, breaks = hist_breaks, col = "skyblue", border = "white",
     main = "Topological Divergence Profile\n(Quantum Bootstrap Cloud vs. Classical Baseline)",
     xlab = "Absolute Robinson-Foulds (RF) Distance",
     ylab = "Frequency (Number of Replicates)",
     xaxt = "n", las = 1)
axis(1, at = 0:max_rf, labels = 0:max_rf)
abline(v = median(rf_distances), col = "darkred", lwd = 2, lty = 2)
legend("topright",
       legend   = paste("Median RF Distance =", median(rf_distances)),
       col      = "darkred", lwd = 2, lty = 2, box.col = "white")
dev.off()
 
 
# =============================================================================
# SECTION 11: SUMMARY METRICS AND NODE SUPPORT TABLES
# =============================================================================
 
res_summary     <- summary(resolution_percentages)
rf_summary      <- summary(rf_distances)
raw_qbs_numeric <- as.numeric(qbs_values)
clean_qbs_pct   <- raw_qbs_numeric[!is.na(raw_qbs_numeric)] * 100
if (length(clean_qbs_pct) == 0) clean_qbs_pct <- 0
 
dataset_summary_matrix <- data.frame(
  Metric_Category = c(
    rep("Topological Resolution (%)", 3),
    rep("Robinson-Foulds (RF) Distance", 3),
    rep("Quantum Bootstrap Support (QBS)", 3)
  ),
  Statistical_Parameter = rep(c("Minimum", "Median", "Maximum"), 3),
  Value = c(
    round(res_summary[["Min."]], 2),
    round(res_summary[["Median"]], 2),
    round(res_summary[["Max."]], 2),
    round(rf_summary[["Min."]], 2),
    round(rf_summary[["Median"]], 2),
    round(rf_summary[["Max."]], 2),
    paste0(round(min(clean_qbs_pct),    1), "%"),
    paste0(round(median(clean_qbs_pct), 1), "%"),
    paste0(round(max(clean_qbs_pct),    1), "%")
  ),
  stringsAsFactors = FALSE
)
 
n_tips            <- length(supported_tree$tip.label)
internal_node_ids <- (n_tips + 1):(n_tips + supported_tree$Nnode)
 
formatted_pct_column <- sapply(raw_qbs_numeric, function(val) {
  if (is.na(val)) return("Root/Unresolved Split")
  paste0(round(val * 100, 1), "%")
})
 
node_support_table <- data.frame(
  Plot_Node_ID           = internal_node_ids,
  Raw_Support_Proportion = supported_tree$node.label,
  Formatted_QBS_Percentage = formatted_pct_column,
  stringsAsFactors = FALSE
)
 
# QBS vs UFBoot delta
uwu <- c(as.numeric(qbs_values) * 100)[-1] - as.numeric(raw_ref$node.label)[-1]
message("QBS vs UFBoot median delta: ", paste0(round(median(uwu), 2), " ± ", round(IQR(uwu), 2)))
 
print(dataset_summary_matrix, row.names = FALSE)
write.csv(dataset_summary_matrix, file = "Quantum_Bootstrap_Dataset_Summary.csv",    row.names = FALSE)
write.csv(node_support_table,     file = "Quantum_Bootstrap_Node_Support_Details.csv", row.names = FALSE)
 
 
# =============================================================================
# SECTION 12: NODAL DEPTH ANALYSIS OF LOW-SUPPORT CONFLICT
# =============================================================================
 
all_node_depths  <- ape::node.depth.edgelength(supported_tree)
internal_depths  <- all_node_depths[(n_tips + 1):length(all_node_depths)]
norm_depths      <- (internal_depths - min(internal_depths)) /
                    (max(internal_depths) - min(internal_depths))
 
support_depth_df <- data.frame(
  Node_ID       = (n_tips + 1):length(all_node_depths),
  QBS           = as.numeric(supported_tree$node.label) * 100,
  Relative_Depth = norm_depths
)
 
weak_nodes <- support_depth_df[!is.na(support_depth_df$QBS) & support_depth_df$QBS < 70, ]
 
if (nrow(weak_nodes) > 0) {
  median_weak_depth <- median(weak_nodes$Relative_Depth)
  message("\n=== LOW-SUPPORT MIDDLENESS AUDIT ===")
  message("Number of low-support nodes (<70%): ", nrow(weak_nodes))
  message("Median Relative Depth of weak nodes: ", round(median_weak_depth, 4))
 
  if (median_weak_depth > 0.4 && median_weak_depth < 0.6) {
    message("Interpretation: Weak support clustered in the MID-SECTION of the tree.")
  } else if (median_weak_depth <= 0.4) {
    message("Interpretation: Weak support clustered deep near the ROOT base.")
  } else {
    message("Interpretation: Weak support clustered out near the peripheral TIPS.")
  }
} else {
  message("All nodes possess high support; no topological weak zones detected.")
}
 
 
# =============================================================================
# SECTION 13: CLASSICAL DISTANCE BOOTSTRAP REFERENCE (CDBS)
# =============================================================================
 
matrix_files <- sort(list.files("bootstrap_inputs", pattern = "matrix_.*\\.csv", full.names = TRUE))
 
classical_dist_boot_list <- lapply(seq_along(matrix_files), function(i) {
  mat <- as.dist(read.csv(matrix_files[[i]], row.names = 1))
  ape::nj(mat)
})
class(classical_dist_boot_list) <- "multiPhylo"
 
cdbs_tree <- phangorn::plotBS(raw_ref, classical_dist_boot_list, type = "none", p = 0)
cdbs_tree <- ape::root(cdbs_tree, "Freycinetia_baueriana", resolve.root = TRUE)
cdbs_tree <- ape::drop.tip(cdbs_tree, "Freycinetia_baueriana")
cdbs_tree <- ape::ladderize(cdbs_tree)
 
rf_distances_cdbs <- sapply(seq_along(classical_dist_boot_list), function(i) {
  TreeDist::RobinsonFoulds(raw_ref_unrooted, classical_dist_boot_list[[i]])
})
 
message("\n=== CLASSICAL DISTANCE BOOTSTRAP RESOLUTION AUDIT ===")
print(summary(sapply(classical_dist_boot_list, function(tr) {
  ape::unroot(tr)$Nnode / max_possible_splits * 100
})))
 
message("\n=== CLASSICAL DISTANCE BOOTSTRAP RF DISTANCES ===")
print(summary(rf_distances_cdbs))
 
message("\n=== CLASSICAL DISTANCE BOOTSTRAP SUPPORT (CDBS) ===")
cdbs_raw <- as.numeric(cdbs_tree$node.label) * 100
print(summary(cdbs_raw[!is.na(cdbs_raw)]))
 
# Node-by-node QBS vs CDBS delta
qbs_numeric  <- as.numeric(qbs_values) * 100
cdbs_numeric <- cdbs_raw
min_nodes    <- min(length(qbs_numeric), length(cdbs_numeric))
qbs_numeric  <- qbs_numeric[1:min_nodes]
cdbs_numeric <- cdbs_numeric[1:min_nodes]
delta_support <- qbs_numeric - cdbs_numeric
 
message("\n=== QBS vs CDBS NODE-BY-NODE DELTA ===")
print(summary(delta_support[!is.na(delta_support)]))
message("Median delta: ", paste0(round(median(delta_support, na.rm = TRUE), 2),
        " ± ", round(IQR(delta_support, na.rm = TRUE), 2)))
 
# Side-by-side QBS vs CDBS comparison tree plot
pdf("QBS_vs_CDBS_Comparison.pdf", width = 14, height = 15)
par(mfrow = c(1, 2), mar = c(5, 8, 4, 2))
 
ape::plot.phylo(supported_tree, cex = 1.0, font = 3, edge.width = 2,
                label.offset = 0.0001, main = "Quantum Bootstrap Support (QBS)")
ape::nodelabels(
  text = paste0(as.character(round(qbs_numeric, 1)), ifelse(is.na(qbs_numeric), "", "%")),
  adj = c(1, -0.5), frame = "none", cex = 0.8, font = 2, col = "darkblue"
)
 
ape::plot.phylo(cdbs_tree, cex = 1.0, font = 3, edge.width = 2,
                label.offset = 0.0001, main = "Classical Distance Bootstrap Support (CDBS)")
ape::nodelabels(
  text = paste0(as.character(round(cdbs_numeric, 1)), ifelse(is.na(cdbs_numeric), "", "%")),
  adj = c(1, -0.5), frame = "none", cex = 0.8, font = 2, col = "darkred"
)
dev.off()
 
# Standalone NJ/CDBS tree
pdf("NJ_tree.pdf", width = 12, height = 15)
par(mar = c(5, 8, 4, 2))
ape::plot.phylo(cdbs_tree, cex = 1.0, font = 3, edge.width = 2,
                label.offset = 0.0001, main = "NJ Support (CDBS)")
ape::nodelabels(
  text = paste0(as.character(round(cdbs_numeric, 1)), ifelse(is.na(cdbs_numeric), "", "%")),
  adj = c(1, -0.5), frame = "none", cex = 0.8, font = 2, col = "darkblue"
)
dev.off()
 
# QBS vs CDBS scatterplot
pdf("QBS_vs_CDBS_Scatterplot.pdf", width = 7, height = 7)
valid_idx <- !is.na(qbs_numeric) & !is.na(cdbs_numeric)
plot(cdbs_numeric[valid_idx], qbs_numeric[valid_idx],
     xlab = "Classical Distance Bootstrap Support (CDBS %)",
     ylab = "Quantum Bootstrap Support (QBS %)",
     main = "QBS vs. CDBS Node Support",
     pch = 19, col = "steelblue", xlim = c(0, 100), ylim = c(0, 100))
abline(0, 1, col = "darkred", lty = 2, lwd = 2)
legend("topleft", legend = "1:1 line", col = "darkred",
       lty = 2, lwd = 2, box.col = "white")
dev.off()
 
# =============================================================================
# END OF PIPELINE
# =============================================================================
 
