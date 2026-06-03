#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#          Quantum Bootstrap Phylogenomics Pipeline
#               By John M. A. Wojahn, PhD, FLS
#                        3rd June 2026
#            Licensed under the GNU Affero v. 3
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

library(ape)
library(phangorn)
library(DECIPHER)
library(Biostrings)
library(treeio)
library(TreeDist)
library(PrideBar)
library(reticulate)

remotes::install("wojahn/PrideBar")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PRE-REQUISITE DEPENDENCY INSTALLATION SCRIPT
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# 1. Install CRAN-hosted core phylogenetics and infrastructure tools
cran_packages <- c("ape", "phangorn", "treeio", "TreeDist", "reticulate", "remotes")
new_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]

if(length(new_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(new_cran, collapse = ", "))
  install.packages(new_cran, repos = "https://cloud.r-project.org")
}

# 2. Install Bioconductor packages required for sequence alignment frameworks
bioc_packages <- c("DECIPHER", "Biostrings")
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]

if(length(new_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(new_bioc, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
}

# 3. Install GitHub-hosted custom tools
if (!requireNamespace("PrideBar", quietly = TRUE)) {
  message("Installing Custom Development Package: PrideBar from GitHub")
  remotes::install_github("wojahn/PrideBar", upgrade = "never")
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 1: TAXON SELECTION AND CORE ALIGNMENT PIPELINE
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

set.seed(35753)

plastomes <- ape::read.FASTA("AlignedOrderedGoodPlastomes_First.fasta")

# Extract the strategic 19-taxon cohort (16 Ingroup, 2 Benstonea, 1 Freycinetia)
plastomes <- plastomes[c(6,8,9,44,30,51,59,13,119,63,84,17,106,68,35,114,113,43,31)]

names(plastomes) <- c(
  "Freycinetia_baueriana", "Benstonea_copelandii",
  "Benstonea_herbacea", "Pandanus_furcatus", "Pandanus_calcis",
  "Pandanus_irregularis", "Pandanus_lacuum", "Pandanus_amaryllifolius",
  "Pandanus_vitiensis", "Pandanus_luzonensis", "Pandanus_polyglossus",
  "Pandanus_aquaticus", "Pandanus_tectorius", "Pandanus_maximus",
  "Pandanus_concretus", "Pandanus_utilis", "Pandanus_tsingycola",
  "Pandanus_eydouxia", "Pandanus_callmanderianus")

ape::write.FASTA(plastomes, "plastomes.fasta")

# Re-align raw profiles using DECIPHER framework
plastomez <- readDNAStringSet("plastomes.fasta")
plastomez <- DECIPHER::RemoveGaps(plastomez)
plastomez_a <- DECIPHER::AlignSeqs(plastomez, processors = 10)
plastomez_a <- DECIPHER::AdjustAlignment(plastomez_a)
Biostrings::writeXStringSet(plastomez_a, "Pandanus_a.fasta")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 2: CLASSICAL REFERENCE COMPUTATION (IQ-TREE)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

iqtree_cmd <- "iqtree -s Pandanus_a.fasta -m MFP -bb 1000 -nt AUTO -redo -o Freycinetia_baueriana -seed 35753 -safe"
system(iqtree_cmd)

# Parse optimal classical configuration
classical_ref <- read.tree("Pandanus_a.fasta.treefile")
classical_ref <- root(classical_ref, "Freycinetia_baueriana", resolve.root = TRUE)
classical_ref <- drop.tip(classical_ref, "Freycinetia_baueriana")
classical_ref <- ladderize(classical_ref)

pdf("IQTREE_tree.pdf", width = 15, height = 20)
plot(classical_ref, show.node.label = FALSE, cex = 1.0, root.edge = TRUE)
nodelabels(classical_ref$node.label, adj = c(1.2, -1.5), frame = "none", cex = 0.8, font = 2, col = "darkblue")
dev.off()

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 3: OPTIMIZED FREQUENTIST RESAMPLING LOOP (100 BOOTSTRAP GENERATION)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

final_alignment <- read.phyDat("Pandanus_a.fasta", format = "fasta")
boot_alignments <- bootstrap.phyDat(final_alignment, FUN = identity, bs = 100)
dir.create("bootstrap_inputs", showWarnings = FALSE)

PrideBar::SetPrideBar(1, 100, 2)
for (i in 1:100)
{
  PrideBar::PrideBar()
  rep_data <- boot_alignments[[i]]

  # Generate initial fast structural layout
  tmp_dm <- dist.ml(rep_data, model = "JC69")
  tmp_tree <- NJ(tmp_dm)
  tmp_fit <- pml(tmp_tree, data = rep_data)

  # Optimize parameters via NNI to fix bootstrap sampling noise before matrix output
  suppressMessages(
    tmp_opt <- optim.pml(tmp_fit, model = "GTR", k = 4,
                         optEdge = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE,
                         rearrangement = "NNI", control = pml.control(trace = 0))
  )

  # Extract the clean, mathematically validated hierarchical matrix
  rep_matrix <- as.matrix(cophenetic(tmp_opt$tree))

  # Save the finalized matrix
  file_name <- paste0("bootstrap_inputs/matrix_", sprintf("%03d", i), ".csv")
  write.csv(rep_matrix, file = file_name, row.names = TRUE)
}
PrideBar::FadePrideBar()

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 4: WRITE AND RUN QUANTUM BACKEND DISCRETE OPTIMIZATION
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# REAL Quantum version
dwave_python_script_content <- c(
  "import os",
  "import sys",
  "import glob",
  "import pandas as pd",
  "import numpy as np",
  "import qa_functions as qaf",
  "from dwave.system import DWaveSampler, EmbeddingComposite",
  "",
  "matrix_dir = 'bootstrap_inputs' if os.path.exists('bootstrap_inputs') else '../bootstrap_inputs'",
  "matrix_files = sorted(glob.glob(f'{matrix_dir}/matrix_*.csv'))",
  "",
  "if not matrix_files:",
  "    print(f\"Error: Cannot find matrix files in '{matrix_dir}'. Execution halted.\")",
  "    sys.exit(1)",
  "",
  "print(f\"Connecting to D-Wave Live QPU Hardware via Leap...\")",
  "try:",
  "    qpu_sampler = DWaveSampler(solver={'topology__type': 'pegasus'})",
  "    qpu_backend = EmbeddingComposite(qpu_sampler)",
  "    print(f\" Connected successfully to: {qpu_sampler.properties['id']}\")",
  "except Exception as e:",
  "    print(f\"Hardware Connection Error: {e}\")",
  "    print(\"Ensure your D-Wave API token is configured via environment variables or 'dwave config'.\")",
  "    sys.exit(1)",
  "",
  "print(f\"Found {len(matrix_files)} bootstrap matrices. Resolving on physical QPU...\")",
  "bootstrap_topologies = []",
  "",
  "for i, file_path in enumerate(matrix_files, 1):",
  "    df = pd.read_csv(file_path, index_col=0)",
  "    matrix = df.to_numpy()",
  "    taxa_labels = list(df.index)",
  "",
  "    max_dist = np.max(matrix) if np.max(matrix) > 0 else 1.0",
  "    similarity_matrix = max_dist - matrix",
  "    np.fill_diagonal(similarity_matrix, 0)",
  "",
  "    int_tags = list(range(len(taxa_labels)))",
  "    ",
  "    native_tree_object = qaf.sa_phylo_tree(",
  "        similarity_matrix, ",
  "        tags=int_tags, ",
  "        sampler=qpu_backend, ",
  "        num_reads=1000, ",
  "        annealing_time=20",
  "    )",
  "",
  "    newick_str = native_tree_object.to_newick(labels=taxa_labels)",
  "    bootstrap_topologies.append(newick_str)",
  "",
  "    if i % 5 == 0:",
  "        print(f\"  -> Processed {i}/100 profiles directly on the QPU.\")",
  "",
  "output_file = 'quantum_bootstrap_trees.tre'",
  "with open(output_file, 'w') as f:",
  "    for tree in bootstrap_topologies:",
  "        f.write(f\"{tree}\\n\")",
  "",
  "print(f\"\\n 100 QPU-annealed topologies populated and saved to '{output_file}'.\")"
)

writeLines(dwave_python_script_content, con = "run_quantum_bootstrap.py")

# SIMULATED Quantum version

# Define the fully corrected Python script content using the local classical simulation backend
dwave_python_script_content <- c(
  "import os",
  "import sys",
  "import glob",
  "import pandas as pd",
  "import numpy as np",
  "import qa_functions as qaf",
  "import neal",
  "",
  "# Verify paths across platforms",
  "matrix_dir = 'bootstrap_inputs' if os.path.exists('bootstrap_inputs') else '../bootstrap_inputs'",
  "matrix_files = sorted(glob.glob(f'{matrix_dir}/matrix_*.csv'))",
  "",
  "if not matrix_files:",
  "    print(f\"Error: Cannot find matrix files in '{matrix_dir}'. Execution halted.\")",
  "    sys.exit(1)",
  "",
  "print(\"Initializing classical Simulated Annealing simulator backend via neal...\")",
  "try:",
  "    # Initialize the local CPU simulator",
  "    sim_backend = neal.SimulatedAnnealingSampler()",
  "    print(\" Classical simulation environment established successfully.\")",
  "except Exception as e:",
  "    print(f\"Backend Initialization Error: {e}\")",
  "    sys.exit(1)",
  "",
  "print(f\"Found {len(matrix_files)} bootstrap matrices. Resolving via local heuristic sweeps...\")",
  "bootstrap_topologies = []",
  "",
  "for i, file_path in enumerate(matrix_files, 1):",
  "    df = pd.read_csv(file_path, index_col=0)",
  "    matrix = df.to_numpy()",
  "    taxa_labels = list(df.index)",
  "",
  "    # Invert matrix to build similarity edge arrays",
  "    max_dist = np.max(matrix) if np.max(matrix) > 0 else 1.0",
  "    similarity_matrix = max_dist - matrix",
  "    np.fill_diagonal(similarity_matrix, 0)",
  "",
  "    # Direct recursive top-down graph partition mapping",
  "    int_tags = list(range(len(taxa_labels)))",
  "    ",
  "    # Pass the local simulator engine along with sample read overrides",
  "    # num_reads=1000: 1000 thermodynamic Monte Carlo trajectories per split",
  "    native_tree_object = qaf.sa_phylo_tree(",
  "        similarity_matrix, ",
  "        tags=int_tags, ",
  "        sampler=sim_backend, ",
  "        num_reads=1000",
  "    )",
  "",
  "    # Convert structural tree object to labels-aware Newick text",
  "    newick_str = native_tree_object.to_newick(labels=taxa_labels)",
  "    bootstrap_topologies.append(newick_str)",
  "",
  "    if i % 5 == 0:",
  "        print(f\"  -> Processed {i}/100 profiles using local simulated annealing.\")",
  "",
  "# Save simulation payloads",
  "output_file = 'quantum_inspired_bootstrap_trees.tre'",
  "with open(output_file, 'w') as f:",
  "    for tree in bootstrap_topologies:",
  "        f.write(f\"{tree}\\n\")",
  "",
  "print(f\"\\n 100 simulated annealing topologies populated and saved to '{output_file}'.\")"
)

# Optional: Write the vector directly to disk from R to ensure a clean overwrite
writeLines(dwave_python_script_content, "run_quantum_bootstrap.py")

# Call Python backend via reticulate
use_condaenv("base", required = TRUE)
py_run_file("run_quantum_bootstrap.py")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 5: CONFIGURATION & INGESTION
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
target_tree_file  <- "Pandanus_a.fasta.treefile"
quantum_boot_file <- "quantum_bootstrap_trees.tre"

raw_ref <- ape::read.tree(target_tree_file)
raw_ref$tip.label <- chartr("-", "_", raw_ref$tip.label)
raw_ref$tip.label <- trimws(raw_ref$tip.label)
ref_labels <- sort(raw_ref$tip.label)

raw_ref_unrooted <- ape::unroot(raw_ref)
max_possible_splits <- length(raw_ref_unrooted$tip.label) - 3

boot_lines <- readLines(quantum_boot_file)
boot_lines <- boot_lines[boot_lines != ""]

clean_boot_list <- lapply(1:length(boot_lines), function(i)
{
  tr <- tryCatch(ape::read.tree(text = boot_lines[i]), error = function(e) NULL)
  if (!is.null(tr)) {
    tr$tip.label <- chartr("-", "_", tr$tip.label)
    tr$tip.label <- trimws(tr$tip.label)
    rogue_taxa <- tr$tip.label[!tr$tip.label %in% ref_labels]
    if (length(rogue_taxa) > 0) {
      tr <- ape::drop.tip(tr, rogue_taxa)
    }
  }
  return(tr)
})
clean_boot_list <- clean_boot_list[!sapply(clean_boot_list, is.null)]
class(clean_boot_list) <- "multiPhylo"

# Compute resolution metrics
resolution_percentages <- sapply(1:length(clean_boot_list), function(i)
{
  boot_tree <- ape::unroot(clean_boot_list[[i]])
  actual_splits <- boot_tree$Nnode
  if (actual_splits > max_possible_splits) actual_splits <- max_possible_splits
  return((actual_splits / max_possible_splits) * 100)
})

message("\n=== QUANTUM BOOTSTRAP RESOLUTION AUDIT ===")
print(summary(resolution_percentages))

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 6: CALCULATE QUANTUM BOOTSTRAP SUPPORT (QBS)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
supported_tree <- phangorn::plotBS(raw_ref, clean_boot_list, type = "none", p = 0)

# Drop outgroup to focus on ingroup support matrix visualization
supported_tree <- ape::root(supported_tree, "Freycinetia_baueriana", resolve.root = TRUE)
supported_tree <- ape::drop.tip(supported_tree, "Freycinetia_baueriana")
supported_tree <- ape::ladderize(supported_tree)
qbs_values <- supported_tree$node.label

# Render primary hybrid support tree
pdf("quantum_bootstrap_hybrid_tree.pdf", width = 12, height = 15)
par(mar = c(5, 8, 4, 2))
ape::plot.phylo(supported_tree, cex = 1.0, font = 3, edge.width = 2, label.offset = 0.0001)
ape::nodelabels(text = paste0(as.character(as.numeric(qbs_values)*100), ifelse(qbs_values == "", "", "%")),
                adj = c(1, -0.5), frame = "none", cex = 0.8, font = 2, col = "darkblue")
dev.off()


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 7: DENSITREE CLOUDOGRAM VISUALIZATION (Deep Outgroup Anchored)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Keep all 19 taxa intact (do not drop Freycinetia_baueriana)
rooted_boot_list <- lapply(1:length(clean_boot_list), function(i)
{
  tr <- clean_boot_list[[i]]

  # Root stably on the absolute outgroup tip across all 100 replicates
  tr_rooted <- ape::root(tr, outgroup = "Freycinetia_baueriana", resolve.root = TRUE)

  tr_bifurcated <- ape::multi2di(tr_rooted)
  tr_bifurcated$edge.length <- rep(1, nrow(tr_bifurcated$edge))
  tr_clean <- ape::ladderize(ape::read.tree(text = ape::write.tree(tr_bifurcated)))
  return(tr_clean)
})
class(rooted_boot_list) <- "multiPhylo"

# Re-ingest the reference tree from Phase 6 but KEEP Freycinetia for the background backbone
classical_ref_rooted <- phangorn::plotBS(raw_ref, clean_boot_list, type = "none", p = 0)
classical_ref_rooted <- ape::root(classical_ref_rooted, outgroup = "Freycinetia_baueriana", resolve.root = TRUE)
classical_ref_rooted <- ape::multi2di(classical_ref_rooted)
classical_ref_rooted$edge.length <- rep(1, nrow(classical_ref_rooted$edge))

pdf("Quantum_Bootstrap_Densitree.pdf", width = 12, height = 8)
par(mar = c(3, 3, 5, 10))

phangorn::densiTree(rooted_boot_list,
                    type = "cladogram",
                    alpha = 0.15,
                    col = "royalblue",
                    consensus = classical_ref_rooted,
                    direction = "rightwards",
                    scaleX = TRUE,
                    cex = 0.9,
                    font = 3,
                    main = "Quantum Bootstrap Densitree\n(Topological Variation Cloud vs. Classical Baseline)")
dev.off()

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 8: ROBINSON-FOULDS (RF) TOPOLOGICAL DISTANCE HISTOGRAM
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

rf_distances <- sapply(1:length(sanitized_boot_list), function(i)
{
  return(TreeDist::RobinsonFoulds(supported_tree, sanitized_boot_list[[i]]))
})

message("\n=== ROBINSON-FOULDS DISTANCE AUDIT ===")
print(summary(rf_distances))

pdf("Quantum_Bootstrap_RF_Histogram.pdf", width = 8, height = 6)
max_rf <- max(rf_distances, na.rm = TRUE)
if(max_rf == 0) max_rf <- 2
hist_breaks <- seq(-0.5, max_rf + 0.5, by = 1)

hist(rf_distances, breaks = hist_breaks, col = "skyblue", border = "white",
     main = "Topological Divergence Profile\n(Quantum Bootstrap Cloud vs. Classical Baseline)",
     xlab = "Absolute Robinson-Foulds (RF) Distance", ylab = "Frequency (Number of Replicates)",
     xaxt = "n", las = 1)
axis(1, at = 0:max_rf, labels = 0:max_rf)
abline(v = median(rf_distances), col = "darkred", lwd = 2, lty = 2)
legend("topright", legend = paste("Median RF Distance =", median(rf_distances)),
       col = "darkred", lwd = 2, lty = 2, box.col = "white")
dev.off()

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 9: MANUSCRIPT METRICS COMPILATION & SUMMARY TABLE EXPORT
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

res_summary <- summary(resolution_percentages)
rf_summary  <- summary(rf_distances)
raw_qbs_numeric <- as.numeric(qbs_values)
clean_qbs_pct   <- raw_qbs_numeric[!is.na(raw_qbs_numeric)] * 100
if(length(clean_qbs_pct) == 0) clean_qbs_pct <- 0

dataset_summary_matrix <- data.frame(
  Metric_Category = c(
    "Topological Resolution (%)", "Topological Resolution (%)", "Topological Resolution (%)",
    "Robinson-Foulds (RF) Distance", "Robinson-Foulds (RF) Distance", "Robinson-Foulds (RF) Distance",
    "Quantum Bootstrap Support (QBS)", "Quantum Bootstrap Support (QBS)", "Quantum Bootstrap Support (QBS)"
  ),
  Statistical_Parameter = c(
    "Minimum", "Median", "Maximum",
    "Minimum", "Median", "Maximum",
    "Minimum Node Support", "Median Node Support", "Maximum Node Support"
  ),
  Value = c(
    round(res_summary[["Min."]], 2), round(res_summary[["Median"]], 2), round(res_summary[["Max."]], 2),
    round(rf_summary[["Min."]], 2), round(rf_summary[["Median"]], 2), round(rf_summary[["Max."]], 2),
    paste0(round(min(clean_qbs_pct), 1), "%"), paste0(round(median(clean_qbs_pct), 1), "%"), paste0(round(max(clean_qbs_pct), 1), "%")
  ),
  stringsAsFactors = FALSE
)

n_tips <- length(supported_tree$tip.label)
internal_node_ids <- (n_tips + 1):(n_tips + supported_tree$Nnode)
formatted_pct_column <- sapply(raw_qbs_numeric, function(val)
{
  if (is.na(val)) return("Root/Unresolved Split")
  return(paste0(round(val * 100, 1), "%"))
})

node_support_table <- data.frame(
  Plot_Node_ID = internal_node_ids,
  Raw_Support_Proportion = supported_tree$node.label,
  Formatted_QBS_Percentage = formatted_pct_column,
  stringsAsFactors = FALSE
)

print(dataset_summary_matrix, row.names = FALSE)
write.csv(dataset_summary_matrix, file = "Quantum_Bootstrap_Dataset_Summary.csv", row.names = FALSE)
write.csv(node_support_table, file = "Quantum_Bootstrap_Node_Support_Details.csv", row.names = FALSE)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PHASE 10: CALCULATE LOW-SUPPORT NODAL DEPTH (MEASURING "MIDDLENESS")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

all_node_depths <- ape::node.depth.edgelength(supported_tree)
internal_depths <- all_node_depths[(n_tips + 1):length(all_node_depths)]
norm_depths <- (internal_depths - min(internal_depths)) / (max(internal_depths) - min(internal_depths))

support_depth_df <- data.frame(
  Node_ID = (n_tips + 1):length(all_node_depths),
  QBS = as.numeric(supported_tree$node.label) * 100,
  Relative_Depth = norm_depths
)

weak_nodes <- support_depth_df[!is.na(support_depth_df$QBS) & support_depth_df$QBS < 70, ]

if (nrow(weak_nodes) > 0)
{
  median_weak_depth <- median(weak_nodes$Relative_Depth)
  message("\n=== LOW-SUPPORT MIDDLENESS AUDIT ===")
  message("Number of low-support nodes (<70%): ", nrow(weak_nodes))
  message("Median Relative Depth of weak nodes: ", round(median_weak_depth, 4))

  if (median_weak_depth > 0.4 && median_weak_depth < 0.6)
  {
    message("Interpretation: Weak support is heavily clustered in the MID-SECTION of the tree topology.")
  } else if (median_weak_depth <= 0.4) {
    message("Interpretation: Weak support is clustered deep near the ROOT base (uwu).")
  } else {
    message("Interpretation: Weak support is clustered out near the peripheral TIPS.")
  }
} else {
  message("All nodes possess high support; no topological weak zones detected.")
}
