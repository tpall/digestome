# Causal DAG of a digester and of what a marker-gene profile observes (manuscript, "Causal structure of
# the readout"). One edge list drives the figure (figures/fig-dag.*) and the checks (dag/dag-checks.txt).
#
#   Rscript manuscript/dag/digester_dag.R
#
# State and community feed back on each other, so the graph is unrolled in time: *_0 is the
# history before the sample, the unsuffixed nodes are the reactor at sampling, *_1 the next interval.
suppressPackageStartupMessages(library(dagitty))
# Checked with dagitty 0.3.4. paths() stops at its `limit` argument, so the limit is set explicitly
# below and a truncated enumeration stops the script (review 2026-10-05, D1).
if (packageVersion("dagitty") < "0.3.4") stop("dagitty >= 0.3.4 expected")
PATH_LIMIT <- 10000
#
# Slices (review Q2): *_0 = the weeks before the sample, about one retention time (30-80 d);
# the unsuffixed nodes = the reactor in the sampling week; *_1 = the 8 weeks after the profile is
# delivered (the estimand window). The 3-month sampling interval is an outer loop: each sample
# starts a new set of slices.
#
# Considered and not drawn:
#   state -> extract  extraction yield depends on the matrix (solids, foam, VFA). Not drawn: if it is
#                     material, the state reaches the profile outside the community.
#   dna_0's own sampling, extraction and reference: drawn only for the current profile; the previous
#                     profile enters the graph only as a cause of the next sample.
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))

nodes <- read.table(header = TRUE, sep = "|", strip.white = TRUE, quote = "", text = "
id        | group      | label
feed      | lever      | Feed mix\\n(protein, fibre, fat)
olr       | lever      | Loading rate
hrt       | lever      | Retention time
temp      | lever      | Temperature
te_dose   | lever      | Trace-element\\ndosing
inoc      | lever      | Inoculum
air       | lever      | Air dosing\\n(desulphurisation)
reactor   | design     | Reactor type,\\nmixing, stages
state_0   | history    | Reactor state\\nbefore
comm_0    | history    | Community\\nbefore
dna_0     | report     | Previous\\nprofile
te        | state      | Trace elements\\nin the digester
acceptor  | state      | Sulfate, nitrate,\\noxygen
state     | state      | Reactor state\\n(VFA, H2, NH3, pH)
comm      | community  | Community\\n(who carries what)
sampling  | measure    | Sampling point\\nand time
extract   | measure    | DNA extraction,\\nsequencing depth
ref       | measure    | Reference genomes\\n(GTDB or own MAGs)
dna       | report     | Marker-gene profile\\n(genome shares)
proc      | observed   | Process data\\n(FOS/TAC, VFA, NH4, pH)
proc_1    | observed   | Process data,\\nnext interval
te_meas   | observed   | Trace-element\\nanalysis
gas       | observed   | Gas yield, CH4 %
cue       | hidden     | What the operator\\nsees on site (foam, smell)
goodday   | measure    | Sample on a\\nstable day
sampled   | selection  | A sample\\nis taken
decision  | lever      | Operator's next\\nchange (load, feed)
state_1   | outcome    | Reactor state,\\nnext interval
")

edges <- read.table(header = TRUE, sep = "|", strip.white = TRUE, quote = "", text = "
from      | to        | why
feed      | state     | protein -> NH3, fibre and fat -> VFA
olr       | state     | more substrate per day -> VFA
hrt       | state     | contact time
temp      | state     | rates; free NH3 share
te        | state     | cofactors of the methanogenic steps
acceptor  | state     | electrons diverted from methane
state_0   | state     | carry-over
comm_0    | state     | the community over the interval removes VFA, H2
feed      | te        | manure brings trace elements, silage little
te_dose   | te        | dosing
feed      | acceptor  | sulfate, nitrate in the feed
air       | acceptor  | oxygen from desulphurisation
comm_0    | comm      | the community persists (weeks to months)
feed      | comm      | substrate selection (hydrolysers follow the feed); immigration from manure
state     | comm      | selection: NH3, VFA, H2 favour some routes
temp      | comm      | selection by temperature
hrt       | comm      | washout of slow growers
te        | comm      | trace-element limitation
acceptor  | comm      | sulfate and nitrate reducers, methane oxidisers
inoc      | comm      | seeding
reactor   | hrt       | design fixes the range of retention times
reactor   | comm      | mixing, solids retention, stages
reactor   | sampling  | where a sample can be taken
comm      | dna       | what the profile reads
sampling  | dna       | point and time in the feeding cycle
extract   | dna       | extraction bias, depth sets the detection floor
ref       | dna       | which genomes reads are mapped to
state     | proc      | titration, VFA, NH4, pH
te        | te_meas   | digester trace-element analysis
comm      | gas       | activity
state     | gas       | inhibition, substrate
proc      | decision  | the operator reads the process data
dna       | decision  | the operator reads the profile
gas       | decision  | the operator watches the gas yield
state     | cue       | foam, smell, colour: the state seen directly
cue       | decision  | a change prompted by what was seen, often unrecorded
state     | sampled   | plants buy a sample when something is wrong
goodday   | sampled   | samples taken on stable days, not only after a disturbance
dna_0     | sampled   | a flagged profile brings the next sample sooner
comm_0    | dna_0     | the previous profile read the earlier community
state_1   | proc_1    | the outcome as the routine log records it
state     | state_1   | carry-over
comm      | state_1   | the community decides how the next change is absorbed
decision  | state_1   | the change itself
")

dag <- dagitty(paste0("dag {\n",
  paste(sprintf("%s [%s]", nodes$id,
    ifelse(nodes$group %in% c("report", "observed", "lever", "measure", "design", "selection"), "", "latent")),
    collapse = "\n"),
  "\n", paste(sprintf("%s -> %s", edges$from, edges$to), collapse = "\n"), "\n}"))
dag <- dagitty(gsub(" \\[\\]", "", as.character(dag)))
exposures(dag) <- "decision"; outcomes(dag) <- "state_1"
stopifnot(isAcyclic(dag))

# ---- drawing -------------------------------------------------------------------------------
style <- c(lever = "fillcolor=\"#E8F1F8\", color=\"#0072B2\"",
           design = "fillcolor=\"#F2F2F2\", color=\"#6D7173\"",
           history = "fillcolor=\"#FFFFFF\", color=\"#6D7173\", style=\"rounded,dashed,filled\"",
           state = "fillcolor=\"#FFFFFF\", color=\"#6D7173\", style=\"rounded,dashed,filled\"",
           community = "fillcolor=\"#FFF4E0\", color=\"#E69F00\", style=\"rounded,dashed,filled\", penwidth=2",
           measure = "fillcolor=\"#F2F2F2\", color=\"#6D7173\"",
           report = "fillcolor=\"#E69F00\", color=\"#9C6A00\", fontcolor=\"#000000\", penwidth=2",
           observed = "fillcolor=\"#E3F2EC\", color=\"#009E73\"",
           outcome = "fillcolor=\"#FBE7E1\", color=\"#D55E00\", style=\"rounded,dashed,filled\"",
           hidden = "fillcolor=\"#FFFFFF\", color=\"#0072B2\", style=\"rounded,dashed,filled\"",
           selection = "fillcolor=\"#FFFFFF\", color=\"#CC79A7\", shape=octagon, style=\"filled\"")
clusters <- list(
  c("Before the sample", "state_0", "comm_0", "dna_0"),
  c("Levers and design", "feed", "olr", "hrt", "temp", "te_dose", "inoc", "air", "reactor"),
  c("Reactor at sampling (not observed directly)", "te", "acceptor", "state", "comm"),
  c("What is measured", "sampling", "extract", "ref", "dna", "proc", "te_meas", "gas", "cue"),
  c("Which samples exist", "goodday", "sampled"),
  c("Next interval", "decision", "state_1", "proc_1"))
node_line <- function(id) {
  n <- nodes[nodes$id == id, ]
  sprintf("  %s [label=\"%s\", %s];", id, n$label, style[[n$group]])
}
dot <- c("digraph G {",
  "  rankdir=LR; newrank=true; nodesep=0.25; ranksep=0.55;",
  "  node [shape=box, style=\"rounded,filled\", fontname=\"Helvetica\", fontsize=11, margin=\"0.12,0.06\"];",
  "  edge [color=\"#6D7173\", arrowsize=0.6];",
  unlist(lapply(seq_along(clusters), function(i) {
    cl <- clusters[[i]]
    c(sprintf("  subgraph cluster_%d { label=\"%s\"; fontname=\"Helvetica\"; fontsize=12; color=\"#BBBBBB\"; style=\"rounded\";",
              i, cl[1]), sapply(cl[-1], node_line), "  }")
  })),
  sprintf("  %s -> %s%s;", edges$from, edges$to,
          ifelse(edges$to == "dna", " [color=\"#9C6A00\", penwidth=1.6]", "")),
  "}")
fig <- file.path(here, "..", "figures", "fig-dag")
writeLines(dot, file.path(here, "dag.dot"))
system2("dot", c("-Tsvg", file.path(here, "dag.dot"), "-o", paste0(fig, ".svg")))
system2("dot", c("-Tpdf", file.path(here, "dag.dot"), "-o", paste0(fig, ".pdf")))
system2("dot", c("-Tpng", "-Gdpi=200", file.path(here, "dag.dot"), "-o", paste0(fig, ".png")))

# ---- checks ----------------------------------------------------------------------------------
out <- c()
say <- function(...) out <<- c(out, sprintf(...))
lat <- latents(dag)
say("Nodes %d, edges %d, latent (not observed): %s", length(names(dag)), nrow(edges), paste(sort(lat), collapse = ", "))

say("\n## 1. What the profile has as causes (its parents and their ancestors)")
say("Parents of dna: %s", paste(parents(dag, "dna"), collapse = ", "))
say("Ancestors of comm: %s", paste(setdiff(ancestors(dag, "comm"), "comm"), collapse = ", "))

say("\n## 2. Effect of the operator's next change on the next state")
say("Parents of decision: %s", paste(parents(dag, "decision"), collapse = ", "))
sets <- adjustmentSets(dag, type = "minimal")
say("As drawn (what the operator saw on site is not recorded): %s",
    if (length(sets)) paste(sapply(sets, function(x) paste0("{ ", paste(x, collapse = ", "), " }")), collapse = "; ")
    else "no adjustment set exists")
rec <- dag; latents(rec) <- setdiff(latents(rec), "cue")
for (x in adjustmentSets(rec, type = "minimal"))
  say("If the reason for each change is recorded: { %s }", paste(x, collapse = ", "))
p <- paths(dag, "decision", "state_1", limit = PATH_LIMIT)
stopifnot(length(p$paths) < PATH_LIMIT)
back <- !grepl("^decision ->", p$paths)
say("Backdoor paths from decision to state_1: %d, of which %d open.", sum(back), sum(back & p$open))
say("Shortest open backdoor paths:")
short <- p$paths[back & p$open]; short <- short[order(nchar(short))]
for (x in head(short, 6)) say("  %s", x)

say("\n## 3. A change between two samples: what besides the community moves the profile")
say("Non-community parents of dna: %s", paste(setdiff(parents(dag, "dna"), "comm"), collapse = ", "))

say("\n## 4. Two plants side by side: open paths between plant-level causes and the profile")
for (x in setdiff(ancestors(dag, "dna"), c("dna", latents(dag)))) {
  n <- sum(grepl(" -> dna$", paths(dag, x, "dna", directed = TRUE, limit = PATH_LIMIT)$paths))
  if (n) say("  %s -> ... -> dna: %d directed path(s)", x, n)
}

say("\n## 5. Implications among observed nodes that involve a measured output")
obs <- setdiff(names(dag), lat)
keep <- Filter(function(ci) all(c(ci$X, ci$Y, ci$Z) %in% obs) &&
                 any(c("dna", "gas", "proc", "te_meas", "proc_1") %in% c(ci$X, ci$Y)),
               impliedConditionalIndependencies(dag))
fmt <- function(ci) trimws(paste(capture.output(print(ci)), collapse = " "))
survives <- sapply(keep, function(ci) dseparated(dag, ci$X, ci$Y, unique(c(ci$Z, "sampled"))))
say("Testable on collected field samples (they still hold given that a sample was taken): %d", sum(survives))
for (ci in keep[survives]) say("  %s", fmt(ci))
say("Not testable on collected field samples: opened by selection on sampled = yes: %d", sum(!survives))
for (ci in keep[!survives]) say("  %s", fmt(ci))

say("\n## 6. Which samples exist")
say("Parents of sampled: %s. Every analysis of collected samples conditions on sampled = yes.",
    paste(parents(dag, "sampled"), collapse = ", "))
say("Without samples on stable days, sampled depends on the state (and the previous profile) alone: the samples")
say("over-represent upset states. Section 5 lists the implications this selection breaks.")

writeLines(out, file.path(here, "dag-checks.txt"))
cat(out, sep = "\n")
