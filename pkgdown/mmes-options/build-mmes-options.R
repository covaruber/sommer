# Builds pkgdown/assets/mmes-options.html, the interactive guide to
# mmes() options. Run from the package root with jsonlite installed:
#   Rscript pkgdown/mmes-options/build-mmes-options.R [extra output path ...]
# Option metadata (capability catalog, defaults, compatibility table and
# fit-time error messages) is read from the working-tree R sources, so rebuild
# the page whenever those change.

stopifnot(file.exists("DESCRIPTION"), requireNamespace("jsonlite", quietly=TRUE))
src <- new.env(parent=baseenv())
for(file in c("R/FUN_mme_control.R", "R/FUN_mme_capabilities.R", "R/FUN_mme_modules.R")){
  sys.source(file, envir=src, keep.source=FALSE)
}

capabilities <- src$.mmes_capabilities()
build <- src$.mmes_control_build
controls <- capabilities$controls
defaults <- build()

# Check the catalog against the public API, including nested PQL fields.
leaf_paths <- function(x, prefix="control"){
  unlist(lapply(names(x), function(name){
    path <- paste(prefix, name, sep=".")
    if(is.list(x[[name]]) && length(x[[name]])) leaf_paths(x[[name]], path) else path
  }), use.names=FALSE)
}
current <- leaf_paths(defaults)
controls <- controls[!controls$path %in% "control.output.pev", ]
catalog_controls <- controls$path[startsWith(controls$path, "control.")]
if(!setequal(current, catalog_controls)){
  stop("Control catalog differs from mmes_control(); update the catalog: ",
    paste(union(setdiff(current, catalog_controls), setdiff(catalog_controls, current)), collapse=", "))
}
for(i in which(startsWith(controls$path, "control."))){
  parts <- strsplit(sub("^control\\.", "", controls$path[i]), ".", fixed=TRUE)[[1]]
  value <- Reduce(function(x, part) x[[part]], parts, defaults)
  controls$default[i] <- paste(deparse(value), collapse="")
}
controls$default[controls$path %in% c("control.convergence.maxTraceProbes",
  "control.convergence.maxLanczosSteps")] <- "NULL (resolves to initial probes / Lanczos steps)"
controls$values[controls$path == "control.execution.memoryPolicy"] <- "speed, limit."
controls$values[controls$path == "control.execution.denseMemoryMB"] <- "Positive finite MB budget."
controls$values[controls$path == "control.execution.gpu"] <- "TRUE or FALSE; requires a CUDA-enabled build."
controls$values[controls$path == "control.execution.gpuMinSize"] <- "Positive integer."
controls$values[controls$path == "control.task.pql.estimateTheta"] <- "TRUE or FALSE."
controls$values[controls$path == "control.output.inverse"] <-
  "0: no fit-time PEVs or inverse; 1: selected inverse and PEVs; 2: full inverse and PEVs."
controls$effect[controls$path == "control.output.inverse"] <-
  "Controls Gaussian fit-time PEV work. Native Laplace requires mode 0 and computes PEV only after fitting through postPEV(); Laplace vcov()/summary() uncertainty is also post-fit. Prediction solve panels use output$panelSize and execution$denseMemoryMB."

public_path <- function(path){
  path <- sub("^control\\.", "", path)
  if(grepl("^task\\.glmm\\.", path)) return(sub("^task\\.glmm\\.", "task$glmm$", path))
  if(grepl("^task\\.pql\\.", path)) return(sub("^task\\.pql\\.", "task$pql$", path))
  sub(".", "$", path, fixed=TRUE)
}
controls$option <- vapply(controls$path, function(path)
  if(grepl("^control\\.", path)) public_path(path) else sub("^(mmes|model)\\.", "", path), "")
controls$kind <- ifelse(grepl("^control\\.", controls$path), "mmes_control()",
  ifelse(grepl("^mmes\\.", controls$path), "mmes() argument", "model formula"))
controls$usage <- vapply(controls$path, function(path){
  if(path == "control.formulation") return('mmes_control(formulation = "direct")')
  if(grepl("^control\\.task\\.pql\\.", path))
    return(sprintf("mmes_control(task = list(pql = list(%s = ...)))", sub("^control\\.task\\.pql\\.", "", path)))
  if(grepl("^control\\.task\\.glmm\\.", path))
    return(sprintf("mmes_control(task = list(glmm = list(%s = ...)))", sub("^control\\.task\\.glmm\\.", "", path)))
  if(grepl("^control\\.", path)){
    parts <- strsplit(sub("^control\\.", "", path), ".", fixed=TRUE)[[1]]
    return(sprintf("mmes_control(%s = list(%s = ...))", parts[1], parts[2]))
  }
  if(grepl("^mmes\\.", path)) return(sprintf("mmes(..., %s = ...)", sub("^mmes\\.", "", path)))
  switch(path,
    model.covariance="random = ~ vsm(usm(env), ism(id))",
    model.rotation="random = ~ vsm(ism(id), Gu = Ai, rotation = TRUE)",
    model.latentRR="random = ~ vsm(rrm(env, nPC = 2), ism(id))",
    model.Gu="random = ~ vsm(ism(id), Gu = Ai)", "")
}, "")
catalog_path <- function(option){
  path <- paste0("control.", gsub("$", ".", option, fixed=TRUE))
  path[option == "formulation"] <- "control.formulation"
  path
}
present <- function(paths){
  unknown <- setdiff(paths, controls$path)
  if(length(unknown)) stop("Unknown guide options: ", paste(unknown, collapse=", "))
  paths
}

# Settings each engine actually reads are derived from the compatibility table.
rules_src <- src$.mmes_compatibility_rules()
rule <- function(id) rules_src[[which(vapply(rules_src, `[[`, "", "id") == id)]]
allowed_under <- function(id, control){
  everything <- src$.mme_control_paths(control)
  bad <- names(rule(id)$bad)
  sort(union(setdiff(everything, rule(id)$settings(control)), bad))
}
direct_allowed <- setdiff(allowed_under("direct-engine", build(formulation="direct")),
  c("task$mode", "task$pql", "task$glmm",
    grep("^diagnostics\\$", src$.mme_control_paths(defaults), value=TRUE)))
solve_allowed <- setdiff(allowed_under("solve-mode", build(task=list(mode="solve"))),
  c("task$mode", "task$glmm", grep("^(diagnostics|execution)\\$", src$.mme_control_paths(defaults), value=TRUE),
    "solve$blockSchur", "solve$blockBanded", "solve$factorSchurMultiFactor", "solve$dFirstML"))
pcg_only <- src$.mme_pcg_paths

# ---------------------------------------------------------------- decision tree
node <- function(id, title, subtitle, question, text, controls, rules=character(), example="")
  list(id=id, title=title, subtitle=subtitle, question=question, text=text,
    controls=present(controls), rules=rules, example=example)
nodes <- list(
  node("model", "Model specification", "Formulas, data, covariance structures",
    "What is the statistical model?",
    "Formulas, covariance constructors, relationship matrices, weights and missing-data policy define the model. They are mmes() arguments rather than numerical options, but they decide which structure paths are eligible further down the tree (see the Structures tab).",
    c(grep("^mmes\\.", controls$path, value=TRUE), "model.covariance", "model.Gu"),
    example='fit <- mmes(Yield ~ Env, random = ~ vsm(ism(Name), Gu = Ai),\n            rcov = ~ units, data = DT, control = mmes_control())'),
  node("task", "Task", "Estimate covariances or solve at known values",
    "Do you need to estimate covariance parameters?",
    "task$mode = 'fit' (default) estimates covariance parameters. task$mode = 'solve' runs a single Henderson PCG solve at covariance parameters supplied through covPar.",
    c("control.task.mode", "control.task.glmm.method", "mmes.covPar"), c("solve-mode"),
    'mmes_control(task = list(mode = "fit"))'),
  node("known", "Known-covariance solve", "task$mode = 'solve' | one PCG solve",
    "Solve the mixed model equations at covPar.",
    "Always Henderson PCG. formulation, estimation modules and non-PCG solvers are not used and stop the call when supplied. Only the settings listed below (plus execution and diagnostics) apply.",
    c(catalog_path(solve_allowed), "mmes.covPar"), c("solve-mode"),
    'fit0 <- mmes(..., control = mmes_control())\nmmes(..., covPar = fit0,\n     control = mmes_control(task = list(mode = "solve"),\n                            solve = list(tol = 1e-10, warmStart = TRUE)))'),
  node("family", "Response family", "Gaussian identity or GLMM",
    "Is the response Gaussian with identity link?",
    "family selects the statistical response model, not a solver. Gaussian identity fits go straight to the equation domain. Other families use PQL by default, or the opt-in Laplace backend for its narrower supported model set.",
    "mmes.family", c("gaussian-pql"), 'mmes(y ~ x, random = ~ id, data = d, family = binomial())'),
  node("pql", "Outer PQL loop", "Non-Gaussian or non-identity link",
    "How should the outer PQL iterations run?",
    "task$pql owns outer PQL convergence; convergence$* owns the inner Gaussian working fits. Changing task$pql on a Gaussian identity fit is an error because it would be ignored.",
    grep("^control\\.task\\.pql\\.", controls$path, value=TRUE), c("gaussian-pql"),
    'mmes_control(task = list(pql = list(maxit = 40L, tol = 1e-6)))'),
  node("laplace", "Native Laplace GLMM", "Binomial or Poisson | ML | native covariance module",
    "Does this GLMM fit the first Laplace backend?",
    "Laplace uses sommer's native likelihood/mode and full-score modules with prepared designs and CovarianceFactor random precisions. It supports binomial logit/probit/cloglog (binary or grouped), Poisson log, NB2 log, Gaussian identity and Gamma log, including familym with conditionally independent observations. Shape/dispersion are estimated per trait by default; configure familyParameters to fix them. It requires ML and Henderson formulation. Structured residual covariance, W, weights formulas, vcc and direct formulation remain unsupported. Native Laplace requires output$inverse=0; use postPEV() after fitting for PEV. vcov(), summary() and prediction uncertainty also run only after fitting. output$panelSize and execution$denseMemoryMB control post-fit solve work. Inspect curvature/boundary diagnostics before inference. Both Gaussian engine objectives remain unchanged.",
    grep("^control\\.task\\.glmm\\.", controls$path, value=TRUE),
    example='mmes_control(task = list(objective = "ml", glmm = list(method = "laplace")))'),
  node("domain", "Equation domain", "formulation = henderson | direct",
    "Coefficient-space (Henderson) or observation-space (direct V) equations?",
    "Henderson (default) solves the mixed model equations in coefficient space and offers solver, representation and curvature choices. Direct inverts the n x n observation covariance V and accepts a much smaller set of options.",
    "control.formulation", c("direct-engine"), 'mmes_control(formulation = "direct")'),
  node("direct", "Direct V-inversion engine", "formulation = 'direct' | ML/REML",
    "Observation-space AI-REML/ML with EM/PXEM recovery.",
    "These are the only settings the direct engine reads (plus diagnostics). Anything else, such as solver choice, representations, curvature variants, PCG estimators or Henderson-only outputs, stops the call when supplied.",
    catalog_path(direct_allowed), c("direct-engine"),
    'mmes_control(formulation = "direct",\n             update = list(recovery = "pxem", blendWeight = 0.3),\n             proposal = list(stepWeight = 0.5))'),
  node("henderson", "Henderson backend", "solve$method",
    "Which linear-algebra backend solves the Henderson equations?",
    "'auto' (default) picks cholmod when any supplied relationship inverse has density > 0.2 and ldlt otherwise. It never picks pcg; PCG must be requested explicitly.",
    "control.solve.method", c("deterministic-solver", "pcg-estimation"),
    'mmes_control(solve = list(method = "cholmod"))'),
  node("pcg", "PCG stochastic REML", "solve$method = 'pcg' | REML only",
    "Iterative solves, SLQ log-determinant and Hutchinson score traces.",
    "Every setting listed here requires solve$method = 'pcg' and stops the call otherwise. PCG is REML-only and rejects EM/PXEM recovery, blending, common-scale, reduced AI, Aitken, hybrid AI, saddle-point equations and factor scores.",
    catalog_path(pcg_only), c("deterministic-solver", "pcg-estimation"),
    'mmes_control(solve = list(method = "pcg", tol = 1e-8, preconditioner = "block-jacobi"),\n             score = list(traceProbes = 32L), logdet = list(lanczosSteps = 40L),\n             update = list(blendWeight = 0))'),
  node("deterministic", "Deterministic factorization", "auto / ldlt / cholmod | ML or REML",
    "Exact factorization with optional structured elimination.",
    "Exact log-determinants and traces. Structured-elimination toggles take 'auto', 'off' or 'force'; 'force' never bypasses structural or memory checks.",
    c("control.task.objective", "control.solve.blockSchur", "control.solve.blockBanded",
      "control.solve.factorSchurMultiFactor", "control.solve.dFirstML",
      "control.score.covarianceChecks", "control.score.traceBatch",
      "control.score.denseCrossproduct"),
    c("deterministic-solver"),
    'mmes_control(solve = list(method = "cholmod", blockSchur = "force"),\n             task = list(objective = "ml"))'),
  node("representation", "Equation representation", "Compression, factor scores, saddle-point, conditional",
    "Can an exact equivalent representation make the equations cheaper?",
    "These are alternatives, not a sequence. Each has its own eligibility and several exclude each other. Conditional equations use sparse latent FA equations with AI or EM updates, without a retained public coefficient matrix or covariance uncertainty. Rotation and latent reduced rank are requested from the model formula.",
    c("control.representation.compression", "control.representation.factorScores",
      "control.representation.fixedRank",
      "control.representation.equations", "control.representation.sparseFactors",
      "model.rotation", "model.latentRR"),
    c("saddle-point", "factor-scores", "conditional-fa"),
    'mmes_control(representation = list(compression = "classes"))'),
  node("curvature", "Curvature", "Exact, hybrid or reduced AI",
    "How is the average-information matrix built?",
    "exact-ai (default) refreshes AI every iteration. hybrid-ai uses damped BFGS between exact refreshes. curvature$reduced assembles AI only in the independent directions of a reducing vcc map.",
    c("control.curvature.method", "control.curvature.reduced"),
    c("accelerated-proposal", "hybrid-aitken"),
    'mmes_control(curvature = list(method = "hybrid-ai"))'),
  node("optimizer", "Update and proposal", "Recovery, blending, extrapolation",
    "How is each covariance update proposed and accepted?",
    "update owns the direction: AI with EM/PXEM recovery, blending or common-scale profiling. Conditional equations support AI or EM; conditional-fa is the legacy conditional EM selector. proposal owns scaling and Aitken extrapolation. vcc ties covariance parameters.",
    c("control.update.method", "control.update.recovery", "control.update.blendWeight", "control.update.commonScale",
      "control.proposal.stepWeight", "control.proposal.extrapolation", "mmes.vcc"),
    c("accelerated-proposal", "blending", "hybrid-aitken", "conditional-fa"),
    'mmes_control(update = list(recovery = "pxem", blendWeight = c(1, 0.5, 0)))'),
  node("convergence", "Convergence", "Iterations and tolerances",
    "When does the optimizer stop?",
    "maxIters and the likelihood and parameter tolerances apply to every estimation path. Adaptive accuracy and stochastic tolerances belong to PCG (see that node).",
    c("control.convergence.maxIters", "control.convergence.likelihoodTol", "control.convergence.parameterTol"),
    example='mmes_control(convergence = list(maxIters = 50L, likelihoodTol = 1e-6))'),
  node("outputs", "Outputs", "Coefficient and covariance uncertainty",
    "Which inference products should be returned?",
    "output$inverse = 0 skips fit-time PEVs and coefficient inversion; 1 computes PEVs; 2 computes the full inverse and PEVs. To defer this work, Gaussian fits need coefficientMatrix = TRUE, while native Laplace fits retain their joint-information inputs automatically; call postPEV(fit, mode = 1L or 2L). There is no separate output$pev control. Covariance-parameter information is a separate request.",
    grep("^control\\.output\\.", controls$path, value=TRUE), c("coefficient-matrix"),
    'mmes_control(output = list(inverse = 1L))'),
  node("policy", "Execution and diagnostics", "Memory, threads, GPU, messages",
    "Which hardware and reporting policies apply?",
    "Execution settings never change the statistical model. GPU requests need a CUDA build and are Henderson-only; direct and solve tasks stay on CPU.",
    grep("^control\\.(execution|diagnostics)\\.", controls$path, value=TRUE),
    example='mmes_control(execution = list(blasThreads = 4L, denseMemoryMB = 4000),\n             diagnostics = list(verbose = FALSE))'))

# Any catalog entry not placed explicitly goes to its module's home node.
home <- c(model="model", task="task", formulation="domain", representation="representation",
  solve="deterministic", logdet="pcg", score="deterministic", curvature="curvature",
  update="optimizer", proposal="optimizer", convergence="convergence", output="outputs",
  execution="policy", diagnostics="policy")
placed <- unique(unlist(lapply(nodes, `[[`, "controls")))
for(i in which(!controls$path %in% placed)){
  target <- home[[controls$module[i]]]
  k <- which(vapply(nodes, `[[`, "", "id") == target)
  nodes[[k]]$controls <- c(nodes[[k]]$controls, controls$path[i])
  message("Placed ", controls$path[i], " under ", target)
}
controls$nodes <- lapply(controls$path, function(path)
  vapply(Filter(function(n) path %in% n$controls, nodes), `[[`, "", "id"))

edges <- list(
  c("model", "task", "", "branch"),
  c("task", "known", "mode = 'solve'", "branch"),
  c("task", "family", "mode = 'fit'", "branch"),
  c("family", "pql", "non-Gaussian", "branch"),
  c("family", "domain", "Gaussian identity", "branch"),
  c("pql", "domain", "inner working fit", "join"),
  c("domain", "direct", "'direct'", "branch"),
  c("domain", "henderson", "'henderson'", "branch"),
  c("henderson", "deterministic", "auto / ldlt / cholmod", "branch"),
  c("henderson", "pcg", "'pcg'", "branch"),
  c("deterministic", "representation", "", "branch"),
  c("representation", "curvature", "", "branch"),
  c("curvature", "optimizer", "", "branch"),
  c("optimizer", "convergence", "", "branch"),
  c("pcg", "convergence", "stochastic AI", "join"),
  c("direct", "convergence", "AI + EM/PXEM", "join"),
  c("convergence", "outputs", "", "branch"),
  c("known", "outputs", "no estimation", "join"),
  c("outputs", "policy", "", "branch"))
edges <- lapply(edges, function(e) list(from=e[1], to=e[2], label=e[3], kind=e[4]))

# ---------------------------------------------------------- compatibility rules
contexts <- list(
  `direct-engine`=list(when="formulation = 'direct' (fit task)", nodes=c("domain", "direct"),
    control=build(formulation="direct")),
  `solve-mode`=list(when="task$mode = 'solve'", nodes=c("task", "known"),
    control=build(task=list(mode="solve"))),
  `deterministic-solver`=list(when="Henderson fit with solve$method 'auto', 'ldlt' or 'cholmod'",
    nodes=c("henderson", "deterministic", "pcg"), control=defaults),
  `pcg-estimation`=list(when="Henderson fit with solve$method = 'pcg'", nodes=c("henderson", "pcg"),
    control=build(solve=list(method="pcg"))),
  `dense-score`=list(when="PCG, compressed, saddle-point or conditional equations",
    nodes=c("henderson", "deterministic", "representation"),
    control=build(solve=list(method="pcg"))),
  `accelerated-proposal`=list(when="proposal$extrapolation = 'aitken' or curvature$method = 'hybrid-ai'",
    nodes=c("curvature", "optimizer"), control=build(curvature=list(method="hybrid-ai"))),
  `saddle-point`=list(when="representation$equations = 'saddle-point'", nodes="representation",
    control=build(representation=list(equations="saddle-point"))),
  `factor-scores`=list(when="representation$factorScores is not 'none'", nodes="representation",
    control=build(representation=list(factorScores="fixed-shape"))),
  blending=list(when="any update$blendWeight > 0", nodes="optimizer",
    control=build(update=list(blendWeight=1))),
  `coefficient-matrix`=list(when="Ordinary equations (not conditional FA), with representation$equations not 'saddle-point' or output$inverse not 0",
    nodes="outputs", control=defaults),
  `conditional-fa`=list(when="representation$equations = 'conditional' or update$method = 'conditional-fa'",
    nodes=c("representation", "optimizer", "outputs"),
    control=build(representation=list(equations="conditional"))))
rule_ids <- vapply(rules_src, `[[`, "", "id")
if(!setequal(rule_ids, names(contexts))) stop("Compatibility rules changed; update contexts in this script: ",
  paste(union(setdiff(rule_ids, names(contexts)), setdiff(names(contexts), rule_ids)), collapse=", "))
conditional <- c(`output$inverse`="only output$inverse = 1",
  `solve$method`="only 'ldlt' or 'cholmod'", `solve$preconditioner`="only 'nystrom'")
rules <- lapply(rules_src, function(r){
  context <- contexts[[r$id]]
  stopifnot(isTRUE(r$when(context$control)))
  bad <- names(r$bad)
  label <- function(path) if(path %in% bad) paste0(path, " (", conditional[[path]], ")") else path
  list(id=r$id, stage="mmes_control() and mmes()", when=context$when, nodes=context$nodes,
    message=r$message, hint=r$hint,
    offending=vapply(unique(r$settings(context$control)), label, "", USE.NAMES=FALSE))
})
rules <- c(rules, list(
  list(id="hybrid-aitken", stage="mmes_control()",
    when="curvature$method = 'hybrid-ai' together with proposal$extrapolation = 'aitken'",
    nodes=c("curvature", "optimizer"),
    message="Hybrid AI and Aitken extrapolation are separate accelerations and cannot be requested together.",
    hint="Choose one of them.", offending=c("curvature$method", "proposal$extrapolation")),
  list(id="laplace-eligibility", stage="mmes()",
    when="task$glmm$method = 'laplace'",
    nodes=c("family", "laplace"),
    message="Native Laplace requires objective='ml', a supported family/link (including grouped trials or familym), and Henderson formulation.",
    hint="Use PQL for unsupported residual covariance, W, weights formulas, vcc or direct formulation.",
    offending="task$glmm$method"),
  list(id="gaussian-pql", stage="mmes() (needs the family)", when="Gaussian identity family",
    nodes=c("family", "pql"),
    message="task$pql controls the outer PQL loop of non-Gaussian fits; a Gaussian identity fit would ignore it.",
    hint="Remove task$pql or use a non-Gaussian family.", offending="task$pql")))

# -------------------------------------------------------------- fit-time checks
module_pattern <- "(formulation|task|representation|solve|logdet|score|curvature|update|proposal|convergence|output|execution|diagnostics)\\$[A-Za-z]+"
fit_checks <- list()
collect <- function(x, file){
  if(!is.call(x)) return(invisible())
  if(identical(x[[1]], as.name("stop"))){
    args <- as.list(x)[-1]
    if(length(names(args))) args <- args[names(args) != "call."]
    text <- paste(vapply(args, function(z) if(is.character(z)) z else "...", ""), collapse="")
    if(grepl(module_pattern, text) && !grepl("^Incompatible", text)){
      options <- unique(regmatches(text, gregexpr(module_pattern, text))[[1]])
      fit_checks[[length(fit_checks) + 1L]] <<- list(file=file, message=text, options=options)
    }
  }
  for(z in as.list(x)[-1]) if(!missing(z)) collect(z, file)
}
for(file in list.files("R", "\\.R$", full.names=TRUE)){
  for(expr in parse(file, keep.source=FALSE)) collect(expr, basename(file))
}
fit_checks <- fit_checks[!duplicated(vapply(fit_checks, `[[`, "", "message"))]

# --------------------------------------------- structure paths and rules of thumb
structure <- function(pattern, example, path, options, requires, nodes)
  list(pattern=pattern, example=example, path=path, options=options, requires=requires, nodes=nodes)
structures <- list(
  structure("Independent residuals", "rcov = ~ units",
    "R^-1 is applied elementwise and log|R| is a sum of log diagonals; no residual factorization.",
    "automatic", "Diagonal residual covariance, optional diagonal W.", "model"),
  structure("Complete repeated residual blocks (multi-trait)", "rcov = ~ vsm(usm(trait), ism(units))",
    "The small block covariance is factored once; its inverse and log-determinant are reused across blocks and traces accumulate block by block.",
    "automatic", "Every block complete. Missing cells or irregular blocks fall back to sparse residual factorization.", "model"),
  structure("Kronecker random covariance", "random = ~ vsm(usm(env), ism(id), Gu = A)",
    "Factor-wise inverse, log-determinant and derivatives (exact algebra). AR(1), compound symmetry and ANTE use native closed forms.",
    "automatic", "Product covariance built from vsm() factors.", "model"),
  structure("Sparse relationship precision (pedigree)", "random = ~ vsm(ism(id), Gu = Ainv)",
    "solve$method = 'auto' selects ldlt (density <= 0.2). output$inverse = 1 uses the Takahashi selected inverse.",
    "solve$method, output$inverse", "attr(Gu, 'inverse') = TRUE and matching level names.", c("henderson", "deterministic")),
  structure("Dense relationship precision (genomic)", "random = ~ vsm(ism(id), Gu = Ginv)",
    "solve$method = 'auto' selects cholmod (density > 0.2). Eligible graphs may use exact block-Schur elimination; GPU kernels can help large dense blocks.",
    "solve$method, solve$blockSchur, execution$denseMemoryMB, execution$gpu", "Workspace within execution$denseMemoryMB.", c("henderson", "deterministic", "policy")),
  structure("Few records relative to random-effect levels", 'mmes_control(formulation = "direct")',
    "Direct inversion works on the n x n observation covariance V, so cost follows the number of records rather than the number of effects.",
    "formulation", "Dense V must fit in memory and stay positive definite.", c("domain", "direct")),
  structure("Very large systems where factor memory is the bottleneck",
    'mmes_control(solve = list(method = "pcg"), update = list(blendWeight = 0))',
    "Iterative PCG with stochastic log-determinant (SLQ) and Hutchinson score traces. REML only.",
    "solve$*, logdet$*, score$traceProbes, convergence$accuracy", "Explicit solve$method = 'pcg'.", "pcg"),
  structure("PCG with factor-wise (tensor) random precision",
    "random = ~ vsm(csm(env), dsm(trait), ism(id)) with solve$method = 'pcg'",
    "The random prior is applied through Kronecker operations instead of being assembled into C each iteration.",
    "representation$matrixFree", "PCG, REML, output$inverse = 0 and compatible factor-wise precision.", "pcg"),
  structure("AR(p), ANTE, SAR or CAR factors", "random = ~ vsm(ar1m(time), ism(id))",
    "Exact sparse factor operands for precision and derivatives.",
    "representation$sparseFactors = 'force'", "Eligible structured factor. It is not a new covariance model.", "representation"),
  structure("Factor-analytic (low rank + diagonal) or compound symmetry", "random = ~ vsm(usm(fa(env, 2)), ism(id))",
    "Equivalent latent factor-score representation of the Henderson equations.",
    "representation$factorScores ('fixed-shape' or 'profile')",
    "Deterministic Henderson REML, output$inverse = 0; no vcc, rotation, compression or positive blending.", "representation"),
  structure("Zero-specific reduced rank", "random = ~ vsm(rrm(env, nPC = 2), ism(id))",
    "Native latent reduced-rank fitting (Henderson) or record-space covariance (direct).",
    "automatic", "Henderson needs a deterministic solver; excludes PCG, compression, hybrid AI, saddle-point, factor scores, reduced AI and task$mode = 'solve'.", c("model", "representation")),
  structure("Balanced layout with one relationship matrix", "random = ~ vsm(ism(id), Gu = A, rotation = TRUE)",
    "Eigenbasis rotation; coefficients are returned in original coordinates. Available on Henderson and direct.",
    "vsm(rotation = TRUE), output$inverse = 0", "Complete balanced Gaussian layout, invariant residuals and block-constant diagonal W.", c("model", "representation")),
  structure("Many records sharing identical design rows", 'mmes_control(representation = list(compression = "classes"))',
    "Exact sufficient-statistic reduction; 'classes' also caches cross-products and derivative operands.",
    "representation$compression", "Gaussian Henderson, independent residuals and fixed diagonal weights.", "representation"),
  structure("Multi-factor latent low-rank layout with dense relationship", "eligible FA layouts",
    "Multi-factor Woodbury/Schur elimination.", "solve$factorSchurMultiFactor",
    "REML, disjoint groups, dense-enough relationship and workspace; blockSchur admission.", "deterministic"),
  structure("Block-banded coefficient graph (bandwidth <= 4)", 'mmes_control(solve = list(method = "cholmod", blockBanded = "force"))',
    "Exact block-banded elimination.", "solve$blockBanded = 'force'", "REML CHOLMOD graph; only 'force' enables it.", "deterministic"),
  structure("ML fit with random coefficients", 'mmes_control(task = list(objective = "ml"))',
    "Eliminates the random block before the fixed-effect Schur complement.", "solve$dFirstML",
    "Deterministic Henderson ML, output$inverse other than 1.", "deterministic"),
  structure("Tied covariance parameters (vcc)", "vcc = data.frame(parameter = c(2L, 3L), group = 1L)",
    "AI is assembled only in independent directions; scores stay in full coordinates.",
    "curvature$reduced, update$commonScale", "Reducing vcc map, deterministic Henderson, Gaussian.", c("curvature", "optimizer")),
  structure("Non-Gaussian response", "family = binomial()",
    "Outer PQL loop around Gaussian working fits.", "task$pql", "Inner fits obey the Gaussian-path restrictions.", c("family", "pql")))

tip <- function(topic, text, nodes) list(topic=topic, text=text, nodes=nodes)
tips <- list(
  tip("Start simple", "Begin with mmes_control() defaults and change one module at a time. Restating a default is always accepted; a non-default setting that the selected path would ignore stops the call and names the offending options.", "model"),
  tip("Solver", "Use solve$method = 'auto' (or 'ldlt') for sparse pedigree-like relationship precision. Try 'cholmod' when the relationship precision is dense or the LDLT factor fills in.", c("henderson", "deterministic")),
  tip("Solver", "Use solve$method = 'pcg' only when factor memory or factorization time is the bottleneck and approximate log-determinants and traces are acceptable. Check sensitivity to solve$tol, score$traceProbes and logdet$lanczosSteps.", "pcg"),
  tip("Solver", "With PCG, leave update$blendWeight at NULL or set it to 0: positive EM/PXEM blending needs a deterministic factorization.", c("pcg", "optimizer")),
  tip("Formulation", "Use formulation = 'direct' when records are few relative to random-effect levels; keep Henderson otherwise. Direct accepts only the short list on its node.", c("domain", "direct")),
  tip("Outputs", "Keep output$inverse = 0 unless you need coefficient uncertainty from the inverse. Mode 1 (selected inverse) suits LDLT; mode 2 (full inverse) can be very expensive in time and memory.", "outputs"),
  tip("Model size", "Unstructured covariance factors make the AI solve grow quadratically in the number of free parameters even when the MME solve is fast. Prefer FA or diagonal structures where they fit the biology.", c("model", "curvature")),
  tip("Model size", "Irregular residual patterns and missing cells remove complete-block shortcuts; the missingness pattern matters as much as the record count.", "model"),
  tip("Convergence", "If a fit oscillates or stalls, try update$recovery = 'pxem' or a blending schedule such as update$blendWeight = c(1, 0.5, 0) before raising convergence$maxIters.", c("optimizer", "convergence")),
  tip("Experimental", "'auto', 'off' and 'force' toggles record a request. 'force' never bypasses mathematical, structural or memory checks. Read fit$executionPlan and fit$engineDiagnostics for what actually ran.", c("deterministic", "representation", "curvature", "policy")),
  tip("Benchmarking", "Hold the model, tolerances, requested outputs, thread settings and iterations constant, and separate data preparation from native fitting time.", "policy"),
  tip("Threads and GPU", "execution$blasThreads sets BLAS threads per fit ('auto' needs flexiblas or RhpcBLASctl). GPU kernels help only large dense blocks; transfers can outweigh the gains.", "policy"))

# ------------------------------------------------------------------- render
data <- list(
  meta=list(version=read.dcf("DESCRIPTION", fields="Version")[1, 1],
    generated=format(Sys.Date()),
    defaults=paste(utils::capture.output(utils::str(unclass(defaults), give.attr=FALSE,
      no.list=TRUE, vec.len=3)), collapse="\n")),
  modules=capabilities$modules, controls=controls, nodes=nodes, edges=edges,
  rules=rules, fitChecks=fit_checks, structures=structures, tips=tips)
json <- jsonlite::toJSON(data, dataframe="rows", auto_unbox=TRUE, na="null", null="null")
json <- gsub("</", "<\\/", json, fixed=TRUE)
template <- readLines("pkgdown/mmes-options/template.html", warn=FALSE, encoding="UTF-8")
marker <- grep("__MMES_OPTIONS_JSON__", template, fixed=TRUE)
stopifnot(length(marker) == 1L)
template[marker] <- sub("__MMES_OPTIONS_JSON__", json, template[marker], fixed=TRUE)
outputs <- c("pkgdown/assets/mmes-options.html", commandArgs(trailingOnly=TRUE))
dir.create("pkgdown/assets", showWarnings=FALSE)
for(output in outputs) writeLines(template, output, useBytes=TRUE)
cat(sprintf("Wrote %s\n%d options | %d nodes | %d rules | %d fit-time checks | %d structure paths | %d tips\n",
  paste(outputs, collapse=", "), nrow(controls), length(nodes), length(rules), length(fit_checks),
  length(structures), length(tips)))
