library(ResTree)
library(GpGp)
library(DiceDesign)

source(file.path(script.dir, "api_contract.R"))

method.labels = c(
  restgp_full = "ResTGP(Full)", restgp_pp = "ResTGP(PP)", ssgv = "sSGV"
)
knob.names = c(restgp_full = "r", restgp_pp = "r", ssgv = "m")

source(file.path(script.dir, "helpers", "timing.R"))
source(file.path(script.dir, "helpers", "data.R"))
source(file.path(script.dir, "helpers", "exact_gp.R"))
source(file.path(script.dir, "helpers", "methods.R"))
source(file.path(script.dir, "helpers", "summaries.R"))
source(file.path(script.dir, "helpers", "run_case.R"))
