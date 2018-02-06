## * .onLoad
.onLoad <- function(lib, pkg="lavaSearch2") {

    # available methods to compute the distribution of the max statistic
    lava::lava.options(search.calcMaxDist = c("integration","boot-naive","boot-residual","boot-wild"),
                       search.statistic = c("Wald","score","LR"),
                       search.p.adjust = c("fastmax", "max", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                       search.calc.quantile.int = TRUE
                       )
}

## * .onAttach
.onAttach <- function(lib, pkg="lavaSearch2") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

lava_categorical2dummy <- get("categorical2dummy", envir = asNamespace("lava"), inherits = FALSE)
lava_estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE)
