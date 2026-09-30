## arrows.phases (2026-09-30): cohort phases, columns of x$x.phase, were
## looked up in x$x.phases and always omitted; the omission warning pasted
## the missing types with sep= instead of collapse=

library(testthat)
library(rcycle)

## a phases object with one segment table and cohort phases
ph <- seq(-pi, pi, length.out = 61)[-61]
st <- t(sapply(0:5, function(i) 2 + cos(ph - i*pi/3)))
rownames(st) <- paste0("c", 1:6)
phases <- get_pseudophase(st, verb = 0)
phases$seg <- data.frame(ID = 1:3, phi = c(-2, 0, 2), amp = 3:1,
                         class = c("a", "b", "c"))

test_that("arrows.phases draws segments and cohort phases", {
    pdf(NULL); on.exit(dev.off())
    plot(1, xlim = c(-pi, pi), ylim = c(0, 1))
    expect_no_warning(arrows.phases(phases, types = c("seg", "phi"), y0 = .5, dy = .1))
    expect_no_warning(arrows.phases(phases, types = "phi", labels.top = 2))
    expect_no_warning(arrows.phases(phases, types = "seg", labels = "class",
                                    labels.top = 1, lxpd = TRUE))
})

test_that("arrows.phases names all missing types in one warning", {
    pdf(NULL); on.exit(dev.off())
    plot(1, xlim = c(-pi, pi), ylim = c(0, 1))
    expect_warning(arrows.phases(phases, types = c("seg", "nope", "none")),
                   "nope;none")
})
