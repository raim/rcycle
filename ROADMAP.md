# rcycle roadmap: division-cycle pseudophase in time

Agreed with the owner on 2026-10-30. Goal: estimate how long cells
spend in each phase of the division cycle, and how variable that is,
from one snapshot of single-cell transcriptomes. This builds on what
`get_pseudophase` already returns (the PC1/PC2 angle `theta`, the rank
phase `phi`, cohort phases) and on the derivative and root machinery
in `R/segments.R` (points of maximal change, used to align phases).

## Why

* With `scale=TRUE` every cell is scaled to unit variance across
  cohorts, so all cells sit on a circle: amplitude is lost and noise
  cells still get an angle.
* The rank phase `phi` spreads cells evenly over the circle, so it
  hides where cells spend longer. `calibrate()` turns `phi` into time
  linearly, which ignores that an exponentially growing population has
  more young than old cells.
* For the Jackson et al. 2020 media, the stored per-cell classes
  (`cdc_class`, cohort of maximal fold change) give a G1 share that
  does not grow with the doubling time (Spearman -0.45 over 9 growing
  media).

## Ergodic rate analysis (ERA)

Kafri et al. 2013 (Box 2): in an asynchronous population growing
exponentially at rate `a = log(2)/T` (doubling time `T`), the speed
along the cycle path `l` is `v(l) = a (2 - F(l))/f(l)`, with `f` the
density of cells along `l` and `F` its cumulative fraction from
division. Integrating `dt = dl/v` gives the time since division

    age = T * (1 - log2(2 - F))

so `F = 0` is age 0 and `F = 1` is age `T`. Only the order of the cells
enters, so `theta` and `phi` give the same age; the density is in `F`.
The same balance for a sub-population that receives a share `s` of the
newborn cells (both cells of each division count, `s = 1/2` for
daughters or mothers in budding yeast), and that holds a fraction `f`
of all cells before it leaves a phase, gives its time in that phase

    tau = -T * log2(1 - f/(2*s))

which is the age formula for `s = 1`.

Assumptions: steady state (balanced, asynchronous exponential growth),
a single path that every cell follows in the same direction, and
division as the origin. Budding yeast divides asymmetrically: mothers
and daughters both restart in G1, daughters stay longer; the age above
is then the population mean, and the group formula splits it.

Pilot (2026-10-30, stored `cdc_phi`, origin at the M/G1 cohort peak,
G1 ends at the S cohort peak): the ERA G1 share rises with doubling
time, Spearman +0.75 over the 9 growing media (YPD 0.22 of T, Proline
0.43), against -0.45 for the `cdc_class` shares.

## Plan

All four parts were approved. Parts 1-3 and the principal curve of
part 4 are implemented; count-aware fitting is open.

1. **Cell age** (`R/phasetime.R`).
   `era_age()` and `era_fraction()` (the formula and its inverse);
   `phase_time()` adds the age of each cell and cohort to a phases
   object, origin at a phase or a cohort peak (`state_phase`), cells
   counted in the direction of increasing phase; `phase_landmarks()`
   gives cohort peaks and points of maximal change (maximum of the
   smoothed `d(theta-phi)/dphi` between two peaks, from
   `get_segments`); `phase_durations()` turns landmarks into ages and
   durations. G1/S is reported both at the S peak and at the point of
   maximal change between the G1 and S peaks. Validation: G1 duration
   against doubling time over the 11 Jackson et al. 2020 media.
2. **Per-cell quality** (`R/quality.R`).
   `project_cells()` puts new cohort state columns into an existing
   PCA (exact for the cells it was built from); `cell_quality()` adds
   the unscaled amplitude (`amplitude * phases$scale`), the share of a
   cell's variance across cohorts in the PC1/PC2 plane (`r2`), and
   flags cells that are not above a Poisson noise reference simulated
   from their own totals (`noise_states()`: cohort counts drawn from
   the population mean fractions, without a cycle); `phase_bootstrap()`
   gives a per-cell phase uncertainty by resampling genes within
   cohorts.
3. **Variability of durations** (`R/phasetime.R`).
   `group_durations()` splits a phase by a cell grouping (G1 by the
   `daughter` and `aged` cohort states) with the group formula above;
   `age_profile()` bins any per-cell value along age (mean, sd, cv),
   for the spread of log total counts as a size proxy, which should
   narrow towards G1 exit if size is controlled there; cells off the
   trajectory (`r2`, radius, `pphi`) against the noise reference.
4. **Geometry**.
   `principal_phase()` (`R/principal.R`) fits a closed principal curve
   (Fourier series in each coordinate, Hastie-Stuetzle iterations)
   through the cells in the first `k` PCs, so more than two PCs, and
   so more cohorts than five cohort means, can define the path; the
   arc-length phase `lambda` and the distance to the curve are
   returned, and `phase_time()` takes `phase="lambda"`.
   Open: count-aware (Poisson) phase fitting with uncertainty, i.e.
   maximum likelihood of each cell's counts along the fitted cohort
   profiles, which would replace both the angle and the bootstrap.

Example: `vignettes/phase_time.Rmd` (Jackson et al. 2020, YPD, and the
11-media check).

## First results (2026-10-30, `inst/scripts/jackson20_phasetime.R`)

* G1 (M/G1 peak to G1S) as a share of T rises with the doubling time:
  Spearman +0.78 over the 9 growing media, +0.90 without noise cells
  (+0.83 for the 6 media in cycle order); to the S peak +0.75 / +0.77.
  YPD: G1 0.17 h of T = 1.05 h.
* Size proxy (YPD): mean log total is lowest from G1 to S; its spread
  is narrowest around the G2 peak, not at G1 exit.
* **The cohort order is wrong in several media.** From the M/G1 peak,
  Proline, Glutamine and MinimalEtOH give M, S, G1, G2 (M right after
  M/G1, G1 after S), so their G1 includes M; YPDDiauxic and CStarve
  are out of order too. YPEtOH, YPDRapa, YPD, Urea, MinimalGlucose and
  AmmoniumSulfate are in cycle order. The principal curve (k = 3 or 4) repairs YPDDiauxic but
  breaks YPDRapa and AmmoniumSulfate, so it is no fix. Next: test the
  cohort order before reading durations (the validation TODO in
  `R/phases.R`), and check whether M and M/G1 should be merged.
* Noise cells (amplitude not above Poisson cells of the same totals):
  13-33% in growing media, 58% in YPDDiauxic, 70% in CStarve.
* Per-cell uncertainty from gene resampling (YPD): median 0.17 of the
  cycle (0.19 h). Single cells are placed only roughly; durations come
  from the population order.
* The `daughter` cohort is high in ~90% of cells right after division
  and again in late M: it marks M/G1, not daughters, so mother and
  daughter G1 (`group_durations`) need another marker, or the daughter
  counts relative to age-matched cells.

## Not here

The Heineike data sets in the pseudophase cache are unpublished and
stay out of public outputs, including this package's vignettes.
