# rcycle

A fast and scalable method to determine a `pseudophase` of cyclic transcription processes from single cell transcriptome (snap shot) data.

![https://vecteezy.com](vignettes/elephant.png)

TODO:

* clearer separate from prcomp, 
    - use get_pseudophase to "decorate" a prcomp result, but perhaps store       in additional tables rather than rotation and x,
* generalize plotPC: 
    - it currently works for prcomp package: keep it that way!
    -
* PWM MODEL: keep here or move to separate package?
* PWM MODEL `kb` (from ~/work/teaching/pwm/ROADMAP.md, 2026-10-02): pulsed
  transcription with basal transcription in the OFF phase and constant
  degradation, `dR/dt = kappa*k + (1-kappa)*k0 - (dr+mu)*R`. Introduced by
  Claude in the D6 analyses (2026-09-30); so far only its closed-form mean,
  in ~/work/teaching/pwm/d6_kdrk0_fit.R. To add: `pwmode_kb`, and `kb` in
  `get_rcycle`, `get_rmean`, `get_ramp`, `get_rates`, `get_times` (+ tests).
  `R - k0/gamma` follows model `k` with ON rate `k - k0` (reversed if
  `k0 > k`), so the mean is `(phi*k + (1-phi)*k0)/gamma`, exact with growth
  since the loss rate is the same in both phases.
