# Sole editable source for solvephrust's built-in catalog.
# Run from the package root: Rscript data-raw/compile_data.R
# Dissociation constants and reaction enthalpies are retained from the original
# implementation, based on Benjamin (2015), Appendices A.1 and A.2.

.solvephrust_default_compounds <- list(
  ph_dependent = list(
    co3 = list(
      constants = list(
        list(k = 10^-6.35, delta_h = 7700),
        list(k = 10^-10.33, delta_h = 14900)
      ), charge = -1L
    ),
    po4 = list(
      constants = list(
        list(k = 10^-2.16, delta_h = -8000),
        list(k = 10^-7.20, delta_h = 4200),
        list(k = 10^-12.35, delta_h = 14700)
      ), charge = -1L
    ),
    ocl = list(
      constants = list(list(k = 10^-7.53, delta_h = 13800)), charge = -1L
    ),
    nh3 = list(
      constants = list(list(k = 10^-9.244, delta_h = 52210)), charge = 1L
    ),
    ch3coo = list(
      constants = list(list(k = 10^-4.757, delta_h = -200)), charge = -1L
    ),
    bo3 = list(
      constants = list(list(k = 10^-9.24, delta_h = -42000)), charge = -1L
    ),
    sio4 = list(
      constants = list(
        list(k = 10^-9.84, delta_h = 25600),
        list(k = 10^-13.2, delta_h = 37000)
      ), charge = -1L
    )
  ),
  ph_independent_charges = list(
    so4 = -2L,
    na = 1L,
    ca = 2L,
    mg = 2L,
    cl = -1L,
    mno4 = -1L,
    no3 = -1L
  )
)

dir.create("R", showWarnings = FALSE)
save(.solvephrust_default_compounds, file = file.path("R", "sysdata.rda"),
     version = 3, compress = "xz")
