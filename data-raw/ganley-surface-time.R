# Generates `ganley_surface_time`.
#
# Transcribed from Ganley, L.C., Brault, S. and Mayo, C.A. (2019) "What we see
# is not what there is: estimating North Atlantic right whale Eubalaena
# glacialis local abundance", Endangered Species Research 38:101-113,
# doi:10.3354/esr00938. Open access under CC-BY.
#
# ONLY values the paper states exactly in its text are here. Specifically:
#
#   Section 3.2: "Percent surface time from focal follows are presented in
#   Fig. 2 (January: 16%; February: 34%; March: 31%; April: 55%)."
#
#   Fig. 2 caption: "Sample size for January = 7, February = 10, March = 22,
#   April = 48."
#
# Deliberately NOT transcribed: the monthly availability estimates. The paper
# reports January exactly (0.27) and gives the range in the abstract
# (0.27-0.85), but the remaining months appear only in Fig. 4A and in Table S1
# of the supplement, which is a separate file. Reading bar heights off a figure
# is not a measurement, and a column mixing exact and eyeballed values is worse
# than a column that stops where the text does.
#
# To add them: fetch the supplement at
# https://www.int-res.com/articles/suppl/n038p101_supp.pdf, take Table S1, and
# say in the documentation that that is where they came from.

ganley_surface_time <- tibble::tibble(
  month = factor(c("Jan", "Feb", "Mar", "Apr"),
                 levels = c("Jan", "Feb", "Mar", "Apr")),
  percent_surface_time = c(16, 34, 31, 55),
  n_follows = c(7L, 10L, 22L, 48L)
)

usethis::use_data(ganley_surface_time, overwrite = TRUE)
