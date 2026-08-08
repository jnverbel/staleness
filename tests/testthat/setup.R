# Every test that draws does so to a null device. Without this, any plot()
# call reached outside an explicitly opened device writes an Rplots.pdf into
# the test directory -- which then gets committed, shows up as a spurious diff
# on every run, and is exactly the sort of stray artefact a CRAN reviewer
# notices. Fixing it once here beats remembering a device in each test.
grDevices::pdf(NULL)
withr::defer(
  {
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
  },
  teardown_env()
)
