# Functions needed for using jupyter lab

publish_gg <- function(gg) {
  "Display a ggplot"
  gg_bundle <- IRdisplay::prepare_mimebundle(gg)
  IRdisplay::publish_mimebundle(gg_bundle$data, gg_bundle$metadata)
}
