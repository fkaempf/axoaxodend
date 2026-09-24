# make_mcns_scene(): base Clio-NG scene for the male CNS (EM, segmentation,
# brain and VNC shells). Copied from pheromone.paper/R/neuroglancer.R.

make_mcns_scene <- function(
  ids,
  title        = "CNS-Quad",
  em_url       = "precomputed://gs://cns-full-clahe",
  malecns_url  = "dvid://https://emdata6-novran.janelia.org/4b2087c0fbe046bfaf0d60bc970e3e5d/segmentation?dvid-service=https://ngsupport-bmcp5imp6q-uk.a.run.app",
  brain_shell_url =
    "precomputed://gs://flyem-cns-roi-7c971aa681da83f9a074a1f0e8ef60f4/brain-shell-smooth-linear",
  vnc_shell_url   =
    "precomputed://gs://flyem-cns-roi-7c971aa681da83f9a074a1f0e8ef60f4/vnc-shell",
  pl_url   = NULL,
  pl_name  = "plshere"
) {
  ## keep exact input order
  ids_chr <- as.character(ids)

  ## EM base layer
  em_layer <- list(
    type   = "image",
    source = list(
      url                     = em_url,
      subsources              = list(default = TRUE),
      enableDefaultSubsources = FALSE
    ),
    tab   = "rendering",
    name  = "em"
  )

  ## malecns segmentation layer
  malecns_layer <- list(
    type   = "segmentation",
    source = list(
      list(
        url                     = malecns_url,
        subsources              = list(default = TRUE, meshes = TRUE),
        enableDefaultSubsources = FALSE
      ),
      "precomputed://https://ngsupport-bmcp5imp6q-uk.a.run.app/neuronjson_segment_properties/emdata6-novran.janelia.org/4b2087c0fbe046bfaf0d60bc970e3e5d/segmentation_annotations/type/group"
    ),
    toolBindings = list(Q = "selectSegments"),
    tab          = "segments",
    segments     = ids_chr,                         # order preserved
    segmentQuery = paste(ids_chr, collapse = " "),  # order preserved
    name         = "malecns"
  )

  ## brain-shell layer
  brain_shell_layer <- list(
    type   = "segmentation",
    source = list(
      url                     = brain_shell_url,
      subsources              = list(default = TRUE, properties = TRUE, mesh = TRUE),
      enableDefaultSubsources = FALSE
    ),
    pick                    = FALSE,
    tab                     = "rendering",
    selectedAlpha           = 0,
    saturation              = 0,
    meshSilhouetteRendering = 7,
    segments                = list("1"),
    colorSeed               = 1336242844,
    segmentDefaultColor     = "#ffffff",
    name                    = "brain-shell"
  )

  ## vnc-shell layer
  vnc_shell_layer <- list(
    type   = "segmentation",
    source = list(
      url                     = vnc_shell_url,
      subsources              = list(default = TRUE, properties = TRUE, mesh = TRUE),
      enableDefaultSubsources = FALSE
    ),
    pick                    = FALSE,
    tab                     = "segments",
    selectedAlpha           = 0,
    saturation              = 0,
    meshSilhouetteRendering = 7,
    segments                = list("1"),
    colorSeed               = 1336242844,
    segmentDefaultColor     = "#ffffff",
    name                    = "vnc-shell"
  )

  layers <- list(
    em_layer,
    malecns_layer,
    brain_shell_layer,
    vnc_shell_layer
  )

  ## optional pl reference layer
  if (!is.null(pl_url)) {
    pl_layer <- list(
      type   = "segmentation",
      source = list(
        url                     = pl_url,
        subsources              = list(default = TRUE, properties = TRUE, mesh = TRUE),
        enableDefaultSubsources = FALSE
      ),
      pick     = FALSE,
      tab      = "segments",
      segments = list("1"),
      name     = pl_name
    )
    layers <- c(layers, list(pl_layer))
  }

  scene <- list(
    title      = title,
    dimensions = list(
      x = list(8e-9, "m"),
      y = list(8e-9, "m"),
      z = list(8e-9, "m")
    ),
    ## orientation + zoom taken from your reference URL
    position             = c(47470.59375, 54212.50390625, 66301.59375),
    crossSectionScale    = 1.1208,
    projectionOrientation= c(
      0,
      0.7071067690849304,
      0.7071067690849304,
      0
    ),
    projectionScale      = 192860.74860894235,
    layers               = layers,
    showAxisLines        = FALSE,
    showSlices           = FALSE,
    prefetch             = FALSE,
    selectedLayer        = list(
      visible = TRUE,
      layer   = "malecns"
    ),
    layout               = "3d"
  )

  scene
}
