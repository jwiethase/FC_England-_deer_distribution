get_pal <- function(n=50, type=c("turbo", "viridis", "inferno", "spectral", "bam")) {
      switch(match.arg(type),
             "turbo" = viridis::viridis_pal(option = "H")(n),
             "viridis" = viridis::viridis_pal(option = "D")(n),
             "spectral" = grDevices::hcl.colors(n, "spectral", rev = F),
             "inferno" = viridis::viridis_pal(option = "B")(n),
             "bam" = grDevices::colorRampPalette(
                   c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", 
                     "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", 
                     "#007FFF", "blue"), space="Lab", bias=0.5)(n))
}
base <- function(type) {
      map_attr = "© <a href='https://www.esri.com/en-us/home'>ESRI</a> © <a href='https://www.google.com/maps/'>Google</a>"
      leaflet(options = leafletOptions(zoomControl = FALSE))  %>% 
            addTiles(urlTemplate = "http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z}&s=Ga",
                     group  =  "Google")  %>% 
            addProviderTiles("CartoDB.Positron", group = "CartoDB")  %>% 
            addProviderTiles("OpenStreetMap", group = "Open Street Map")  %>% 
            addProviderTiles('Esri.WorldImagery', group = "ESRI")  %>% 
            addTiles(urlTemplate = "", attribution = map_attr)  %>% 
            addLayersControl(
                  baseGroups = c("Open Street Map", "ESRI", "CartoDB", "Google"),
                  options = layersControlOptions(collapsed = FALSE)) %>% 
            setView(-1.264, 52.531, 6) %>%
            htmlwidgets::onRender("function(el, x) {
        L.control.zoom({ position: 'bottomright' }).addTo(this)
    }")
}

addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft","topleft"),
                                  pal, values, na.label = "NA", bins = 7, colors, 
                                  opacity = 0.5, labels = NULL, labFormat = labelFormat(), 
                                  title = NULL, className = "info legend", layerId = NULL, 
                                  group = NULL, data = getMapData(map), decreasing = FALSE) {
      
      position <- match.arg(position)
      type <- "unknown"
      na.color <- NULL
      extra <- NULL
      if (!missing(pal)) {
            if (!missing(colors)) 
                  stop("You must provide either 'pal' or 'colors' (not both)")
            if (missing(title) && inherits(values, "formula")) 
                  title <- deparse(values[[2]])
            values <- evalFormula(values, data)
            type <- attr(pal, "colorType", exact = TRUE)
            args <- attr(pal, "colorArgs", exact = TRUE)
            na.color <- args$na.color
            if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] == 
                0) {
                  na.color <- NULL
            }
            if (type != "numeric" && !missing(bins)) 
                  warning("'bins' is ignored because the palette type is not numeric")
            if (type == "numeric") {
                  cuts <- if (length(bins) == 1) 
                        pretty(values, bins)
                  else bins   
                  if (length(bins) > 2) 
                        if (!all(abs(diff(bins, differences = 2)) <= 
                                 sqrt(.Machine$double.eps))) 
                              stop("The vector of breaks 'bins' must be equally spaced")
                  n <- length(cuts)
                  r <- range(values, na.rm = TRUE)
                  cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
                  n <- length(cuts)
                  p <- (cuts - r[1])/(r[2] - r[1])
                  extra <- list(p_1 = p[1], p_n = p[n])
                  p <- c("", paste0(100 * p, "%"), "")
                  if (decreasing == TRUE){
                        colors <- pal(rev(c(r[1], cuts, r[2])))
                        labels <- rev(labFormat(type = "numeric", cuts))
                  }else{
                        colors <- pal(c(r[1], cuts, r[2]))
                        labels <- rev(labFormat(type = "numeric", cuts))
                  }
                  colors <- paste(colors, p, sep = " ", collapse = ", ")
            }
            else if (type == "bin") {
                  cuts <- args$bins
                  n <- length(cuts)
                  mids <- (cuts[-1] + cuts[-n])/2
                  if (decreasing == TRUE){
                        colors <- pal(rev(mids))
                        labels <- rev(labFormat(type = "bin", cuts))
                  }else{
                        colors <- pal(mids)
                        labels <- labFormat(type = "bin", cuts)
                  }
            }
            else if (type == "quantile") {
                  p <- args$probs
                  n <- length(p)
                  cuts <- quantile(values, probs = p, na.rm = TRUE)
                  mids <- quantile(values, probs = (p[-1] + p[-n])/2, na.rm = TRUE)
                  if (decreasing == TRUE){
                        colors <- pal(rev(mids))
                        labels <- rev(labFormat(type = "quantile", cuts, p))
                  }else{
                        colors <- pal(mids)
                        labels <- labFormat(type = "quantile", cuts, p)
                  }
            }
            else if (type == "factor") {
                  v <- sort(unique(na.omit(values)))
                  colors <- pal(v)
                  labels <- labFormat(type = "factor", v)
                  if (decreasing == TRUE){
                        colors <- pal(rev(v))
                        labels <- rev(labFormat(type = "factor", v))
                  }else{
                        colors <- pal(v)
                        labels <- labFormat(type = "factor", v)
                  }
            }
            else stop("Palette function not supported")
            if (!any(is.na(values))) 
                  na.color <- NULL
      }
      else {
            if (length(colors) != length(labels)) 
                  stop("'colors' and 'labels' must be of the same length")
      }
      legend <- list(colors = I(unname(colors)), labels = I(unname(labels)), 
                     na_color = na.color, na_label = na.label, opacity = opacity, 
                     position = position, type = type, title = title, extra = extra, 
                     layerId = layerId, className = className, group = group)
      invokeMethod(map, data, "addLegend", legend)
}

range01 <- function(x){(x - min(x))/(max(x) - min(x))}
