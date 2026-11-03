#' @name gapfrac_fisheye
#' @aliases gapfrac_fisheye list.lenses
#' @title Derive angular gap fraction from a classified fisheye image
#' @importFrom terra rast values as.data.frame metags plot
#' @importFrom dplyr group_by nth right_join summarise mutate relocate ungroup select filter case_when bind_rows
#' @importFrom dismo circles
#' @importFrom tidyr pivot_wider
#' @importFrom graphics segments
#'
#' @param img.bw SpatLayer. A binary, single-layer fisheye image generated from [binarize_fisheye()].
#' @param maxVZA Numeric. The maximum Zenith angle (in degrees) corresponding to the image radius. Default= 90.
#' @param lens Character. The lens type for fisheye-lens correction. A list of lenses is available by typing *list.lenses*. If missing, it is assumed equidistant.
#' @param startVZA Numeric. The minimum Zenith angle (in degrees) considered in the analysis. Default is 0.
#' @param endVZA Numeric. The maximum Zenith angle (in degrees) considered in the analysis. Default is 70.
#' @param nrings Numeric. The number of equiangular zenith rings considered in the analysis. Default is 7.
#' @param nseg Numeric. The number of azimuth segments considered in the analysis. Default is 8.
#' @param display Logical. If set to TRUE, it desplays the zenith rings and azimuth segments overlaid on the fisheye image.
#' @param message Logical. If set to TRUE, it reports the circular mask used in the analysis.
#'
#' @return A dataframe of gap fraction (GF) for zenith rings (rows) and azimuth segments (columns).
#'
#' @description
#' The function calculates the gap fraction for a number of zenith annuli (rings) and azimuth sectors (segments).
#' A list of lens is available for correcting for lens distorsion. Type 'list.lenses'.
#'
#' @author
#' Francesco Chianucci
#'
#' @seealso
#' Lens correction functions have been retrieved from the following sources:
#'
#' Pekin and Macfarlane 2009: \doi{10.3390/rs1041298}
#'
#' Paul Bourke:<http://www.paulbourke.net/dome/fisheyecorrect/>
#'
#' Hemisfer: <https://www.schleppi.ch/patrick/hemisfer/index.php>
#'
#' @examples
#' \donttest{
#' c.im<-system.file('extdata/circular_coolpix4500+FC-E8_chestnut.jpg',package='hemispheR')
#'
#' #List of lenses for fisheye projection correction:
#' list.lenses
#'
#' #Zenith rings similar to LAI-2000/2200:
#' c.im |>
#'  import_fisheye(circ.mask=camera_fisheye('Coolpix4500+FC-E8')) |>
#'   binarize_fisheye() |>
#'   gapfrac_fisheye(lens='FC-E8',nrings=5,nseg=8,endVZA=75,display=TRUE)
#'
#' #The hinge angle method close to 1 radian (57):
#' c.im |>
#'  import_fisheye(circ.mask=camera_fisheye('Coolpix4500+FC-E8')) |>
#'   binarize_fisheye() |>
#'   gapfrac_fisheye(lens='FC-E8',nrings=1,nseg=8,startVZA=55,endVZA=60,display=TRUE)
#' }
#' @export
gapfrac_fisheye <- function(img.bw,
                            maxVZA=90,
                            lens='equidistant',
                            startVZA=0,
                            endVZA=70,
                            nrings=7,
                            nseg=8,
                            message=FALSE,
                            display=FALSE){

  r <- ring <- theta <- alpha.to <- value <- alpha.from <- id <- NULL

  if (length(setdiff(unique(terra::values(img.bw)),
                     c(NA,0, 1))) > 0){
    stop("Error: please select a binary (0,1) fisheye image, which can be classified using the binarize_fisheye() function")
  }
  if (length(setdiff(c(0,1),
                     unique(terra::values(img.bw))))>0){
    stop("Error: only 0 or 1 is contained in the binary image, calculation of gap fraction is not possible")
  }




  # p<-terra::as.points(img.bw, values=TRUE)
  # imgdf <- cbind(terra::crds(p), terra::values(p))
  imgdf <- terra::as.data.frame(img.bw, xy = TRUE)
  base::names(imgdf) <- c('x', 'y', 'value')

  #calculate back xc,yc,rc:
  xc=round(mean(imgdf$x))
  yc=round(mean(imgdf$y))

  if(length(unique(terra::values(img.bw))) == 2){
    rc = round(sqrt(xc^2 + yc^2)) - 2
    if(message == TRUE){
      message(paste('Used parameters xc, yc and rc are', xc, yc, rc, sep = ", "))
    }}

  if(length(unique(terra::values(img.bw))) == 3){
    rc=round((max(imgdf$x) - min(imgdf$x)) / 2)
    if(message == TRUE){
      message(paste('Used parameters xc, yc and rc are', xc, yc, rc, sep = ", "))
    }}

  #convert cartesian to polar coordinates
  #origin:
  imgdf$dx<-imgdf$x-mean(imgdf$x)
  imgdf$dy<-imgdf$y-mean(imgdf$y)
  imgdf$r<-round(sqrt(imgdf$dx^2+imgdf$dy^2))

  #azimuth theta:
  imgdf$theta <- NA
  # imgdf$theta <- atan2(imgdf$dy,imgdf$dx) #East-conterclockwise
  imgdf$theta <- base::atan2(imgdf$dx, imgdf$dy) #North-clockwise count
  imgdf$theta[imgdf$theta<0] <- imgdf$theta[imgdf$theta<0] + 2 * pi #rescale range to 0-2 pi instead of -pi, pi
  imgdf$theta <- imgdf$theta * 180 / pi #convert to degree


  startVZA=startVZA
  endVZA=endVZA
  maxVZA=maxVZA

  # no. zenith rings
  nrings = nrings
  calibration_file <- system.file(
    "data", "lens_data",
    package = "hemispheR"
  )

  lens_table <- read.csv(calibration_file)
  list.lenses <- lens_table$lens
  
  if (length(setdiff(lens, list.lenses)) > 0 | is.null(lens)){
    stop(cat('Unknown lens type. Allowed lenses are:\n"', sep = '', paste0(list.lenses, collapse = '",\n"'), '".'))
  }
  lens_params <- lens_table[lens_table$lens == lens, ]
  if (nrow(lens_params) == 0) {
    stop("Unknown lens type.")
  }
  
  VZAbins <- seq(startVZA, endVZA, length.out = nrings + 1)
  
  x <- if (lens_params$scale == "radians") {
    VZAbins * pi / 180
  } else {
    VZAbins / maxVZA
  }
  
  coef <- as.vector(lens_params[paste0("a", 1:5)])
  
  rset <- round(
    rc * rowSums(
      sweep(
        outer(x, seq_along(coef), "^"),
        2, coef, "*"
      )
     )
  )
  rset <- switch(
    lens_params$model,
    polynomial = round(rc * rowSums(
      sweep(outer(x, seq_along(coef), "^"), 2, coef, "*")
    )),
    equisolid = round(rc * sin(VZAbins * pi / 360) /
                      sin(maxVZA * pi / 360)),
    stereographic = round(rc * tan(VZAbins * pi / 360) /
                          tan(maxVZA * pi / 360)),
    orthographic = round(rc * sin(VZAbins * pi / 180) /
                         sin(maxVZA * pi / 180))
  )
  
  #8 azimuth segments
  nseg = nseg

  #define the sector of the image based on the rdf:
  options(dplyr.summarise.inform=FALSE)

  rdf <- imgdf |>
    dplyr::mutate(ring=cut(r,rset,include.lowest=TRUE,labels=seq(startVZA+((endVZA-startVZA)/2/nrings),endVZA-((endVZA-startVZA)/2/nrings),(endVZA-startVZA)/nrings)))|>
    # dplyr::mutate(ring=cut(r,rset,include.lowest=TRUE,labels=seq(setVZA/nrings,setVZA,setVZA/nrings)))|>
    tidyr::drop_na(ring) |>
    dplyr::mutate(ring=as.numeric(as.character(ring))) |>
    dplyr::mutate(alpha.to=cut(theta,seq(0,2*pi,2*pi/nseg)*180/pi,include.lowest=TRUE,labels=rep(seq(2*pi/nseg,2*pi,2*pi/nseg))*180/pi)) |>
    dplyr::mutate(alpha.to=trunc(as.numeric(as.character(alpha.to)))) |>
    dplyr::group_by(ring,alpha.to) |>
    dplyr::summarise(GF=mean(value,na.rm=TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(alpha.from=alpha.to-floor(as.numeric(as.character(((2*pi/nseg)*180/pi)))))|>
    # dplyr::mutate(GF=dplyr::case_when(
    #   GF==0~dplyr::nth(sort(unique(GF)),2),
    #   # GF==0~0.00004530,
    #   TRUE~GF
    # )) |>
    dplyr::relocate(ring,alpha.from)

  rdfw <- rdf |>
    tidyr::pivot_wider(id_cols=ring,names_from=c('alpha.from','alpha.to'),values_from='GF',names_prefix='GF') |>
    dplyr::mutate(id=base::names(img.bw),lens=lens) |>
    dplyr::mutate(circ=paste(xc,yc,rc, sep='_')) |>
    dplyr::relocate(id) |>
    dplyr::ungroup()


  # Create your default data frame
  metadf <- data.frame(channel = NA,
                       stretch = NA,
                       gamma = NA,
                       zonal = NA,
                       thd = NA,
                       method = NA,
                       stringsAsFactors = FALSE)

  tags <- tryCatch(
    terra::metags(img.bw),
    error = function(e) {
      message("Unable to retrieve metadata using terra::metags(). This may be due to a change in the terra package or the image format.")
      return(NULL)
    }
  )

  if (is.null(tags)) {
    message("Metadata is not available or could not be read from the image. Please ensure you are using a supported version of the 'terra' package or a valid image file.")
    metadata <- data.frame(name = c('channel',
                                    'stretch',
                                    'gamma',
                                    'zonal',
                                    'thd',
                                    'method'),
                           value=NA)
  } else {
    metadata <- terra::metags(img.bw)
  }




  # Fill in values only if the name is a valid column in metadf
  for (w in seq_len(nrow(metadata))) {
    if (metadata$name[w] %in% names(metadf)) {
      metadf[[metadata$name[w]]] <- metadata$value[w]
    }
  }


  rdfw <- cbind(rdfw,metadf)

  if (display == TRUE){
    terra::plot(img.bw,col=c('black','white'), main='applied rings & segments',legend=FALSE)
    #TODO
    for (j in 1:length(rset)){
      dd <- dismo::circles(cbind(xc,yc), rset[j], lonlat=F)
      plot(dd,  add=TRUE, border='yellow', lwd=3)
    }
    # slices<-seq(pi/nseg,2*pi,2*pi/nseg)
    slices <-seq(0,2*pi,2*pi/nseg)
    for (i in 1:length(slices)){
      graphics::segments(xc,yc,xc+max(rset)*sin(slices[i]),yc+max(rset)*cos(slices[i]),col='yellow',lwd=3)
    }
  }
  return(rdfw)
}
