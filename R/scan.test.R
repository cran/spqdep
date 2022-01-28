#'
#' @title Compute the scan test.
#'
#' @description This function compute the scan test for Bernoulli and Multinomial categorical spatial process.
#'
#' @param data an (optional) data frame or a sf object containing the variable to testing for.
#' @param formula a symbolic description of the factor (optional).
#' @param fx a factor (optional).
#' @param coor (optional) coordinates of observations.
#' @param case Only for Bernoulli distribution. A element of factor, there are cases and non-cases for testing for cases versus non-cases
#' @param nv Maximum windows size, default nv = N/2. The algorithm scan for clusters of geographic size between 1
#' and the upper limit (nv) defined by the user.
#' @param nsim Number of permutations.
#' @param alternative Only for Bernoulli spatial process. A character string specifying the type of cluster, must be one
#' of "High" (default), "Both" or "Low".
#' @param distr distribution of the spatial process: "bernoulli" for two levels or "multinomial" for three or more levels.
#' @param windows a string to select the type of cluster "circular" (default) of "elliptic".
#' @param control List of additional control arguments.
#' @usage scan.test(formula = NULL, data = NULL, fx = NULL, coor = NULL, case = NULL, nv = NULL,
#' nsim = NULL, distr = NULL, windows = "circular", alternative = "High", control = list())
#' @details
#'
#' Two alternative sets of arguments can be included in this function to compute the scan test:
#'
#'   \tabular{ll}{
#'     \code{Option 1} \tab A factor (fx) and coordinates (coor). \cr
#'     \code{Option 2} \tab A sf object (data) and the formula to specify the factor.
#'     The function consider the coordinates of the centroids of the elements of th sf object. \cr
#'     }
#'
#'  The spatial scan statistics are widely used in epidemiology, criminology or ecology.
#'  Their purpose is to analyze the spatial distribution of points or geographical
#'  regions by testing the hypothesis of spatial randomness distribution on the
#'  basis of different distributions (e.g. Bernoulli, Poisson or Normal distributions).
#'  The \code{scan.test} function obtain the scan statistic for the Bernoulli and Multinomial distribution.
#'
#'  To test independence in a spatial process, under the null, the type of windows is irrelevant but under the alternative the elliptic
#'  windows can to identify with more precision the cluster.
#'
#'  For big data sets (N >>) the windows = "elliptic" can be so slowly
#'
#'
#'
#' @return A object of the \emph{htest} and \emph{scantest} class
#'   \tabular{ll}{
#'     \code{method} \tab The type of test applied ().\cr
#'     \code{fx} \tab Factor included as input to get the scan test.\cr
#'     \code{MLC} \tab Observations included into the Most Likelihood Cluster (MLC).\cr
#'     \code{statistic} \tab Value of the scan test (maximum Log-likelihood ratio). \cr
#'     \code{N} \tab Total number of observations.\cr
#'     \code{nn} \tab Windows used to get the cluster.\cr
#'     \code{nv} \tab Maximum number of observations into the cluster.\cr
#'     \code{data.name} \tab A character string giving the name of the factor.\cr
#'     \code{coor} \tab coordinates.\cr
#'     \code{alternative} \tab Only for bernoulli spatial process. A character string describing the alternative hypothesis select by the user.\cr
#'     \code{p.value} \tab p-value of the scan test.\cr
#'     \code{cases.expect} \tab Expected cases into the MLC.\cr
#'     \code{cases.observ} \tab Observed cases into the MLC.\cr
#'     \code{nsim} \tab Number of permutations.\cr
#'     }
#' @section Control arguments:
#'   \tabular{ll}{
#'     \code{seedinit} \tab Numerical value for the seed (only for boot version). Default value seedinit=123 \cr
#'       }

#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Kulldorff M, Nagarwalla N. (1995).
#'     Spatial disease clusters: Detection and Inference.
#'       \emph{Statistics in Medicine}. 14:799-810
#'     \item Jung I, Kulldorff M, Richard OJ (2010).
#'     A spatial scan statistic for multinomial data.
#'       \emph{Statistics in Medicine}. 29(18), 1910-1918
#'      \item Páez, A., López-Hernández, F. A., Ortega-García, J. A., & Ruiz, M. (2016).
#'      Clustering and co-occurrence of cancer types: A comparison of techniques with an application to pediatric cancer in Murcia, Spain.
#'      \emph{Spatial Analysis in Health Geography}, 69-90.
#'
#'   }
#' @seealso
#' \code{\link{local.sp.runs.test}}, \code{\link{dgp.spq}}, \code{\link{Q.test}},
#' @export
#'
#' @examples
#'
#' # Case 1: Scan test bernoulli
#' data(provinces_spain)
#' sf::sf_use_s2(FALSE)
#' provinces_spain$Male2Female <- factor(provinces_spain$Male2Female > 100)
#' levels(provinces_spain$Male2Female) = c("men","woman")
#' formula <- ~ Male2Female
#' scan <- scan.test(formula = formula, data = provinces_spain, case="men",
#' nsim = 99, distr = "bernoulli")
#' print(scan)
#' summary(scan)
#' plot(scan, sf = provinces_spain)
#'
#' ## With maximum number of neighborhood
#' scan <- scan.test(formula = formula, data = provinces_spain, case = "woman", nv = 5,
#' nsim = 99, distr = "bernoulli")
#' print(scan)
#' plot(scan, sf = provinces_spain)
#'
#' \donttest{
#' ## With elliptic windows
#' scan <- scan.test(formula = formula, data = provinces_spain, case = "men", nv = 20,
#' nsim = 99, distr = "bernoulli", windows ="elliptic")
#' print(scan)
#' scan <- scan.test(formula = formula, data = provinces_spain, case = "men", nv = 15,
#' nsim = 99, distr = "bernoulli", windows ="elliptic", alternative = "Low")
#' print(scan)
#' plot(scan, sf = provinces_spain)
#'
#'
#' # Case 2: scan test multinomial
#' data(provinces_spain)
#' provinces_spain$Older <- cut(provinces_spain$Older, breaks = c(-Inf,19,22.5,Inf))
#' levels(provinces_spain$Older) = c("low","middle","high")
#' formula <- ~ Older
#' scan <- scan.test(formula = formula, data = provinces_spain, nsim = 99, distr = "multinomial")
#' print(scan)
#' plot(scan, sf = provinces_spain)
#'
#' # Case 3: scan test multinomial
#' data(FastFood.sf)
#' sf::sf_use_s2(FALSE)
#' formula <- ~ Type
#' scan <- scan.test(formula = formula, data = FastFood.sf, nsim = 99,
#' distr = "multinomial", windows="elliptic", nv = 50)
#' print(scan)
#' summary(scan)
#' plot(scan, sf = FastFood.sf)
#'
#' # Case 4: DGP two categories
#' N <- 500
#' cx <- runif(N)
#' cy <- runif(N)
#' listw <- spdep::knearneigh(cbind(cx,cy), k = 10)
#' p <- c(1/2,1/2)
#' rho <- 0.5
#' fx <- dgp.spq(p = p, listw = listw, rho = rho)
#' scan <- scan.test(fx = fx, nsim = 9, case = "A", nv = 30, coor = cbind(cx,cy),
#' distr = "bernoulli",windows="elliptic")
#' print(scan)
#' plot(scan)
#'
#' # Case 5: DGP three categories
#' N <- 200
#' cx <- runif(N)
#' cy <- runif(N)
#' listw <- spdep::knearneigh(cbind(cx,cy), k = 10)
#' p <- c(1/3,1/3,1/3)
#' rho <- 0.5
#' fx <- dgp.spq(p = p, listw = listw, rho = rho)
#' scan <- scan.test(fx = fx, nsim = 19, coor = cbind(cx,cy), nv = 30,
#' distr = "multinomial", windows = "elliptic")
#' print(scan)
#' plot(scan)
#' }


scan.test <- function(formula = NULL, data = NULL, fx = NULL, coor = NULL, case = NULL,
                      nv = NULL, nsim = NULL, distr = NULL, windows = "circular",
                      alternative = "High", control = list()) {

  if (is.null(distr))
    stop("Select a distribution, bernoulli or multinomial")

  coor.input <- coor
  distr <- match.arg(distr, c("bernoulli", "multinomial"))
  windows <- match.arg(windows, c("circular", "elliptic"))

  # Selecciona los argumentos. Bien con (formula + data) o bien incluye la variable (fx)
  if (!is.null(formula) && !is.null(data)) {
    if (inherits(data, "Spatial")) data <- as(data, "sf")
    mfx <- get_all_vars(formula, data)[,1]
    data.name <- names(get_all_vars(formula, data))
  } else if (!is.null(fx)) {
    mfx <- fx
    if (is.null(names(fx))) data.name <- "fx"
  } else stop("data wrong")

  if (!is.factor(mfx))
    stop(paste(deparse(substitute(fx)), "is not a factor"))

  # Si se introduce un objeto sf. Para controlar errores en el metodo plot
  sf = FALSE
  if (inherits(data, "sf")) sf = TRUE

  if (distr == "bernoulli"){
  alternative <- match.arg(alternative, c("Low", "High", "Both"))
  if (is.null(case)) {stop ("case argument must be an element of the factor")}
  if (length(unique(mfx))!=2) {stop ("The factor mut be have 2 levels for bernoulli")}
  case <- match.arg(case, unique(mfx))
  }
  if (distr == "multinomial"){
    if (length(unique(mfx)) < 3) {stop ("The factor must be have almost 3 levels for multinomial")}
    case <- match.arg(case, unique(mfx))
  }

  # Controls
  con <- list(seedinit = 123)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  seedinit <- con$seedinit

  ## Scan
  N <- length(mfx)
  if (is.null(nv))  nv <- trunc(N/2)
  if (!is.null(nv) && nv > trunc(N/2)) stop("nv must be lower than N/2")
  ## Previus

  if (is.null(coor) &&
      sum(inherits(data, "sf")) > 0){
  coor <- suppressWarnings(st_coordinates(st_centroid(data)))
  }

  cx <- coor[,1]
  cy <- coor[,2]

  if (windows=="circular"){
    nn <- suppressWarnings(cbind(1:N, knearneigh(cbind(cx, cy),
                                k = (nv-1))$nn))
  }
  if (windows=="elliptic"){
    nn <- suppressWarnings(nn_ellipse(coor = cbind(cx,cy), nv = nv, p = 30)$ellipses)
  }

  # XF <- matrix(mfx[nn], ncol = nv, nrow = N)
  XF <- matrix(mfx[nn], ncol = nv, nrow = dim(nn)[1])

  #####################################
  ## Obtaining the scan statisitic
  #####################################
  ## Bernoulli
  if (distr == "bernoulli"){
  oz <- t(apply(XF == case, 1 , cumsum))
  nz <- t(matrix(rep(1:nv,dim(nn)[1]),nrow = nv))
  O <- sum(mfx== case)
  oznz <- oz/nz
  OozNnz <- (O-oz)/(N-nz)
  OozNnz.1 <- 1 - OozNnz
  OozNnz.1[OozNnz.1 < 0] <- 0
  a <- log(oznz)
  a[a==-Inf]<-0
  b <- log(1-oznz)
  b[b==-Inf]<-0
  c <- log(OozNnz)
  c[c==-Inf]<-0
  d <- log(OozNnz.1)
  d[d==-Inf] <- 0
  if (alternative == "Both"){
  lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)
  }
  if (alternative == "High"){
    lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)*(oznz>OozNnz)
  }
  if (alternative == "Low"){
    lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)*(oznz<OozNnz)
  }
  # lnlz <- oz*log(oz/nz) + (nz-oz)*log(1-oz/nz) + (O-oz)*log((O-oz)/(N-nz)) + (N-O-nz+oz)*log(1-(O-oz)/(N-nz))
  lnlz0 <- exp(O*log(O) + (N-O)*log(N-O) - N*log(N))
  lnlz <- log(lnlz/lnlz0)
  lnlz[lnlz==-Inf] <- 0
  lnlz[is.na(lnlz)] <- 0

  a <- which(lnlz == max(lnlz), arr.ind = TRUE)
  if (dim(a)[1] > 1) {
  warning(paste0("A total of ",dim(a)[1], " clusters has the same value of the statistic. Report as MLC only the first"))
  a <- a[1,]
  }
  MLC <- nn[a[1],1:a[2]]
  loglik <- max(lnlz)
  cases.observ <- sum(mfx[MLC]== case)
  cases.expect <- a[2]*(O/N)
  }

  ## Multinomial
  if (distr == "multinomial"){
  lnlz <- 0
  case <- unique(mfx)
  CZ <- t(matrix(rep(1:nv,dim(nn)[1]),nrow = nv))
  for (f in 1:length(case)){
  Ck <- sum(mfx == case[f])
  CkZ <- t(apply(XF == case[f] , 1 , cumsum))
  a <- log(CkZ/CZ)
  a[a==-Inf] <- 0
  b <- log((Ck-CkZ)/(N-CZ))
  b[b==-Inf] <- 0
  c <- log(Ck/N)
  c[c==-Inf] <- 0
  lnlz <- lnlz + (CkZ*a + (Ck-CkZ)*b - Ck*c)
  }
  lnlz[is.na(lnlz)] <- 0
  a <- which(lnlz == max(lnlz), arr.ind = TRUE)
  if (dim(a)[1] > 1) {
    warning(paste0("A total of ",dim(a)[1], " clusters has the same value of the statistic. Report as MLC only the first"))
    a <- a[1,]
  }

  MLC <- nn[a[1],1:a[2]]
  loglik <- max(lnlz)
  cases.observ <- addmargins(table(mfx[MLC]))
  cases.expect <- addmargins(table(mfx)*length(MLC)/N)
  }


  #####################################
  ## Scan mc
  #####################################
  if (distr == "bernoulli"){
  if (!is.null(seedinit)) set.seed(seedinit)
  scan.mc <- rep(0,nsim)
  for (f in 1:nsim){
  fxp <- mfx[sample(N)]
  XF <- matrix(fxp[nn], ncol = nv, nrow = dim(nn)[1])
  oz <- t(apply(XF==case, 1 , cumsum))
  oznz <- oz/nz
  OozNnz <- (O-oz)/(N-nz)
  OozNnz.1 <- 1 - OozNnz
  OozNnz.1[OozNnz.1 < 0] <- 0
  a <- log(oznz)
  a[a==-Inf] <- 0
  b <- log(1-oznz)
  b[b==-Inf] <- 0
  c <- log(OozNnz)
  c[c==-Inf] <- 0
  d <- log(OozNnz.1)
  d[d==-Inf] <- 0
  if (alternative == "Both"){
    lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)
  }
  if (alternative == "High"){
    lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)*(oznz > OozNnz)
  }
  if (alternative == "Low"){
    lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)*(oznz < OozNnz)
  }
  lnlz <- log(lnlz/lnlz0)
  lnlz[lnlz==-Inf] <- 0
  lnlz[is.na(lnlz)] <- 0
  scan.mc[f] <- max(lnlz)
  }
  }

  if (distr == "multinomial"){
  if (!is.null(seedinit)) set.seed(seedinit)
  scan.mc <- rep(0,nsim)
  case <- unique(mfx)
  CZ <- t(matrix(rep(1:nv,dim(nn)[1]),nrow = nv))
  for (k in 1:nsim){
  fxp <- mfx[sample(N)]
  XF <- matrix(fxp[nn], ncol = nv, nrow = dim(nn)[1])
  lnlz <- 0
    for (f in 1:length(case)){
      Ck <- sum(mfx == case[f])
      CkZ <- t(apply(XF == case[f] , 1 , cumsum))
      a <- log(CkZ/CZ)
      a[a==-Inf] <- 0
      b <- log((Ck-CkZ)/(N-CZ))
      b[b==-Inf] <- 0
      c <- log(Ck/N)
      c[c==-Inf] <- 0
      lnlz <- lnlz + (CkZ*a + (Ck-CkZ)*b - Ck*c)
    }
    lnlz[is.na(lnlz)] <- 0
    scan.mc[k] <- max(lnlz)
  }
  }
  p.value <- sum(scan.mc > loglik)/(nsim + 1)
  names(loglik) <- "scan-loglik"

  if (distr == "bernoulli"){
  vec <- matrix(c(length(MLC), cases.expect, cases.observ))
  rownames(vec) <- c("Total observations in the MLC =  ","Expected cases in the MLC =","Observed cases in the MLC =")
  colnames(vec) <- ""
  scan <- list(method = paste("Scan test. Distribution: ",distr),
              fx = mfx, MLC = MLC, statistic = loglik, N = N, estimate = vec, nn = nn, nv = nv, coor = coor.input,
              p.value = p.value, nsim = nsim, data.name = data.name, distr = distr,
              case = case,
              alternative = alternative,
              cases.expect = cases.expect,
              cases.observ = cases.observ, sf = sf)
  }

  if (distr == "multinomial"){
    vec <- round(rbind(cases.expect,cases.observ),2)
    # names(vec) <- c("Number observ inside MLC","Expected cases in MLC","Observed cases in MLC")
    scan <- list(method = paste("Scan test. Distribution: ",distr),
                 fx = mfx, MLC = MLC, statistic = loglik, N = N, estimate = vec, nn = nn, nv = nv, coor = coor.input,
                 p.value = p.value, nsim = nsim, data.name = data.name, distr = distr,
                 cases.expect = cases.expect,
                 cases.observ = cases.observ, sf = sf)
  }


  class(scan) <- c("htest","scantest")
  scan

}