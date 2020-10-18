#' @importFrom colorRamps blue2red
#' @importFrom GISTools add.alpha
#' @importFrom graphics axis
#' @importFrom graphics layout
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics text
#' @importFrom graphics legend
#' @importFrom stats aggregate
#' @importFrom stats predict


#' @title Explainability and matchplot for PDP
#'
#' @description Computes explainability and matchplots for a partial dependence function.
#'
#' @param model   A model with corresponding predict function that returns numeric values.
#' @param x       Data frame.
#' @param vnames  Character vector of the variable set for which the patial dependence function is to be computed.
#' @param viz     Logical specifying whether a matchplot should be created.
#' @param ...     Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @return Numeric value with the explainability of the partial dependence function for the variables specified in \code{vnames}.
#'
#' @export
#'
#' @examples
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' xpy(boston.rf, boston, c("lstat", "rm"))
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname xpy
xpy <- function(model, x, vnames, viz = TRUE, ...){

  xs    <- x[, names(x) %in% vnames, drop = FALSE]
  xs    <- data.frame(ID = 1:nrow(xs), xs)
  xrest <- x[, !names(x) %in% vnames, drop = FALSE]

  xx    <- merge(xs, xrest, by.x = NULL, by.y = NULL)
  ID <- xx[,1]

  xx$yhat <- predict(model, xx[,-1], ...)
  pdps <- aggregate(xx$yhat, list(ID), mean)
  pred <- predict(model, x, ...)

  avpred <- mean(pred)

  #cbind(pdps$x, pred, avpred)
  if(viz){
    rnge <- range(c(range(pred), range(pdps$x)))
    plot(pdps$x, pred, xlim = rnge, ylim = rnge, xlab = "PDP", ylab = "Prediction", pch = 4, main = "PDP vs. Predictions")
    lines(rnge, rnge, lty = "dotted")
  }

  ASE <- mean((pred - pdps$x)^2)
  ASE0 <- mean((pred - avpred)^2)
  xty <- 1 - ASE / ASE0
  xty
}


#' @title Forward variable selection for PDP explanation
#'
#' @description Computes forward variable selection for partial dependence function based on explainability.
#'
#' @param model   A model with corresponding predict function that returns numeric values.
#' @param x       Data frame.
#' @param target  Character specifying the name of the (numeric) target variable (must be contained in data).
#' @param ...     Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @return Object of class \code{vsexp} containing the following elements:
#' @return \item{selection.order}{Vector of variable names in order of their entrance to the PD function during the variable selection process.}
#' @return \item{explainability}{Explainabilities after the different iterations of the variable selection.}
#' @return \item{details}{A data frame containing the explainabilities for all variables (rows) if they were entered in the model at each step (columns).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' vs <- fw.xpy(boston.rf, boston, "cmedv")
#' vs
#'
#' plot(vs)
#'}
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname fw.xpy
fw.xpy <- function(model, x, target, ...){

  # Initialization and selection of first variable
  n <- 1
  cat("Step", n, "\n")

  sel   <- NULL
  trace <- NULL
  nms <- nms.full <- names(x)[-which(names(x) == target)]
  xpys <- rep(NA, length(nms))
  names(xpys) <- nms

  for(v in nms) xpys[which(names(xpys) == v)] <- xpy(model, x, v, viz = F, ...)

  sel   <- c(sel, which.max(xpys))
  trace <- c(trace, max(xpys, na.rm = T))
  print(xpys)
  cat("\n", nms.full[sel], max(xpys, na.rm = T), "\n\n")

  # forward selection variables such that explainability is maximized
  while(length(nms) > 1){
    n <- n + 1
    cat("Step", n, "\n")

    nms <- nms.full[-sel]
    xpys <- cbind(xpys, NA)
    for(v in nms) xpys[which(rownames(xpys) == v), ncol(xpys)] <- xpy(model, x, c(names(sel), v), viz = F, ...)

    sel <- c(sel, which.max(xpys[,ncol(xpys)]))
    colnames(xpys) <- c(paste("Step", 1:n))
    trace <- c(trace, max(xpys, na.rm = T))

    print(xpys)
    cat("\n", nms.full[sel], max(xpys[,ncol(xpys)], na.rm=T), "\n\n")
  }

  res <- list(selection.order = sel, explainability = trace, details = xpys)
  class(res) <- "vsexp"
  return(res)
}


#' @export
plot.vsexp <- function(x, ...){
  plot(0:length(x$explainability), c(0,x$explainability), type = "l", xaxt = "n", xlab = "", ylim = c(0, 1), ylab = "explainability")
  axis(1, at = 1:length(x$selection.order), labels = names(x$selection.order), las = 2)
}


#' @export
print.vsexp <- function(x, ...) print(cbind(x$selection.order, x$explainability))



#' @title Explanation gap visualization
#'
#' @description Visualization of 2D PDP vs. unexplained residual predictions.
#'
#' @param model   A model with corresponding predict function that returns numeric values.
#' @param x       Data frame.
#' @param vnames  Character vector of the variable set for which the patial dependence function is to be computed.
#' @param type    Character, either \code{"pdp"}, \code{"gap"} or \code{"both"} specifying the meaning of the colours
#' of scatter: either the value of the PD function or the differences between PD and the model's predcitions.
#' In case of \code{"both"} two plots are created.
#' @param depth   Integer specifiying the number colours in the heat map.
#' @param alpha   Numeric value for alpha blending of the points in the scatter plot.
#' @param right   Position where to place the legend relative to the range of the x axis.
#' @param top     Position where to place the legend relative to the range of the y axis.
#' @param digits  Nuber of digits for rounding in the legend.
#' @param ...     Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' pdp2d(boston.rf, boston, c("lstat", "rm"), type = "both")
#' }
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname pdp2d
pdp2d <- function(model, x, vnames, type = "pdp", depth = 21, alpha = 2/3, right = 0.8, top = 0.95, digits = 1, ...){
  if(sum(names(x) %in% vnames) != 2) stop("You should use 2 Variables in order to compute scatterplots!")

  xs    <- x[, names(x) %in% vnames, drop = FALSE]
  xs    <- data.frame(ID = 1:nrow(xs), xs)
  xrest <- x[, !names(x) %in% vnames, drop = FALSE]

  xx    <- merge(xs, xrest, by.x = NULL, by.y = NULL)
  ID <- xx[,1]

  xx$yhat <- predict(model, xx[,-1], ...)
  pdps <- aggregate(xx$yhat, list(ID), mean)
  pred <- predict(model, x, ...)

  if(type == "both") par(mfrow = c(1,2))
  if(type != "gap"){
    minmax <- range(pdps$x)
    cols <- seq(minmax[1], minmax[2], length.out = depth)
    cols <- sapply(pdps$x, function(z) which.min(abs(z-cols)))
    cols <- colorRamps::blue2red(depth)[cols]
    cols <- GISTools::add.alpha(cols, alpha)
    plot(x[[vnames[1]]], x[[vnames[2]]], pch = 16, col = cols, xlab = vnames[1], ylab = vnames[2], main = "PD(x)")

    # legend
    vals   <- seq(minmax[1], minmax[2], length.out = depth)
    cols   <- colorRamps::blue2red(depth)
    steps  <- seq(1,depth, length.out = 5)
    xpos   <- range(xs[,3]) %*% c(1-right, right)
    ypos   <- range(xs[,2]) %*% c(1-top, top)
    legend(xpos, ypos, round(vals[steps],1),  fill = cols[steps], bty = "n", cex = 0.7)

  }
  if(type != "pdp"){
    diffs <- pdps$x - pred
    cols <- seq(-max(abs(range(diffs))), max(abs(range(diffs))), length.out = depth)
    cols <- sapply(diffs, function(z) which.min(abs(z-cols)))
    cols <- colorRamps::blue2red(depth)[cols]
    cols <- GISTools::add.alpha(cols, alpha)
    plot(x[[vnames[1]]], x[[vnames[2]]], pch = 16, col = cols, xlab = vnames[1], ylab = vnames[2], main = "PD(x) - f(x)")

    # legend
    vals   <- seq(-max(abs(range(diffs))), max(abs(range(diffs))), length.out = depth)
    cols   <- colorRamps::blue2red(depth)
    steps  <- seq(1,depth, length.out = 5)
    xpos   <- range(xs[,3]) %*% c(1-right, right)
    ypos   <- range(xs[,2]) %*% c(1-top, top)
    legend(xpos, ypos, round(vals[steps],1),  fill = cols[steps], bty = "n", cex = 0.7)

  }

}


#' @title Scatterplot matrix of 2D partial dependence plots
#'
#' @description Creates a scatterplot matrix of 2D partial dependence plots.
#'
#' @param model   A model with corresponding predict function that returns numeric values.
#' @param x       Data frame.
#' @param vnames  Character vector of the variable set for which the patial dependence function is to be computed.
#' @param depth   Integer specifiying the number colours in the heat map.
#' @param alpha   Numeric value for alpha blending of the points in the scatter plot.
#' @param ...     Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' PDPmatrix(boston.rf, boston, vnames = c("lstat", "rm", "lon", "nox"))
#' }
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname PDPmatrix
PDPmatrix <- function(model, x, vnames, depth = 21, alpha = 2/3, ...){
  n <- length(vnames)
  mat <- matrix(0,n,n)
  k <- 1
  for(i in 1:n){
    for(j in i:n){
      mat[i,j] <- k
      k <- k+1
    }
  }
  layout(mat, rep(2,n), rep(2,n))

  pred <- predict(model, x, ...)
  minmax <- range(pred)
  colrefs <- seq(minmax[1], minmax[2], length.out = depth)

  for(i in 1:(n-1)){

    plot(5, 5, type="n", axes=FALSE, ann=FALSE, xlim=c(0, 10), ylim = c(0,10))
    text(5,5, vnames[i], cex = 2)

    for(j in (i+1):n){
      xs    <- x[, names(x) %in% vnames[c(i,j)], drop = FALSE]
      xs    <- data.frame(ID = 1:nrow(xs), xs)
      xrest <- x[, !names(x) %in% vnames[c(i,j)], drop = FALSE]

      xx    <- merge(xs, xrest, by.x = NULL, by.y = NULL)
      ID <- xx[,1]

      xx$yhat <- predict(model, xx[,-1], ...)
      pdps <- aggregate(xx$yhat, list(ID), mean)

      cols <- sapply(pdps$x, function(z) which.min(abs(z-colrefs)))
      cols <- colorRamps::blue2red(depth)[cols]
      cols <- GISTools::add.alpha(cols, alpha)
      par(mar=c(1,1,1,1))
      plot(x[[vnames[j]]], x[[vnames[i]]], pch = 16, col = cols, xlab = vnames[1], ylab = vnames[2], main = "")
    }
  }

  plot(5, 5, type="n", axes=FALSE, ann=FALSE, xlim=c(0, 10), ylim = c(0,10))
  text(5,5, vnames[n], cex = 2)
}


