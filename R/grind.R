# grindR: Phase plane analysis in R
# Original author: Rob de Boer, Utrecht University
# Packaged by: Bram van Dijk, Utrecht University

grind_version <- "23-07-2026"

#' grindR: Phase plane analysis and parameter estimation for ODE models
#'
#' @keywords internal
#' @importFrom graphics arrows contour legend lines mtext par plot points segments
#' @importFrom stats median na.omit quantile sd setNames
#' @importFrom deSolve ode dede
#' @importFrom rootSolve steady
#' @importFrom FME modFit modCost
"_PACKAGE"

# `s`, `p`, `model` and `data` are supplied by the user in the global
# environment and used as default argument values; declare them so R CMD check
# does not flag them as undefined globals.
utils::globalVariables(c("s", "p", "model", "data"))

# Package-level options and mutable state -----------------------------------
#
# grindR reproduces the classic sourced-script idioms WITHOUT writing to the
# user's global environment (which CRAN forbids). Two mechanisms:
#
#  1. User-facing options (the colour palette, plot fonts, legend size) are
#     exported package data, so they exist after library(grindR). Functions
#     read them through getters -- .pal(), .sizeLegend(), .fontMain(),
#     .fontSub() -- that PREFER a copy in the user's global environment if one
#     exists. So the classic idiom keeps working:
#
#         colors[1] <- "blue"      # lands in the user's session
#         plane(...); run(...)     # getters read it back -> palette updated
#
#     The package only READS .GlobalEnv, never writes it, so this is CRAN-safe.
#
#  2. Internal cross-call state (which axes plane()/continue() last used) lives
#     in a private package environment, pkg_env, mutated in place. This
#     replaces the original <<- globals; users never see or touch it.
#
# The argument-routing lists (args_*, methods_*) are read-only constants
# defined at the bottom of this file, after all functions exist.

#' Default colour palette used by grindR
#'
#' Nullclines, trajectories and legends are coloured by state-variable index
#' using this palette. Override it in your session with a plain assignment,
#' e.g. \code{colors[1] <- "blue"}, and grindR will use your version.
#'
#' @format A character vector of colour names/numbers.
#' @export
colors <- c("blue","darkorange","darkgreen","red","darkmagenta","gold",
            "darkorchid","aquamarine","deeppink","gray", seq(2,991))

#' Length of the default \code{\link{colors}} palette
#' @export
ncolors <- length(colors)

#' Legend size, as a fraction of R's default (passed as \code{cex})
#' @export
sizeLegend <- 0.75

#' Font for plot titles (1 = plain, 2 = bold); see \code{\link[graphics]{par}}
#' @export
font.main <- 1

#' Font for plot subtitles (1 = plain, 2 = bold); see \code{\link[graphics]{par}}
#' @export
font.sub <- 1

# Private cross-call state (internal, not exported) -------------------------
# Mutated in place with pkg_env$foo <- ...; replaces the original <<- globals.
pkg_env <- new.env(parent = emptyenv())
# plane() axis state
pkg_env$x_plane    <- 1;  pkg_env$xmin_plane      <- -0.001; pkg_env$xmax_plane      <- 1.05
pkg_env$y_plane    <- 2;  pkg_env$ymin_plane      <- -0.001; pkg_env$ymax_plane      <- 1.05
pkg_env$log_plane  <- ""; pkg_env$addone_plane    <- FALSE
# continue() axis state
pkg_env$x_continue <- 1;  pkg_env$xmin_continue   <- 0;      pkg_env$xmax_continue   <- 1
pkg_env$y_continue <- 1;  pkg_env$ymin_continue   <- 0;      pkg_env$ymax_continue   <- 1
pkg_env$log_continue <- "";pkg_env$addone_continue <- FALSE

# Option getters (internal) -------------------------------------------------
# Return `name` from the user's global environment if they defined it there,
# otherwise the package default. This is what makes user overrides such as
# `colors[1] <- "blue"` visible while keeping the package CRAN-safe (read-only
# access to .GlobalEnv; inherits = FALSE so we read the user's copy, not the
# package's own exported default sitting further up the search path).
.gget <- function(name, default) {
  if (exists(name, envir = globalenv(), inherits = FALSE))
    return(get(name, envir = globalenv()))
  default
}
.pal        <- function() .gget("colors",     colors)
.sizeLegend <- function() .gget("sizeLegend", sizeLegend)
.fontMain   <- function() .gget("font.main",  font.main)
.fontSub    <- function() .gget("font.sub",   font.sub)

# Helper: coordinate vector for log/linear axes ----------------------------

plane_coord <- function(Log, Min, Max, npixels) {
  if (Log) return(10^seq(log10(Min), log10(Max), length.out = npixels))
  return(seq(Min, Max, length.out = npixels))
}

# plane() ------------------------------------------------------------------

#' Draw a phase plane with nullclines, vector field and trajectories
#'
#' \code{plane()} plots a two-dimensional phase plane for two state variables.
#' It draws the nullcline of each variable (the curve where its derivative is
#' zero); intersections of nullclines are equilibria. Optionally it overlays a
#' vector field (\code{vector}) and/or a phase portrait of trajectories
#' (\code{portrait}). The chosen axes are remembered, so a subsequent
#' \code{run(traject = TRUE)} or \code{plane(add = TRUE)} draws onto the same
#' plane.
#'
#' @param xmin,xmax Limits of the horizontal axis (the \code{x} variable).
#' @param ymin,ymax Limits of the vertical axis (the \code{y} variable).
#' @param xlab,ylab Axis labels; default to the state-variable names.
#' @param log Which axes to draw on a log scale: \code{""}, \code{"x"},
#'   \code{"y"} or \code{"xy"}.
#' @param npixels Grid resolution used to compute the nullclines. Higher is
#'   smoother but slower.
#' @param state Named numeric vector of state values. Non-axis variables are
#'   held at these values (see \code{zero}). Defaults to the global \code{s}.
#' @param parms Named numeric vector of parameters. Defaults to the global
#'   \code{p}.
#' @param odes The model function \code{f(t, state, parms)}. Defaults to the
#'   global \code{model}.
#' @param x,y State variable (index or name) on the horizontal and vertical
#'   axis. Default \code{1} and \code{2}.
#' @param time Time at which the derivatives are evaluated (for
#'   non-autonomous models).
#' @param grid Number of points per axis for the vector field and phase
#'   portrait.
#' @param show Names or indices of the variables whose nullclines are drawn
#'   (default: the two axis variables).
#' @param addone If \code{TRUE}, plot \code{variable + 1} so a log axis can
#'   include zero.
#' @param portrait If \code{TRUE}, integrate and draw trajectories from a grid
#'   of starting points.
#' @param vector Vector field to overlay: \code{0} none (default), \code{1}
#'   sign-only arrows (one per axis), \code{2} a single arrow in the flow
#'   direction, \code{3} sign-only segments without arrowheads.
#' @param add If \code{TRUE}, draw onto the existing phase plane (reusing its
#'   axes and limits) instead of starting a new one.
#' @param legend If \code{TRUE}, draw a legend of the variable colours.
#' @param zero If \code{TRUE}, set all non-axis state variables to zero before
#'   computing the plane; if \code{FALSE}, hold them at their \code{state}
#'   values.
#' @param lwd,col,pch Line width of the nullclines, colour of the phase-portrait
#'   trajectories, and plotting character of their starting points.
#' @param vectorlen Scaling factor for the length of the vectors.
#' @param arrowsize Arrowhead length as a fraction of the vector's shaft
#'   (so heads scale with the arrows); clamped to a visible range.
#' @param ... Additional arguments passed on to \code{run} (for the phase
#'   portrait) or to \code{plot}.
#'
#' @return None; \code{plane} is called for the plot it draws.
#'
#' @examples
#' model <- function(t, state, parms) {
#'   with(as.list(c(state, parms)), {
#'     dR <- b * R * (1 - R / K) - a * R * N / (R + h)
#'     dN <- c * a * R * N / (R + h) - delta * N
#'     list(c(dR, dN))
#'   })
#' }
#' p <- c(b = 1, K = 2, a = 1, c = 1, delta = 0.3, h = 1)
#' s <- c(R = 1, N = 1)
#'
#' plane(xmax = 2, ymax = 2)               # nullclines only
#' plane(xmax = 2, ymax = 2, vector = 1)   # add sign-only arrows
#' run(traject = TRUE)                     # overlay a trajectory from s
#'
#' @seealso \code{\link{run}}, \code{\link{newton}}, \code{\link{continue}}
#'
#' @export
plane <- function(xmin=-0.001, xmax=1.05, ymin=-0.001, ymax=1.05, xlab="", ylab="", log="", npixels=500, state=s, parms=p, odes=model, x=1, y=2, time=0, grid=5, show=NULL, addone=FALSE, portrait=FALSE, vector=0, add=FALSE, legend=TRUE, zero=TRUE, lwd=2, col="black", pch=20, vectorlen=1, arrowsize=0.3, ...) {
  # Make a phase plane with nullclines and/or phase portrait
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
  }
  dots_run <- if (!is.null(dots)) dots[names(dots) %in% args_run] else NULL
  if (add) {
    x <- pkg_env$x_plane
    y <- pkg_env$y_plane
    xmin <- pkg_env$xmin_plane; xmax <- pkg_env$xmax_plane
    ymin <- pkg_env$ymin_plane; ymax <- pkg_env$ymax_plane
    log <- pkg_env$log_plane; addone <- pkg_env$addone_plane
  } else {
    if (!is.numeric(x)) x <- index(x,names(state))
    if (!is.numeric(y)) y <- index(y,names(state))
    pkg_env$x_plane <- x
    pkg_env$y_plane <- y
    pkg_env$xmin_plane <- xmin; pkg_env$xmax_plane <- xmax
    pkg_env$ymin_plane <- ymin; pkg_env$ymax_plane <- ymax
    pkg_env$log_plane <- log; pkg_env$addone_plane <- addone
  }
  ishows <- if (!is.null(show)) index(show, names(state)) else c(x, y)
  nvar <- length(state)
  if (zero) state[1:nvar] <- rep(0,nvar)
  lvec <- 15                         # length of vector
  logx <- ifelse(grepl('x',log), TRUE, FALSE)
  logy <- ifelse(grepl('y',log), TRUE, FALSE)
  xc <- plane_coord(logx,xmin,xmax,npixels)
  yc <- plane_coord(logy,ymin,ymax,npixels)
  if (xlab == "") xlab <- names(state)[x]
  if (ylab == "") ylab <- names(state)[y]
  if (addone) {
    if (logx) xlab <- paste(xlab,"+ 1")
    if (logy) ylab <- paste(ylab,"+ 1")
  }
  if (!add) {
    do.call('plot',c(list(1,1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",log=log,font.main=.fontMain(),font.sub=.fontSub()),dots[names(dots) %in% args_plot]))
    mtext(xlab, side = 1, line = 2.5, col = .pal()[x], font = 2, cex = 1)
    mtext(ylab, side = 2, line = 2.5, col = .pal()[y], font = 2, cex = 1)
  }
  
  npixels2 <- npixels^2
  vstate <- as.list(state)
  vparms <- as.list(parms)
  vparms <- lapply(vparms,rep,vparms,npixels2)
  vstate <- lapply(vstate,rep,vstate,npixels2)
  vstate[[x]] <- rep.int(xc, npixels)
  vstate[[y]] <- rep.int(yc, rep.int(npixels, npixels))
  if (addone & logx) vstate[[x]] <- vstate[[x]] - 1
  if (addone & logy) vstate[[y]] <- vstate[[y]] - 1
  dvstate <- odes(time,vstate,vparms)[[1]]
  dim(dvstate) <- c(npixels,npixels,nvar)
  for (i in ishows) 
    contour(xc,yc,dvstate[,,i],levels=0,drawlabels=FALSE,add=TRUE,col=.pal()[i],lwd=lwd)
  
  if (portrait | vector) {
    
    pin <- par("pin") # get the dimension of the plot in R in inches

    dx <- if (logx) (log10(xmax)-log10(xmin))/grid else (xmax-xmin)/grid
    dy <- if (logy) (log10(ymax)-log10(ymin))/grid else (ymax-ymin)/grid

    # Isotropic shaft length in inches. It tracks the grid-cell size
    # (min(pin) * 0.5/grid) so vectors stay proportional to the field under
    # par(mfrow) and coarse grids, but is floored at min(pin)/lvec so they do
    # not shrink to nothing at fine grids. Both terms scale with the panel, so
    # the result is robust to the plotting layout.
    shaft_in <- min(pin) * max(0.5/lvec, 0.4/grid) * vectorlen
    # Arrowhead length scales with the shaft (arrowsize = head-to-shaft ratio),
    # clamped so heads stay visible but never dwarf the arrow.
    ahead    <- max(0.02, min(0.15, 0.8*shaft_in * arrowsize))
    vx <- if (logx) 1 + 3.32*grid*dx/lvec else shaft_in * (xmax-xmin)/pin[1]
    vy <- if (logy) 1 + 3.32*grid*dy/lvec else shaft_in * (ymax-ymin)/pin[2]
    
    for (i in seq(1,grid)) {
      state[x] <- ifelse(logx, 10^((i-1)*dx + dx/2 + log10(xmin)),
                         (i-1)*dx + dx/2 + xmin)
      for (j in seq(1,grid,1)) {
        state[y] <- ifelse(logy, 10^((j-1)*dy + dy/2 + log10(ymin)),
                           (j-1)*dy + dy/2 + ymin)
        if (portrait) {
          points(state[x],state[y],pch=pch)
          nsol <- do.call('run',c(list(state=state,parms=parms,odes=odes,timeplot=FALSE,table=TRUE),dots_run))
          lines(cbind(nsol[x+1],nsol[y+1]),col=col)
        }
        if (vector) {
          dv  <- odes(time, state, parms)[[1]]
          dvx <- dv[x]
          dvy <- dv[y]

          if (vector == 3) {
            # Two orthogonal segments, sign-only, no arrowhead
            dt <- sign(c(dvx, dvy))
            segments(
              x0 = state[x], y0 = state[y],
              x1 = if (logx) state[x] * vx^dt[1] else state[x] + vx * dt[1],
              y1 = state[y]
            )
            segments(
              x0 = state[x], y0 = state[y],
              x1 = state[x],
              y1 = if (logy) state[y] * vy^dt[2] else state[y] + vy * dt[2]
            )

          } else if (vector == 2) {
            # Single arrow in true (dx, dy) direction, normalised length
            mag <- sqrt(dvx^2 + dvy^2)
            if (mag > 0) {
              norm_dx <- 1.4*dvx / mag
              norm_dy <- 1.4*dvy / mag
            } else {
              norm_dx <- 0
              norm_dy <- 0
            }
            x1 <- if (logx) state[x] * vx^norm_dx else state[x] + vx * norm_dx
            y1 <- if (logy) state[y] * vy^norm_dy else state[y] + vy * norm_dy
            arrows(
              x0 = state[x], y0 = state[y],
              x1 = x1, y1 = y1,
              lwd = 1.2, length = ahead
            )

          } else {
            # Default (vector == 1, or TRUE): two orthogonal arrows, sign-only
            dt <- sign(c(dvx, dvy))
            arrows(
              x0 = state[x], y0 = state[y],
              x1 = if (logx) state[x] * vx^dt[1] else state[x] + vx * dt[1],
              y1 = state[y],
              lwd = 1.2, length = ahead
            )
            arrows(
              x0 = state[x], y0 = state[y],
              x1 = state[x],
              y1 = if (logy) state[y] * vy^dt[2] else state[y] + vy * dt[2],
              lwd = 1.2, length = ahead
            )
          }
        }
      }
    }
  }
  if(!add){
    if (legend)
      legend("topright",legend=names(state)[ishows],col=.pal()[ishows],lty=1,lwd=lwd,cex=.sizeLegend())
  }
}

# run() --------------------------------------------------------------------

#' Numerically integrate an ODE model over time
#'
#' \code{run()} solves a system of ordinary differential equations (via
#' \code{\link[deSolve]{ode}}), optionally plots the time course, and returns
#' the final state. Because it returns the end point as a named vector, the
#' result can be fed back into \code{run()}, \code{\link{plane}} or
#' \code{\link{newton}}.
#'
#' @param tmax Final time of the integration.
#' @param tstep Time increment between successive output points.
#' @param state Named numeric vector of initial values for the state
#'   variables. Defaults to the global \code{s}.
#' @param parms Named numeric vector of parameters. Defaults to the global
#'   \code{p}.
#' @param odes The model: a function \code{f(t, state, parms)} returning the
#'   derivatives as \code{list(c(...))}. Defaults to the global \code{model}.
#' @param ymin,ymax Lower and upper limits of the y-axis (\code{ymax = NULL}
#'   auto-scales).
#' @param log Which axes to draw on a log scale: \code{""}, \code{"x"},
#'   \code{"y"} or \code{"xy"}.
#' @param xlab,ylab Axis labels for the time plot.
#' @param tmin Start time of the integration.
#' @param draw Function used to draw the time course, typically
#'   \code{\link[graphics]{lines}} or \code{\link[graphics]{points}}.
#' @param times Optional explicit vector of output times. If supplied it
#'   overrides \code{tmin}, \code{tmax} and \code{tstep}.
#' @param show Names or indices of the state variables to plot (default: all).
#' @param arrest Times, or names of parameters holding times, at which the
#'   integrator is forced to stop exactly (useful at discontinuities). Cannot
#'   be combined with \code{events}.
#' @param events A \code{\link[deSolve]{events}} specification passed to the
#'   solver for discrete changes to the state.
#' @param after A string of R code evaluated after each time step; may modify
#'   \code{state} and \code{parms}. Cannot be combined with \code{solution} or
#'   \code{delay}.
#' @param tweak A string of R code evaluated once after integration, operating
#'   on the result data frame \code{nsol}.
#' @param timeplot If \code{TRUE}, draw the time course.
#' @param traject If \code{TRUE}, add the trajectory to the current phase plane
#'   (see \code{\link{plane}}) instead of drawing a time plot.
#' @param table If \code{TRUE}, return the full time series as a data frame
#'   instead of only the final state.
#' @param add If \code{TRUE}, add to the existing plot rather than starting a
#'   new one.
#' @param legend If \code{TRUE}, draw a legend on the time plot.
#' @param solution If \code{TRUE}, treat \code{odes} as an explicit solution
#'   \code{f(t)} evaluated directly, rather than integrating derivatives.
#' @param delay If \code{TRUE}, integrate delay differential equations with
#'   \code{\link[deSolve]{dede}}.
#' @param lwd,col,pch Line width, colour and plotting character for the
#'   trajectory / time plot.
#' @param ... Additional arguments passed on to the solver
#'   (e.g. \code{method}, \code{atol}, \code{rtol}) or to \code{plot}.
#'
#' @return By default, a named numeric vector with the state at \code{tmax}.
#'   If \code{table = TRUE}, a data frame with a \code{time} column followed by
#'   one column per state variable.
#'
#' @examples
#' model <- function(t, state, parms) {
#'   with(as.list(c(state, parms)), {
#'     dR <- b * R * (1 - R / K) - d * R
#'     list(c(dR))
#'   })
#' }
#' p <- c(b = 1, d = 0.1, K = 2)
#' s <- c(R = 0.1)
#'
#' run(tmax = 50)                     # integrate and plot the time course
#' run(tmax = 50, timeplot = FALSE)   # just the final state
#' run(table = TRUE)                  # the full time series as a data frame
#'
#' @seealso \code{\link{plane}}, \code{\link{newton}}, \code{\link{fit}}
#'
#' @export
run <- function(tmax=100, tstep=1, state=s, parms=p, odes=model,
                ymin=0, ymax=NULL, log="", xlab="Time", ylab="Density",
                tmin=0, draw=lines, times=NULL, show=NULL,
                arrest=NULL, events=NULL, after=NULL, tweak=NULL,
                timeplot=TRUE, traject=FALSE, table=FALSE,
                add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE,
                lwd=2, col="black", pch=20, ...) {
  if (delay & (solution | !is.null(after)))
    stop("Don't use solution or after with delay equations")
  if (delay) args_run <- args_run_dde
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run, args_plot)])
    if (length(unknown) > 0) warning(paste("Unknown argument(s):", unknown, sep=" "))
    dots_run <- dots[names(dots) %in% args_run]
  } else dots_run <- NULL
  nvar <- length(state)
  if (is.null(times)) {
    times <- seq(tmin, tmax, by=tstep)
  } else {
    times <- sort(times)
    tmin  <- min(times)
    tmax  <- max(times)
  }
  if (!is.null(arrest)) {
    if (!is.null(events)) stop("Don't combine the option arrest with events")
    if (!is.numeric(arrest)) arrest <- sort(as.numeric(parms[arrest]))
    nearby   <- deSolve::nearestEvent(arrest, times)
    nearby   <- nearby[nearby < arrest]
    lennear  <- length(nearby)
    if (lennear == 1 && nearby[1] == 0) lennear <- 0
    if (lennear > 0) {
      if (nearby[1] == 0) nearby <- nearby[2:lennear]
      arrest <- sort(unique(c(nearby, arrest)))
    }
    events <- list(func=dummyEvent, time=arrest)
    times  <- deSolve::cleanEventTimes(times, arrest, eps=.Machine$double.eps*10)
    times  <- sort(c(times, arrest))
  }
  if (solution) {
    if (!is.null(after)) stop("Don't combine the option after with solution")
    nsol <- sapply(times, odes, state, parms)
    if (is.list(nsol)) {
      nsol <- unlist(nsol)
      if (nvar > 1) dim(nsol) <- c(nvar, length(times))
    }
    if (nvar > 1) nsol <- data.frame(times, t(nsol))
    else          nsol <- data.frame(times, nsol)
    names(nsol) <- c("time", names(state))
  } else {
    if (is.null(after)) {
      nsol <- as.data.frame(
        do.call(if (!delay) 'ode' else 'dede',
                c(list(times=times, func=odes, y=state, parms=parms,
                        events=events), dots_run))
      )
    } else {
      keep <- state
      nsol <- t(sapply(seq(length(times)-1), function(i) {
        t <- times[i+1]
        f <- do.call('ode', c(list(times=c(times[i],t), func=odes,
                                    y=state, parms=parms), dots_run))
        dim(f) <- c(2, nvar+1)
        state[1:nvar] <- f[2, 2:(nvar+1)]
        eval(parse(text=after))
        parms <<- parms
        state <<- state
      }))
      if (nvar > 1) {
        nsol <- as.data.frame(cbind(times, rbind(as.numeric(keep), nsol)))
      } else {
        nsol <- as.data.frame(cbind(times, c(keep, nsol)))
      }
      names(nsol) <- c("time", names(state))
      state <- keep
    }
  }
  if (!is.null(tweak)) eval(parse(text=tweak))
  if (timeplot & !traject)
    do.call('timePlot', c(list(data=nsol, tmin=tmin, tmax=tmax,
                                ymin=ymin, ymax=ymax, log=log, add=add,
                                xlab=xlab, ylab=ylab, show=show, draw=draw,
                                lwd=lwd, legend=legend,
                                font.main=.fontMain(), font.sub=.fontSub()),
                           dots[names(dots) %in% args_plot]))
  if (traject) {
    points(nsol[1, pkg_env$x_plane+1], nsol[1, pkg_env$y_plane+1], pch=pch)
    lines(nsol[, pkg_env$x_plane+1], nsol[, pkg_env$y_plane+1], lwd=lwd, col=col)
  }
  if (table) return(nsol)
  f <- state
  f[1:length(f)] <- as.numeric(nsol[nrow(nsol), 2:(nvar+1)])
  return(f)
}

# newton() -----------------------------------------------------------------

#' Find and classify a steady state
#'
#' \code{newton()} locates a steady state (equilibrium) of the model near a
#' starting guess using Newton-Raphson root-finding
#' (\code{\link[rootSolve]{steady}}), then classifies its stability. The
#' Jacobian at the equilibrium is computed numerically by finite differences
#' (\code{\link[rootSolve]{jacobian.full}}) and its eigenvalues by
#' \code{\link[base]{eigen}}; the sign of the dominant eigenvalue determines
#' stability. For two-variable systems the point is labelled a stable/unstable
#' node, spiral, or saddle.
#'
#' @param state Named numeric vector used as the starting guess for the
#'   steady state. Defaults to the global \code{s}.
#' @param parms Named numeric vector of parameters. Defaults to the global
#'   \code{p}.
#' @param odes The model function \code{f(t, state, parms)}. Defaults to the
#'   global \code{model}.
#' @param time Time at which the derivatives are evaluated (for
#'   non-autonomous models).
#' @param positive If \code{TRUE}, restrict the search to non-negative state
#'   values.
#' @param jacobian If \code{TRUE}, also print the Jacobian matrix.
#' @param vector If \code{TRUE}, also print the eigenvectors.
#' @param plot If \code{TRUE}, mark the equilibrium on the current phase plane
#'   (filled circle if stable, open if unstable); see \code{\link{plane}}.
#' @param silent If \code{TRUE}, suppress all printing and instead return the
#'   full result as a list (see Value).
#' @param addone If \code{TRUE}, offset the plotted point by 1, to match a
#'   phase plane drawn with \code{plane(addone = TRUE)}.
#' @param ... Additional arguments passed on to
#'   \code{\link[rootSolve]{steady}}.
#'
#' @return If a steady state is found: by default the equilibrium as a named
#'   numeric vector (with the stability classification and eigenvalues printed).
#'   If \code{silent = TRUE}, a list with components \code{state} (the
#'   equilibrium), \code{jacobian}, \code{values} (eigenvalues) and
#'   \code{vectors} (eigenvectors). Returns \code{NULL} if the solver does not
#'   converge.
#'
#' @examples
#' model <- function(t, state, parms) {
#'   with(as.list(c(state, parms)), {
#'     dR <- b * R * (1 - R / K) - a * R * N / (R + h)
#'     dN <- c * a * R * N / (R + h) - delta * N
#'     list(c(dR, dN))
#'   })
#' }
#' p <- c(b = 1, K = 2, a = 1, c = 1, delta = 0.3, h = 1)
#' s <- c(R = 1, N = 1)
#'
#' newton(c(R = 0.5, N = 0.5))            # find and classify a steady state
#' eq <- newton(s, silent = TRUE)         # return state, Jacobian, eigen-data
#'
#' @seealso \code{\link{plane}}, \code{\link{continue}}
#'
#' @export
newton <- function(state=s, parms=p, odes=model, time=0,
                   positive=FALSE, jacobian=FALSE, vector=FALSE,
                   plot=FALSE, silent=FALSE, addone=FALSE, ...) {
  One <- ifelse(addone, 1, 0)
  q   <- rootSolve::steady(y=state, func=odes, parms=parms,
                            time=time, positive=positive, ...)
  if (attr(q, "steady")) {
    equ <- q$y
    equ <- ifelse(abs(equ) < 1e-8, 0, equ)
    jac <- rootSolve::jacobian.full(y=equ, func=odes, parms=parms)
    eig <- eigen(jac)
    dom <- max(Re(eig$values))
    if (!silent) {
      print(equ)
      if (length(equ) == 2) {
        if (is.complex(eig$values[1])) {
          cat(ifelse(dom < 0, "Stable spiral point, ", "Unstable spiral point, "))
        } else if (prod(eig$values) > 0) {
          cat(ifelse(dom < 0, "Stable node, ", "Unstable node, "))
        } else cat("Saddle point (unstable), ")
      } else {
        cat(ifelse(dom < 0, "Stable point, ", "Unstable point, "))
      }
      cat("eigenvalues:\n")
      print(eig$values)
    }
    if (vector)   { cat("Eigenvectors:\n"); print(eig$vectors) }
    if (jacobian) { cat("Jacobian:\n");     print(jac) }
    if (plot)
      points(equ[pkg_env$x_plane] + One, equ[pkg_env$y_plane] + One,
             pch=ifelse(dom < 0, 19, 1))
    if (silent) return(list(state=equ, jacobian=jac,
                             values=eig$values, vectors=eig$vectors))
    return(equ)
  }
  cat("No convergence: start closer to a steady state")
  return(NULL)
}

# continue() ---------------------------------------------------------------

#' Trace a steady state along a parameter (bifurcation diagram)
#'
#' \code{continue()} follows a steady state as one parameter is slowly varied
#' and plots one state variable against that parameter, producing a bifurcation
#' diagram. Starting from an equilibrium, it uses natural-parameter
#' continuation: at each step it nudges the parameter and re-solves for the
#' steady state (\code{\link[rootSolve]{steady}}) using the previous point as
#' the initial guess. Each branch segment is coloured by stability (the sign of
#' the dominant eigenvalue of the finite-difference Jacobian), and detected
#' bifurcations and turning points are reported to the console.
#'
#' Natural-parameter continuation can lose a branch at folds (where the branch
#' turns back in the parameter) or where two branches pass close together; if a
#' branch disappears unexpectedly, try a smaller \code{step}, restart from a
#' point on the missing branch, or use \code{add = TRUE} to trace it separately.
#'
#' @param state Named numeric vector at (or very near) a steady state, used as
#'   the starting point. Defaults to the global \code{s}.
#' @param parms Named numeric vector of parameters. Defaults to the global
#'   \code{p}.
#' @param odes The model function \code{f(t, state, parms)}. Defaults to the
#'   global \code{model}.
#' @param step Initial continuation step, as a fraction of \code{xmax} (or a
#'   multiplicative factor when the x-axis is logarithmic). Reduced
#'   automatically in difficult regions.
#' @param x Parameter to vary (index or name); the horizontal axis. Default 1.
#' @param y State variable to plot (index or name); the vertical axis.
#'   Default 2.
#' @param time Time at which the derivatives are evaluated (for
#'   non-autonomous models).
#' @param xmin,xmax Range of the parameter (horizontal axis) to scan.
#' @param ymin,ymax Range of the plotted variable (vertical axis).
#' @param xlab,ylab Axis labels; default to the parameter and variable names.
#' @param log Which axes to draw on a log scale: \code{""}, \code{"x"},
#'   \code{"y"} or \code{"xy"}.
#' @param col Length-3 vector of branch colours indexed by the sign of the
#'   dominant eigenvalue: stable (\code{-}), neutral (\code{0}) and unstable
#'   (\code{+}).
#' @param lwd Length-3 vector of line widths, indexed as \code{col}.
#' @param addone If \code{TRUE}, plot \code{variable + 1} so a log axis can
#'   include zero.
#' @param positive If \code{TRUE}, restrict the steady-state search to
#'   non-negative state values.
#' @param nvar If \code{TRUE}, colour branches by the number of non-zero
#'   (surviving) state variables instead of by stability.
#' @param add If \code{TRUE}, draw onto the existing bifurcation diagram
#'   (reusing its axes), e.g. to trace another branch.
#' @param ... Additional arguments passed on to
#'   \code{\link[rootSolve]{steady}}.
#'
#' @return \code{NULL}, invisibly; \code{continue} is called for the diagram it
#'   draws (and the bifurcation/turning points it prints).
#'
#' @examples
#' model <- function(t, state, parms) {
#'   with(as.list(c(state, parms)), {
#'     dR <- b * R * (1 - R / K) - a * R * N / (R + h)
#'     dN <- c * a * R * N / (R + h) - delta * N
#'     list(c(dR, dN))
#'   })
#' }
#' p <- c(b = 1, K = 2, a = 1, c = 1, delta = 0.3, h = 1)
#' s <- c(R = 1, N = 1)
#'
#' f <- newton(c(R = 0.5, N = 0.5), silent = TRUE)$state  # a steady state
#' continue(f, x = "K", xmax = 3)                         # vary K, plot N
#'
#' @seealso \code{\link{newton}}, \code{\link{plane}}
#'
#' @export
continue <- function(state=s, parms=p, odes=model, step=0.01,
                     x=1, y=2, time=0,
                     xmin=0, xmax=1, ymin=0, ymax=1.05,
                     xlab="", ylab="", log="",
                     col=c("red","black","blue"), lwd=c(2,1,1),
                     addone=FALSE, positive=FALSE, nvar=FALSE,
                     add=FALSE, ...) {
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_steady, args_plot)])
    if (length(unknown) > 0) warning(paste("Unknown argument(s):", unknown, sep=" "))
    dots_steady <- dots[names(dots) %in% args_steady]
  } else dots_steady <- NULL
  if (add) {
    x <- pkg_env$x_continue; y <- pkg_env$y_continue
  } else {
    if (!is.numeric(x)) x <- index(x, names(parms))
    if (!is.numeric(y)) y <- index(y, names(state))
    pkg_env$x_continue <- x; pkg_env$y_continue <- y
  }
  p0  <- parms[x]
  q0  <- do.call('steady', c(list(y=state, func=odes, parms=parms,
                                   time=time, positive=positive), dots_steady))
  if (!attr(q0, "steady"))
    stop("No convergence: start closer to a steady state")
  cat("Starting at", names(parms[x]), "=", parms[x], "with:\n")
  print(q0$y)
  bary <- q0$y[y]
  if (!add) {
    if (missing(xmax) & parms[x] >= 1) xmax <- 2*parms[x]
    if (missing(xmin) & parms[x] <  0) xmin <- 2*parms[x]
    if (!missing(xmin) & xmin >= parms[x]) stop("xmin should be smaller than parameter")
    if (!missing(xmax) & xmax <= parms[x]) stop("xmax should be larger than parameter")
    if (missing(ymax) & bary >= 1.05) ymax <- 2*bary
    if (missing(ymin) & bary <  0)    ymin <- 2*bary
    if (!missing(ymin) & ymin >= bary & !addone) stop("ymin should be smaller than y-variable")
    if (!missing(ymax) & ymax <= bary)           stop("ymax should be larger than y-variable")
    if (xlab == "") xlab <- names(p0)
    if (ylab == "") {
      ylab <- names(state)[y]
      if (addone) ylab <- paste(ylab, "+ 1")
    }
    do.call('plot', c(list(1, 1, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
                            xlab=xlab, ylab=ylab, log=log,
                            font.main=.fontMain(), font.sub=.fontSub()),
                      dots[names(dots) %in% args_plot]))
    pkg_env$xmin_continue <- xmin; pkg_env$xmax_continue <- xmax
    pkg_env$ymin_continue <- ymin; pkg_env$ymax_continue <- ymax
    pkg_env$log_continue  <- log;  pkg_env$addone_continue <- addone
  } else {
    xmin <- pkg_env$xmin_continue; xmax <- pkg_env$xmax_continue
    ymin <- pkg_env$ymin_continue; ymax <- pkg_env$ymax_continue
    log  <- pkg_env$log_continue;  addone <- pkg_env$addone_continue
    if (ymin >= bary & !addone) stop("Initial point below minimum of current y-axis")
    if (ymax <= bary)           stop("Initial point above maximum of current y-axis")
  }
  logx <- ifelse(grepl('x', log), TRUE, FALSE)
  COL <- function(s, i) {
    if (!nvar) return(col[i])
    return(col[length(s[s > 1e-9])])
  }
  FUN <- function(lastState, lastDom, step) {
    lastP <- p0
    preLastState <- lastState
    nok <- 0
    while (xmin < lastP & lastP < xmax &
           ymin < lastState[y]+One & lastState[y] < ymax) {
      parms[x] <- ifelse(logx, lastP*(1 + step), lastP + step)
      q <- do.call('steady', c(list(y=lastState, func=odes, parms=parms,
                                     time=time, positive=positive), dots_steady))
      newState <- q$y
      if (attr(q,"steady") &
          sum(abs(newState-lastState))/(1e-9+sum(abs(lastState))) < 0.1) {
        jac <- rootSolve::jacobian.full(y=newState, func=odes, parms=parms)
        dom <- sign(max(Re(eigen(jac)$values)))
        if (dom != lastDom)
          cat("Bifurcation at", names(parms[x]), "=", parms[x], "\n")
        lines(c(if(logx) parms[x]/(1+step) else parms[x]-step, parms[x]),
              c(lastState[y]+One, newState[y]+One),
              col=COL(lastState, dom+2), lwd=lwd[dom+2])
        preLastState <- lastState
        lastState    <- newState
        lastDom      <- dom
        lastP        <- parms[x]
        if (nok > 10 & abs(step) < actualStep)
          step <- sign(step)*min(2*abs(step), actualStep)
        nok <- nok + 1
      } else {
        nok <- 0
        if (abs(step) > actualStep/100) {
          step <- step/2
        } else {
          parms[x]   <- lastP
          predState  <- lastState + 5*(lastState-preLastState)
          q <- do.call('steady', c(list(y=predState, func=odes, parms=parms,
                                         time=time, positive=positive), dots_steady))
          newState <- q$y
          if (attr(q,"steady") &
              sum(abs(newState-lastState))/(1e-9+sum(abs(lastState))) > 0.001) {
            cat("Turning point at", names(parms[x]), "=", parms[x], "\n")
            jac <- rootSolve::jacobian.full(y=newState, func=odes, parms=parms)
            dom <- sign(max(Re(eigen(jac)$values)))
            middle <- (lastState[y]+newState[y])/2
            lines(c(parms[x],parms[x]),
                  c(lastState[y]+One, middle+One),
                  col=COL(lastState, lastDom+2), lwd=lwd[lastDom+2])
            lines(c(parms[x],parms[x]),
                  c(middle+One, newState[y]+One),
                  col=COL(newState, dom+2), lwd=lwd[dom+2])
            step         <- -step
            preLastState <- lastState
            lastState    <- newState
            lastDom      <- dom
            lastP        <- parms[x]
          } else {
            cat("Final point at", names(parms[x]), "=", parms[x], "\n")
            cat("If this looks wrong try changing the step size\n")
            break
          }
        }
      }
    }
  }
  One     <- ifelse(addone, 1, 0)
  orgWarn <- getOption("warn")
  options(warn = -1)
  jac        <- rootSolve::jacobian.full(y=q0$y, func=odes, parms=parms)
  dom        <- sign(max(Re(eigen(jac)$values)))
  actualStep <- if(logx) step else step*xmax
  FUN(lastState=q0$y, lastDom=dom,  actualStep)
  FUN(lastState=q0$y, lastDom=dom, -actualStep)
  options(warn = orgWarn)
  return(NULL)
}

# fit() --------------------------------------------------------------------

#' Estimate parameters by fitting a model to data
#'
#' \code{fit()} estimates free parameters and/or initial state values of an ODE
#' model by fitting it to one or more time-series datasets, minimising the sum
#' of squared residuals with \code{\link[FME]{modFit}} (via \code{\link{cost}}).
#' It can fit parameters shared across datasets (\code{free}), parameters that
#' take a separate value per dataset (\code{differ}), and per-dataset values
#' held fixed (\code{fixed}); it can fit on a log scale, draw the fit against
#' the data, and bootstrap confidence intervals.
#'
#' @param datas A data frame, or a list of data frames, of observations. The
#'   first column is time; the remaining columns are observed variables, matched
#'   by name to the state variables. Defaults to the global \code{data}.
#' @param state Named numeric vector of initial state values. Defaults to the
#'   global \code{s}.
#' @param parms Named numeric vector of parameters. Defaults to the global
#'   \code{p}.
#' @param odes The model function \code{f(t, state, parms)}. Defaults to the
#'   global \code{model}.
#' @param free Names of parameters and/or initial state values to estimate
#'   (shared across all datasets). Defaults to all when neither \code{free} nor
#'   \code{differ} is given.
#' @param who Deprecated alias for \code{free}.
#' @param differ Names of parameters/states to estimate separately for each
#'   dataset (one value per dataset); a character vector, or a named list of
#'   starting values.
#' @param fixed Named list of per-dataset values that are held fixed (not
#'   estimated); one value per dataset.
#' @param tmin,tmax Time range of the fit plot.
#' @param ymin,ymax Vertical range of the fit plot (\code{NULL} auto-scales).
#' @param log Which axes to draw on a log scale: \code{""}, \code{"x"},
#'   \code{"y"} or \code{"xy"}.
#' @param xlab,ylab Axis labels of the fit plot.
#' @param bootstrap Number of bootstrap resamples for confidence intervals
#'   (\code{0} = none).
#' @param show Names or indices of the variables to plot.
#' @param fun Optional function applied to both the data and the model output
#'   before the cost is computed (e.g. a transform).
#' @param costfun Cost function to minimise; defaults to \code{\link{cost}}.
#' @param logpar If \code{TRUE}, estimate parameters on a log scale (keeping
#'   them positive).
#' @param lower,upper Lower and upper bounds on the estimated parameters.
#' @param initial If \code{TRUE}, take each dataset's initial state from its
#'   first row (which must be at time 0) rather than estimating it.
#' @param add If \code{TRUE}, add the fit plot to an existing one.
#' @param timeplot If \code{TRUE}, plot the fitted model against the data.
#' @param legend If \code{TRUE}, draw a legend.
#' @param main,sub Plot title and subtitle; may be one per dataset.
#' @param pchMap Optional mapping of data columns to plotting characters.
#' @param ... Additional arguments passed on to \code{\link[FME]{modFit}},
#'   \code{run} or \code{plot} (e.g. \code{method}).
#'
#' @return The \code{\link[FME]{modFit}} object (a list with the estimates in
#'   \code{$par}, the residual sum of squares in \code{$ssr}, and so on). If
#'   \code{bootstrap > 0}, a \code{$bootstrap} matrix of resampled estimates is
#'   added. The estimates and SSR are also printed.
#'
#' @examples
#' model <- function(t, state, parms) {
#'   with(as.list(c(state, parms)), list(c(r * N * (1 - N / K))))
#' }
#' p <- c(r = 0.5, K = 10)
#' s <- c(N = 0.1)
#'
#' set.seed(1)
#' data <- run(tmax = 20, table = TRUE, timeplot = FALSE)   # synthetic data
#' data$N <- data$N + rnorm(nrow(data), sd = 0.2)
#' fit(data, free = c("r", "K"))
#'
#' @seealso \code{\link{run}}, \code{\link{cost}}
#'
#' @export
fit <- function(datas=data, state=s, parms=p, odes=model,
                free=NULL, who=NULL, differ=NULL, fixed=NULL,
                tmin=0, tmax=NULL, ymin=NULL, ymax=NULL,
                log="", xlab="Time", ylab="Density",
                bootstrap=0, show=NULL, fun=NULL, costfun=cost,
                logpar=FALSE, lower=-Inf, upper=Inf,
                initial=FALSE, add=FALSE, timeplot=TRUE,
                legend=TRUE, main=NULL, sub=NULL, pchMap=NULL, ...) {
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_fit, args_run, args_plot)])
    if (length(unknown) > 0) warning(paste("Unknown argument(s):", unknown, sep=" "))
    dots_run <- dots[names(dots) %in% args_run]
    dots_fit <- dots[names(dots) %in% c(args_fit, args_run)]
    if ("method" %in% names(dots)) {
      if (!(dots$method %in% c(methods_run, methods_fit)))
        stop(paste("Unknown method:", dots$method))
      if (!(dots$method %in% methods_fit)) {
        dots_fit[["run_method"]] <- dots$method
        dots_fit[["method"]]     <- NULL
      }
      if (!(dots$method %in% methods_run))
        dots_run[["method"]] <- NULL
    }
  }
  if (is.null(free) & !is.null(who)) free <- who
  if (!is.null(fun)) fun <- match.fun(fun)
  if (is.data.frame(datas)) datas <- list(datas)
  nsets   <- length(datas)
  all     <- c(state, parms); allNms <- names(all)
  if (initial) totp <- parms else totp <- all
  isVar   <- setNames(c(rep(TRUE,length(state)), rep(FALSE,length(parms))), allNms)
  if (is.null(free) & is.null(differ)) free <- allNms
  if (initial) free <- idrop(names(state), free)
  ifree   <- index(free, names(totp))
  if (!is.null(fixed)) {
    if (!is.list(fixed)) stop("fixed should be a list")
    if (length(intersect(names(fixed), free))   > 0) stop("fixed should not overlap with free")
    if (length(intersect(names(fixed), differ)) > 0) stop("fixed should not overlap with differ")
  }
  if (is.null(differ)) {
    guess <- setNames(rep(0, length(free)), free)
  } else {
    if (!is.list(differ)) {
      ldiff <- makelist(differ, state=state, parms=parms, nsets=nsets)
    } else {
      ldiff <- differ; differ <- names(ldiff)
    }
    free  <- idrop(differ, free)
    guess <- setNames(rep(0, length(free)+nsets*length(differ)),
                      c(free, rep(differ, nsets)))
  }
  lenfree <- length(free)
  lendiff <- length(differ)
  VarsFree <- free[isVar[free]]
  ParsFree <- free[!isVar[free]]
  if (length(VarsFree) > 0) guess[VarsFree] <- state[VarsFree]
  if (length(ParsFree) > 0) guess[ParsFree] <- parms[ParsFree]
  if (!is.null(differ)) {
    for (inum in seq(lendiff))
      for (iset in seq(nsets))
        guess[lenfree+inum+(iset-1)*lendiff] <- ldiff[[differ[inum]]][iset]
    if (length(lower) == (lenfree+lendiff))
      lower <- c(lower, rep(lower[(lenfree+1):(lenfree+lendiff)], nsets-1))
    if (length(upper) == (lenfree+lendiff))
      upper <- c(upper, rep(upper[(lenfree+1):(lenfree+lendiff)], nsets-1))
  }
  if (logpar) {
    guess <- log(guess)
    if (length(lower) > 1) lower <- log(lower)
    if (length(upper) > 1) upper <- log(upper)
  }
  logy <- ifelse(grepl('y', log), TRUE, FALSE)
  f <- do.call('modFit', c(list(f=costfun, p=guess,
                                 datas=datas, odes=odes,
                                 state=state, parms=parms,
                                 free=free, differ=differ, fixed=fixed,
                                 fun=fun, logpar=logpar,
                                 ParsFree=ParsFree,
                                 lower=lower, upper=upper,
                                 initial=initial, isVar=isVar), dots_fit))
  found <- f$par
  if (logpar) found <- exp(found)
  cat("SSR:", f$ssr, " Estimates:\n")
  print(found)
  if (logpar) { cat("Log values free parameters:\n"); print(f$par) }
  if (timeplot) {
    tmaxn <- ifelse(is.null(tmax),
                    max(unlist(lapply(seq(nsets), function(i) max(datas[[i]][1])))), tmax)
    ymaxn <- ifelse(is.null(ymax),
                    max(unlist(lapply(seq(nsets), function(i) max(na.omit(datas[[i]][2:ncol(datas[[i]])]))))), ymax)
    yminn <- ifelse(is.null(ymin),
                    min(unlist(lapply(seq(nsets), function(i) min(na.omit(datas[[i]][2:ncol(datas[[i]])]))))), ymin)
    for (iset in seq(nsets)) {
      data <- datas[[iset]]
      if (length(VarsFree) > 0) state[VarsFree] <- found[VarsFree]
      if (length(ParsFree) > 0) parms[ParsFree] <- found[ParsFree]
      if (!is.null(fixed))
        for (inum in seq(length(fixed))) {
          name <- names(fixed)[inum]
          if (isVar[name]) state[match(name,names(state))] <- fixed[[inum]][iset]
          else             parms[match(name,names(parms))] <- fixed[[inum]][iset]
        }
      if (!is.null(differ))
        for (i in seq(lendiff)) {
          value <- found[lenfree+i+(iset-1)*lendiff]
          if (isVar[differ[i]]) state[match(differ[i],names(state))] <- value
          else                  parms[match(differ[i],names(parms))] <- value
        }
      if (initial) {
        if (data[1,1] > 0) stop("Data doesn't start at time=0")
        state[1:length(state)] <- unlist(data[1, 2:ncol(data)])
      }
      tmaxi <- ifelse(is.null(tmax), ifelse(add, tmaxn, max(data[,1])), tmax)
      nsol  <- do.call('run', c(list(tmax=tmaxi, state=state, parms=parms,
                                      odes=odes, table=TRUE, timeplot=FALSE), dots_run))
      ymaxi <- ifelse(is.null(ymax), ifelse(add, ymaxn,
                      max(na.omit(data[2:ncol(data)]), nsol[2:ncol(nsol)])), ymax)
      ymini <- ifelse(is.null(ymin), ifelse(add, yminn,
                      min(na.omit(data[2:ncol(data)]), nsol[2:ncol(nsol)])), ymin)
      solnames <- names(nsol)[2:ncol(nsol)]
      colnames <- names(data)[2:ncol(data)]
      imain <- main[min(length(main), iset)]
      isub  <- sub[min(length(sub),  iset)]
      if (is.null(show)) {
        timePlot(nsol, tmin=tmin, tmax=tmaxi, ymin=ymini, ymax=ymaxi, log=log,
                 main=imain, sub=isub, add=ifelse(iset>1, add, FALSE),
                 xlab=xlab, ylab=ylab, font.main=.fontMain(), font.sub=.fontSub(),
                 legend=legend)
        timePlot(data, draw=points, add=TRUE, legend=FALSE, lwd=1.5,
                 colMap=index(colnames, solnames), pchMap=pchMap)
      } else {
        for (i in show) {
          timePlot(nsol, tmin=tmin, tmax=tmaxi, ymin=ymini, ymax=ymaxi, log=log,
                   show=i, main=imain, sub=isub,
                   add=ifelse(i != show[1], add, FALSE),
                   xlab=xlab, ylab=ylab, font.main=.fontMain(), font.sub=.fontSub(),
                   legend=legend)
          if (i %in% colnames)
            timePlot(data, draw=points, add=TRUE, legend=FALSE, lwd=1.5,
                     show=i, colMap=index(colnames, solnames), pchMap=pchMap)
        }
      }
    }
  }
  if (bootstrap == 0) return(f)
  imat <- sapply(seq(bootstrap), function(i) {
    samples <- lapply(seq(nsets), function(j) {
      datas[[j]][sample(nrow(datas[[j]]), replace=TRUE), ]
    })
    ifit <- do.call('modFit', c(list(f=costfun, p=f$par,
                                      datas=samples, odes=odes,
                                      state=state, parms=parms,
                                      free=free, differ=differ, fixed=fixed,
                                      fun=fun, logpar=logpar,
                                      ParsFree=ParsFree,
                                      lower=lower, upper=upper,
                                      initial=initial, isVar=isVar), dots_fit))
    ifit$par
  })
  if (length(found) == 1)
    imat <- matrix(imat, nrow=1, dimnames=list(names(found[1]), NULL))
  print(apply(imat, 1, function(i)
    c(mean=mean(i), sd=sd(i), median=median(i), quantile(i, c(.025, .975)))))
  f$bootstrap <- t(imat)
  return(f)
}

# cost() -------------------------------------------------------------------

#' Cost function minimised by fit()
#'
#' \code{cost()} measures the mismatch between the model and the data for a
#' candidate parameter vector, returned as an \code{\link[FME]{modCost}} object.
#' It runs the model at the observation times and compares it with the data
#' across all datasets. This is the default \code{costfun} used by
#' \code{\link{fit}}; you rarely call it directly, but you can supply your own
#' via \code{fit(costfun = ...)}.
#'
#' @param datas List of observation data frames (as prepared by \code{fit}).
#' @param odes The model function \code{f(t, state, parms)}.
#' @param state Named numeric vector of initial state values.
#' @param parms Named numeric vector of parameters.
#' @param guess Named numeric vector of the parameter values currently being
#'   tried (supplied by \code{\link[FME]{modFit}}).
#' @param free Names of the shared free parameters/states.
#' @param differ Names of the per-dataset parameters/states.
#' @param fixed Named list of per-dataset fixed values.
#' @param fun Optional function applied to the data and model output before the
#'   cost is computed.
#' @param logpar If \code{TRUE}, \code{guess} is on a log scale.
#' @param ParsFree Names of the free parameters that are model parameters
#'   rather than state variables (set internally by \code{fit}).
#' @param initial If \code{TRUE}, take each dataset's initial state from its
#'   first row.
#' @param isVar Logical vector marking which names are state variables (set
#'   internally by \code{fit}).
#' @param ... Additional arguments passed on to \code{run} or
#'   \code{\link[FME]{modCost}}.
#'
#' @return An \code{\link[FME]{modCost}} object giving the residuals and cost
#'   used by \code{\link[FME]{modFit}}.
#'
#' @seealso \code{\link{fit}}
#'
#' @export
cost <- function(datas, odes, state, parms, guess, free, differ, fixed,
                 fun, logpar, ParsFree, initial, isVar, ...) {
  dots <- list(...)
  if (!is.null(dots)) {
    dots_run <- dots[names(dots) %in% args_run]
    dots_fit <- dots[names(dots) %in% args_fit]
    if ("run_method" %in% names(dots)) dots_run[["method"]] <- dots$run_method
  }
  if (!is.null(fun)) fun <- match.fun(fun)
  VarsFree <- free[isVar[free]]
  ParsFree <- free[!isVar[free]]
  lenfree  <- length(free)
  lendiff  <- length(differ)
  if (length(VarsFree) > 0) state[VarsFree] <- guess[VarsFree]
  if (length(ParsFree) > 0) parms[ParsFree] <- guess[ParsFree]
  nsets    <- length(datas)
  totcost  <- NULL
  for (iset in seq(nsets)) {
    data <- datas[[iset]]
    if (initial) {
      if (data[1,1] > 0) stop("Data doesn't start at time=0: data[", iset, "]")
      state[1:length(state)] <- unlist(data[1, 2:ncol(data)])
    }
    if (!is.null(fixed))
      for (inum in seq(length(fixed))) {
        name <- names(fixed)[inum]
        if (isVar[name]) state[match(name,names(state))] <- fixed[[inum]][iset]
        else             parms[match(name,names(parms))] <- fixed[[inum]][iset]
      }
    if (!is.null(differ))
      for (i in seq(lendiff)) {
        value <- guess[lenfree+i+(iset-1)*lendiff]
        if (isVar[differ[i]]) state[match(differ[i],names(state))] <- value
        else                  parms[match(differ[i],names(parms))] <- value
      }
    times <- sort(unique(data[,1]))
    if (!(0 %in% times)) times <- c(0, times)
    if (logpar) {
      if (iset == 1) parms[ParsFree]  <- exp(parms[ParsFree])
      parms[differ] <- exp(parms[differ])
    }
    nsol <- do.call('run', c(list(times=times, state=state, parms=parms,
                                   odes=odes, timeplot=FALSE, table=TRUE), dots_run))
    if (!is.null(fun)) {
      data[2:ncol(data)]  <- fun(data[2:ncol(data)])
      nsol[2:ncol(nsol)]  <- fun(nsol[2:ncol(nsol)])
    }
    totcost <- do.call('modCost', c(list(model=nsol, obs=data, cost=totcost), dots_fit))
  }
  return(totcost)
}

# timePlot() ---------------------------------------------------------------

#' Plot a time series
#'
#' \code{timePlot()} plots the columns of a data frame (time in the first
#' column, one or more variables in the rest) against time, colouring each
#' variable with the grindR palette. It is used internally by \code{\link{run}}
#' and \code{\link{fit}} to draw time courses, but can be called directly on any
#' such data frame, e.g. the result of \code{run(table = TRUE)}.
#'
#' @param data A data frame whose first column is time and whose remaining
#'   columns are the variables to plot.
#' @param tmin,tmax Time range (\code{tmax = NULL} uses the data range).
#' @param ymin,ymax Vertical range (\code{ymax = NULL} auto-scales).
#' @param log Which axes to draw on a log scale: \code{""}, \code{"x"},
#'   \code{"y"} or \code{"xy"}.
#' @param xlab,ylab Axis labels.
#' @param show Names or indices of the variables to plot (default: all).
#' @param legend If \code{TRUE}, draw a legend.
#' @param draw Drawing function, typically \code{\link[graphics]{lines}} or
#'   \code{\link[graphics]{points}}.
#' @param lwd Line width.
#' @param add If \code{TRUE}, add to the existing plot.
#' @param main,sub Plot title and subtitle.
#' @param colMap Optional mapping of columns to palette colour indices, used to
#'   keep colours consistent between model and data.
#' @param pchMap Optional mapping of columns to plotting characters.
#' @param ... Additional arguments passed on to \code{plot}.
#'
#' @return None; \code{timePlot} is called for the plot it draws.
#'
#' @examples
#' model <- function(t, state, parms) {
#'   with(as.list(c(state, parms)), list(c(r * N * (1 - N / K))))
#' }
#' p <- c(r = 0.5, K = 10)
#' s <- c(N = 0.1)
#'
#' out <- run(tmax = 20, table = TRUE, timeplot = FALSE)
#' timePlot(out)
#'
#' @seealso \code{\link{run}}, \code{\link{fit}}
#'
#' @export
timePlot <- function(data, tmin=0, tmax=NULL, ymin=0, ymax=NULL,
                     log="", xlab="Time", ylab="Density",
                     show=NULL, legend=TRUE, draw=lines, lwd=2,
                     add=FALSE, main=NULL, sub=NULL,
                     colMap=NULL, pchMap=NULL, ...) {
  colnames <- names(data)[2:ncol(data)]
  ivar <- seq(ncol(data)-1)
  if (!is.null(show)) {
    ivar <- sort(index(show, colnames))
    data <- data[, c(1, ivar+1)]
  }
  if (!is.null(draw)) draw <- match.fun(draw)
  if (is.null(tmax)) tmax <- max(data[,1])
  if (is.null(ymax)) ymax <- max(na.omit(data[, 2:ncol(data)]))
  logx <- ifelse(grepl('x', log), TRUE, FALSE)
  logy <- ifelse(grepl('y', log), TRUE, FALSE)
  if (tmin == 0 & logx) tmin <- min(data[,1])
  if (ymin == 0 & logy) {
    ymin <- min(na.omit(data[, 2:ncol(data)]))
    if (ymin <= 0) ymin <- min(1, ymax/100)
  }
  if (!add)
    plot(1, 1, type='n', xlim=c(tmin,tmax), ylim=c(ymin,ymax),
         log=log, xlab=xlab, ylab=ylab, main=main, sub=sub,
         font.main=.fontMain(), font.sub=.fontSub(), ...)
  for (i in seq(ncol(data)-1)) {
    j <- ifelse(is.null(colMap), ivar[i], colMap[ivar[i]])
    k <- ifelse(is.null(pchMap), j, pchMap[ivar[i]])
    draw(data[,1], data[,i+1], col=.pal()[min(j,length(.pal()))], lwd=lwd, pch=k)
  }
  if (legend)
    legend("topright", legend=colnames[ivar], col=.pal()[ivar],
           lty=1, lwd=lwd, cex=.sizeLegend(),
           pch=ifelse(identical(draw, lines), NA, ivar))
}

# dummyEvent() -------------------------------------------------------------

#' A no-op event function
#'
#' Returns the state unchanged. Used internally by \code{\link{run}} to force
#' the solver to stop exactly at specified times (the \code{arrest} argument)
#' without altering the state.
#'
#' @param t Current time.
#' @param state Named numeric vector of state values.
#' @param parms Named numeric vector of parameters.
#' @return The \code{state}, unchanged.
#' @keywords internal
dummyEvent <- function(t, state, parms) return(state)

# Helper functions ---------------------------------------------------------

#' Match names to their positions
#'
#' Returns the positions of \code{strings} within \code{names}. Used internally
#' to turn variable or parameter names into indices.
#'
#' @param strings Character vector of names to look up.
#' @param names Character vector of reference names.
#' @param error If \code{TRUE}, stop when a string is not found; otherwise
#'   ignore unmatched strings.
#' @return Integer positions of the matched strings, or \code{NULL} if none
#'   match.
#' @keywords internal
index <- function(strings, names, error=TRUE) {
  hit <- strings %in% names
  if (error & length(strings[!hit] > 0)) stop("Unknown: ", paste(strings[!hit], collapse=", "))
  m <- match(strings[hit], names)
  if (length(m) > 0) return(m)
  return(NULL)
}

#' Drop names from a set
#'
#' Returns \code{names} with any elements listed in \code{strings} removed.
#'
#' @param strings Character vector of names to remove.
#' @param names Character vector to remove them from.
#' @return The remaining names, or \code{NULL} if all were removed.
#' @keywords internal
idrop <- function(strings, names) {
  hit <- strings %in% names
  m   <- match(strings[hit], names)
  if (length(m) == length(names)) return(NULL)
  if (length(m) > 0) return(names[-m])
  return(names)
}

#' Build a list of repeated starting values
#'
#' For the given state/parameter names, builds a list in which each value is
#' repeated \code{nsets} times. Used by \code{\link{fit}} to set up per-dataset
#' (\code{differ}) starting values.
#'
#' @param strings Character vector of state/parameter names.
#' @param state Named numeric vector of state values. Defaults to the global
#'   \code{s}.
#' @param parms Named numeric vector of parameters. Defaults to the global
#'   \code{p}.
#' @param nsets Number of datasets (times to repeat each value).
#' @return A named list of starting values, one element per name in
#'   \code{strings}.
#' @keywords internal
makelist <- function(strings, state=s, parms=p, nsets=1) {
  all <- c(state, parms)
  nms <- names(all)
  hit <- strings %in% nms
  if (length(strings[!hit] > 0))
    stop("Unknown: ", paste(c(strings[!hit]), collapse=", "))
  lst <- as.list(all[match(strings[hit], nms)])
  return(lapply(lst, rep, lst, nsets))
}

# Argument-routing constants -----------------------------------------------
# Defined here, after all functions, because they call formals(run) etc.
# These are read-only, so keeping them as namespace objects is fine.

args_plot    <- unique(c(names(c(formals(graphics::plot.default),
                                  formals(graphics::axis),
                                  formals(graphics::axTicks))),
                          "xaxp","yaxp"))
args_fit     <- unique(names(c(formals(FME::modFit), formals(FME::modCost))))
args_run     <- unique(names(c(formals(run), formals(deSolve::ode), formals(deSolve::lsoda))))
args_run_dde <- unique(names(c(formals(run), formals(deSolve::ode),
                                formals(deSolve::dede), formals(deSolve::lsoda))))
args_steady  <- unique(names(c(formals(rootSolve::steady), formals(rootSolve::stode))))
methods_run  <- as.character(formals(deSolve::ode)$method)
methods_fit  <- as.character(formals(FME::modFit)$method)

# Startup message ----------------------------------------------------------
# No global-environment seeding: user-facing options (colors, fonts, legend
# size) are exported package data read through getters, and internal state
# lives in pkg_env. See the top of this file.

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("grindR (", grind_version, ") loaded")
}
