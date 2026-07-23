# Phase plane analysis in R: grind.R 
# Developed by: Rob de Boer, Utrecht University
# Maintained by: Bram van Dijk, Utrecht University
grind_version <- "23-07-2026"

library(deSolve)     # run() calls the ode() function
library(rootSolve)   # newton() and continue() call steady()
library(coda)        # required by FME
library(FME)         # fit() calls modFit() and modCost()

options(stringsAsFactors=FALSE)
colors <- c("red","blue","darkgreen","darkorange","darkmagenta","gold","darkorchid","aquamarine","deeppink","gray",seq(2,991))
#colors <- seq(2,101)  # Use standard R colors
ncolors <- length(colors)
sizeLegend <- 0.75    # legend size is 75% of the default value in R
font.main <- 1        # plain (default is bold: 2)
font.sub  <- 1        # plain (default is bold: 2)

x_plane <- 1; xmin_plane <- -0.001; xmax_plane <- 1.05
y_plane <- 2; ymin_plane <- -0.001; ymax_plane <- 1.05
log_plane <- ""; addone_plane <- FALSE

plane_coord <- function(Log,Min,Max,npixels) {
  if (Log) return(10^seq(log10(Min),log10(Max),length.out=npixels))
  return(seq(Min,Max,length.out=npixels))
}

plane <- function(xmin=-0.001, xmax=1.05, ymin=-0.001, ymax=1.05, xlab="", ylab="", log="", npixels=500, state=s, parms=p, odes=model, x=1, y=2, time=0, grid=5, eps=NULL, show=NULL, addone=FALSE, portrait=FALSE, vector=FALSE, add=FALSE, legend=TRUE, zero=TRUE, lwd=2, col="black", pch=20, vectorstyle="arrows", vectorlen=1, arrowsize=0.3, ...) {
  # Make a phase plane with nullclines and/or phase portrait
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
  }
  dots_run <- if (!is.null(dots)) dots[names(dots) %in% args_run] else NULL
  if (!is.null(eps)) print("Option eps is no longer supported")
  if (add) {
    x <- x_plane
    y <- y_plane
    xmin <- xmin_plane; xmax <- xmax_plane
    ymin <- ymin_plane; ymax <- ymax_plane
    log <- log_plane; addone <- addone_plane
  } else {
    if (!is.numeric(x)) x <- index(x,names(state))
    if (!is.numeric(y)) y <- index(y,names(state))
    x_plane <<- x
    y_plane <<- y
    xmin_plane <<- xmin; xmax_plane <<- xmax
    ymin_plane <<- ymin; ymax_plane <<- ymax 
    log_plane <<- log; addone_plane <<- addone
  }
  ishows <- if (!is.null(show)) index(show, names(state)) else c(x, y)
  nvar <- length(state)
  if (zero) state[1:nvar] <- rep(0,nvar)
  lvec <- 25                         # length of vector
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
    do.call('plot',c(list(1,1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",log=log,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
    mtext(xlab, side = 1, line = 2.5, col = colors[x], font = 2, cex = 1)
    mtext(ylab, side = 2, line = 2.5, col = colors[y], font = 2, cex = 1)
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
    contour(xc,yc,dvstate[,,i],levels=0,drawlabels=FALSE,add=TRUE,col=colors[i],lwd=lwd)
  
  if (portrait | vector) {
    
    pin <- par("pin") # get the dimension of the plot in R in inches

    dx <- if (logx) (log10(xmax)-log10(xmin))/grid else (xmax-xmin)/grid
    dy <- if (logy) (log10(ymax)-log10(ymin))/grid else (ymax-ymin)/grid

    # Shaft length tracks the grid-cell size but is floored at min(pin)/lvec, so
    # vectors stay proportional to the field (robust to par(mfrow) and coarse
    # grids) without vanishing at fine grids. Heads scale with the shaft.
    shaft_in <- min(pin) * max(1/lvec, 0.5/grid) * vectorlen
    ahead    <- max(0.02, min(0.15, shaft_in * arrowsize))
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
            
          } else if (vector == 1) {
            # Two orthogonal arrows, sign-only, with arrowheads
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
            
          } else if (vector == 2) {
            # Single arrow in true (dx, dy) direction, normalised length
            mag <- sqrt(dvx^2 + dvy^2)
            if (mag > 0) {
              norm_dx <- dvx / mag
              norm_dy <- dvy / mag
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
          }
          else{
            stop("Unknown vectorstyle: ", vectorstyle)
          }
        }
      }
    }
  }
  if(!add){
    if (legend)
      legend("topright",legend=names(state)[ishows],col=colors[ishows],lty=1,lwd=lwd,cex=sizeLegend)
  }
}

dummyEvent <-function(t, state , parms) return(state)

run <- function(tmax=100, tstep=1, state=s, parms=p, odes=model, ymin=0, ymax=NULL, log="", xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL, show=NULL, arrest=NULL, events=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE, table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2, col="black", pch=20, ...) {   
  # run model and make a table, time plot, or trajectory
  if (delay & (solution | !is.null(after))) stop("Don't use solution or after with delay equations")
  if (delay) args_run <- args_run_dde
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
    dots_run <- dots[names(dots) %in% args_run]
  } else dots_run <- NULL
  nvar <- length(state)
  #if (!is.numeric(x)) x <- index(x,names(state))
  #if (!is.numeric(y)) y <- index(y,names(state))
  #if (!is.null(show)) ishows <- index(show, names(state))
  #else ishows <- seq(1,nvar)
  if (is.null(times)) {
    times <- seq(tmin,tmax,by=tstep)
  } else { 
    times <- sort(times)
    tmin <- min(times)
    tmax <- max(times)
  }
  if (!is.null(arrest)) {
    if (!is.null(events)) stop("Don't combine the option arrest with events")
    if (!is.numeric(arrest)) arrest <- sort(as.numeric(parms[arrest]))
    nearby <- nearestEvent(arrest,times)  # Find nearby times
    nearby <- nearby[nearby < arrest]     # Pick those before arrest
    lennear <- length(nearby)
    if (lennear == 1 && nearby[1] == 0) lennear <- 0
    if (lennear > 0) {                    # Add to arrest for safety
      if (nearby[1] == 0) nearby <- nearby[2:lennear]
      arrest <- sort(unique(c(nearby,arrest)))   
    }
    events <- list(func=dummyEvent,time=arrest)
    times <- cleanEventTimes(times,arrest,eps = .Machine$double.eps*10)
    #times <- cleanEventTimes(times,arrest,eps=min(arrest)/10000)
    times <- sort(c(times,arrest))
  }
  if (solution) {                              # Run in one go
    if (!is.null(after)) stop("Don't combine the option after with solution")
    nsol <- sapply(times,odes,state,parms)
    if (is.list(nsol)) {
      nsol <- unlist(nsol)
      if (nvar > 1) dim(nsol) <- c(nvar,length(times))
    }
    if (nvar > 1) nsol <- data.frame(times,t(nsol))
    else nsol <- data.frame(times,nsol)
    names(nsol) <- c("time",names(state))
  } else {
    if (is.null(after)){                      # Run in one go
      nsol <- as.data.frame(
        do.call(if(!delay) 'ode' else 'dede', 
                c(list(times = times, func = odes, y = state, parms = parms, events = events), dots_run))
      )
    } else {                                  # After: Run in individual steps
      keep <- state
      nsol <- t(sapply(seq(length(times)-1),function(i){
        t <- times[i+1]
        f <- do.call('ode',c(list(times=c(times[i],t),func=odes,y=state,parms=parms),dots_run))
        dim(f) <- c(2,nvar+1)
        state[1:nvar] <- f[2,2:(nvar+1)]
        eval(parse(text=after))
        parms <<- parms
        state <<- state
      }
      ))
      if (nvar > 1) {
        nsol <- as.data.frame(cbind(times,rbind(as.numeric(keep),nsol)))
      } else { 
        nsol <- as.data.frame(cbind(times,c(keep,nsol)))
      }
      names(nsol) <- c("time",names(state))
      state <- keep
    }
  }
  if (!is.null(tweak)) eval(parse(text=tweak))
  if (timeplot & !traject)
    do.call('timePlot',c(list(data=nsol,tmin=tmin,tmax=tmax,ymin=ymin,ymax=ymax,log=log,add=add,xlab=xlab,ylab=ylab,show=show,draw=draw,lwd=lwd,legend=legend,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
  if (traject) {
    points(nsol[1,x_plane+1],nsol[1,y_plane+1],pch=pch)
    lines(nsol[,x_plane+1],nsol[,y_plane+1],lwd=lwd,col=col)
  }
  if (table) return(nsol)
  f <- state
  f[1:length(f)] <- as.numeric(nsol[nrow(nsol),2:(nvar+1)])
  return(f)
}

newton <- function(state=s, parms=p, odes=model, time=0, positive=FALSE, jacobian=FALSE, vector=FALSE, plot=FALSE, silent=FALSE, addone=FALSE, ...) {  
  # find a steady state
  # if (!is.numeric(x)) x <- index(x,names(state))
  # if (!is.numeric(y)) y <- index(y,names(state))
  One <- ifelse(addone, 1, 0)
  q <- steady(y=state,func=odes,parms=parms,time=time,positive=positive, ...)
  if (attr(q,"steady")) {
    equ <- q$y
    equ <- ifelse(abs(equ) < 1e-8, 0, equ)
    jac <- jacobian.full(y=equ,func=odes,parms=parms)
    eig <- eigen(jac)
    dom <- max(Re(eig$values))
    if (!silent) {
      print(equ)
      if (length(equ) == 2) {
        if (is.complex(eig$values[1])) {
          cat(ifelse(dom < 0, "Stable spiral point, ", "Unstable spiral point, "))
        } else if(prod(eig$values) > 0) {
          cat(ifelse(dom < 0, "Stable node, ", "Unstable node, "))
        } else cat("Saddle point (unstable), ")
      } else {
        cat(ifelse(dom < 0, "Stable point, ", "Unstable point, "))
      }
      cat("eigenvalues:\n")
      print(eig$values)
    }
    if (vector) {cat("Eigenvectors:\n"); print(eig$vectors)}
    if (jacobian) {cat("Jacobian:\n"); print(jac)}
    if (plot) 
      points(equ[x_plane] + One, equ[y_plane] + One, pch = ifelse(dom < 0, 19, 1))
    if (silent) return(list(state=equ,jacobian=jac,values=eig$values,vectors=eig$vectors))
    return(equ)
  }
  cat("No convergence: start closer to a steady state")
  return(NULL)
}

x_continue <- 1; xmin_continue <- 0; xmax_continue <- 1
y_continue <- 1; ymin_continue <- 0; ymax_continue <- 1
log_continue <- ""; addone_continue <- FALSE

continue <- function(state=s, parms=p, odes=model, step=0.01, x=1, y=2, time=0, xmin=0, xmax=1,ymin=0, ymax=1.05, xlab="", ylab="", log="", col=c("red","black","blue"), lwd= c(2,1,1), addone=FALSE, positive=FALSE, nvar=FALSE, add=FALSE, ...) {  
  # continue a steady state
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_steady,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
    dots_steady <- dots[names(dots) %in% args_steady]
  }else dots_steady <- NULL
  if (add) {
    x <- x_continue
    y <- y_continue
  } else {
    if (!is.numeric(x)) x <- index(x,names(parms))
    if (!is.numeric(y)) y <- index(y,names(state))
    x_continue <<- x
    y_continue <<- y
  }
  
  p0 <- parms[x]
  q0 <- do.call('steady',c(list(y=state,func=odes,parms=parms,time=time,positive=positive),dots_steady))
  if (!attr(q0,"steady"))
    stop("No convergence: start closer to a steady state")
  cat("Starting at",names(parms[x]),"=",parms[x],"with:\n")
  print(q0$y)
  bary <- q0$y[y]
  if (!add) {
    if (missing(xmax) & parms[x] >= 1) xmax <- 2*parms[x]
    if (missing(xmin) & parms[x] < 0) xmin <- 2*parms[x]
    if (!missing(xmin) & xmin >= parms[x]) stop("xmin should be smaller than parameter")
    if (!missing(xmax) & xmax <= parms[x]) stop("xmax should be larger than parameter")
    if (missing(ymax) & bary >= 1.05) ymax <- 2*bary
    if (missing(ymin) & bary < 0) ymin <- 2*bary
    if (!missing(ymin) & ymin >= bary & !addone) stop("ymin should be smaller than y-variable")
    if (!missing(ymax) & ymax <= bary) stop("ymax should be larger than y-variable")
    if (xlab == "") xlab <- names(p0)
    if (ylab == "") {
      ylab <- names(state)[y]
      if (addone) ylab <- paste(ylab,"+ 1")
    }
    do.call('plot',c(list(1,1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,log=log,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
    xmin_continue <<- xmin; xmax_continue <<- xmax
    ymin_continue <<- ymin; ymax_continue <<- ymax 
    log_continue <<- log; addone_continue <<- addone
  } else {
    xmin <- xmin_continue; xmax <- xmax_continue
    ymin <- ymin_continue; ymax <- ymax_continue
    log <- log_continue; addone <- addone_continue
    if (ymin >= bary & !addone) stop("Initial point below minimum of current y-axis")
    if (ymax <= bary) stop("Initial point above maximum of current y-axis")
  }
  logx <- ifelse(grepl('x',log), TRUE, FALSE)
  COL <- function(s,i) {
    if (!nvar) return(col[i])
    return(col[length(s[s>1e-9])])
  }
  FUN <- function(lastState,lastDom,step) {
    lastP <- p0
    preLastState <- lastState
    nok <- 0
    while (xmin < lastP & lastP < xmax & ymin < lastState[y]+One & lastState[y] < ymax) {
      parms[x] <- ifelse(logx, lastP*(1 + step), lastP + step)
      q <- do.call('steady',c(list(y=lastState,func=odes,parms=parms,time=time,positive=positive),dots_steady))
      newState <- q$y  # should be steady state and closeby
      if (attr(q,"steady") & sum(abs(newState-lastState))/(1e-9+sum(abs(lastState))) < 0.1) {
        jac <- jacobian.full(y=newState,func=odes,parms=parms)
        dom <- sign(max(Re(eigen(jac)$values)))
        if (dom != lastDom) cat("Bifurcation at",names(parms[x]),"=",parms[x],"\n")
        lines(c(if(logx) parms[x]/(1 + step) else parms[x] - step, parms[x]),
              c(lastState[y] + One, newState[y] + One),
              col = COL(lastState, dom + 2),
              lwd = lwd[dom + 2])
        preLastState <- lastState
        lastState <- newState
        lastDom <- dom
        lastP <- parms[x]
        if (nok > 10 & abs(step) < actualStep) step <- sign(step)*min(2*abs(step),actualStep)
        nok <- nok + 1
      }else{
        nok <- 0
        if (abs(step) > actualStep/100) {
          step <- step/2
        } else { # Go back one step, overpredict, call steady, and turn
          parms[x] <- lastP
          predState <- lastState + 5*(lastState-preLastState)
          q <- do.call('steady',c(list(y=predState,func=odes,parms=parms,time=time,positive=positive),dots_steady))
          newState <- q$y  # should be steady state and not the same
          if (attr(q,"steady") & sum(abs(newState-lastState))/(1e-9+sum(abs(lastState))) > 0.001) {
            cat("Turning point point at",names(parms[x]),"=",parms[x],"\n")
            jac <- jacobian.full(y=newState,func=odes,parms=parms)
            dom <- sign(max(Re(eigen(jac)$values)))
            middle <- (lastState[y]+newState[y])/2
            lines(c(parms[x],parms[x]),c(lastState[y]+One,middle+One), col=COL(lastState,lastDom+2),lwd=lwd[lastDom+2])
            lines(c(parms[x],parms[x]),c(middle+One,newState[y]+One), col=COL(newState,dom+2),lwd=lwd[dom+2])
            step <- -step
            preLastState <- lastState
            lastState <- newState
            lastDom <- dom
            lastP <- parms[x]
          }else{
            cat("Final point at",names(parms[x]),"=",parms[x],"\n")
            cat("If this looks wrong try changing the step size\n")
            break
          }
        }
      }
    }
  }
  One <- ifelse(addone, 1, 0)
  orgWarn <- getOption("warn")
  options(warn = -1)
  jac <- jacobian.full(y=q0$y,func=odes,parms=parms)
  dom <- sign(max(Re(eigen(jac)$values)))
  actualStep <- if(logx) step else step*xmax
  FUN(lastState=q0$y,lastDom=dom,actualStep)
  FUN(lastState=q0$y,lastDom=dom,-actualStep)
  options(warn = orgWarn)
  return(NULL)
}

fit <- function(datas=data, state=s, parms=p, odes=model, free=NULL, who=NULL, differ=NULL, fixed=NULL, tmin=0, tmax=NULL, ymin=NULL, ymax=NULL, log="", xlab="Time", ylab="Density", bootstrap=0, show=NULL, fun=NULL, costfun=cost, logpar=FALSE, lower=-Inf, upper=Inf, initial=FALSE, add=FALSE, timeplot=TRUE, legend=TRUE, main=NULL, sub=NULL, pchMap=NULL, ...) {
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_fit,args_run,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
    dots_run <- dots[names(dots) %in% args_run]
    dots_fit <- dots[names(dots) %in% c(args_fit,args_run)]
    if ("method" %in% names(dots)) {
      if (!(dots$method %in% c(methods_run,methods_fit))) 
        stop(paste("Unknown method:",dots$method))
      if (!(dots$method %in% methods_fit)) {
        dots_fit[["run_method"]] <- dots$method  # Used in cost
        dots_fit[["method"]] <- NULL
      }
      if (!(dots$method %in% methods_run)) 
        dots_run[["method"]] <- NULL
    }
  }
  if (is.null(free) & !is.null(who)) free <- who  # for compatability
  if (!is.null(fun)) fun <- match.fun(fun)
  if (is.data.frame(datas)) datas <- list(datas)
  nsets <- length(datas)
  all <- c(state,parms);  allNms <- names(all)
  if (initial) totp <- parms else totp <- all
  isVar <- setNames(c(rep(TRUE,length(state)),rep(FALSE,length(parms))),allNms) 
  if (is.null(free) & is.null(differ)) free <- allNms               
  if (initial) free <- idrop(names(state),free)   # remove state from free
  ifree <- index(free,names(totp))                # Check if free is correct
  if (!is.null(fixed)) {
    if (!is.list(fixed)) stop("fixed should be a list")
    if (length(intersect(names(fixed),free)) > 0) stop("fixed should not overlap with free")
    if (length(intersect(names(fixed),differ)) > 0) stop("fixed should not overlap with differ")
  }
  if (is.null(differ)) { 
    guess <- setNames(rep(0,length(free)),free)
  } else {
    if (!is.list(differ)) {
      ldiff <- makelist(differ,state=state,parms=parms,nsets=nsets)
    } else {
      ldiff <- differ; differ <- names(ldiff)
    }
    free <- idrop(differ,free)                      # remove differ from free
    guess <- setNames(rep(0,length(free)+nsets*length(differ)), c(free,rep(differ,nsets)))
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
      lower <- c(lower,rep(lower[(lenfree+1):(lenfree+lendiff)],nsets-1))
    if (length(upper) == (lenfree+lendiff))
      upper <- c(upper,rep(upper[(lenfree+1):(lenfree+lendiff)],nsets-1))
  }
  if (logpar) {
    guess <- log(guess)
    if (length(lower) > 1) lower <- log(lower)
    if (length(upper) > 1) upper <- log(upper)
  }
  logy <- ifelse(grepl('y',log), TRUE, FALSE)
  f <- do.call('modFit',c(list(f=costfun,p=guess,datas=datas,odes=odes,state=state,parms=parms,free=free,differ=differ,fixed=fixed,fun=fun,logpar=logpar,ParsFree=ParsFree,lower=lower,upper=upper,initial=initial,isVar=isVar),dots_fit))
  found <- f$par
  if (logpar) found <- exp(found)
  cat("SSR:",f$ssr," Estimates:\n")
  print(found)
  if (logpar) {
    cat("Log values free parameters:\n")
    print(f$par)
  }
  if (timeplot) {         # We are done, now plot the results 
    tmaxn <- ifelse(is.null(tmax),max(unlist(lapply(seq(nsets),function(i){max(datas[[i]][1])}))),tmax)
    ymaxn <- ifelse(is.null(ymax),max(unlist(lapply(seq(nsets),function(i){max(na.omit(datas[[i]][2:ncol(datas[[i]])]))}))),ymax)
    yminn <- ifelse(is.null(ymin),min(unlist(lapply(seq(nsets),function(i){min(na.omit(datas[[i]][2:ncol(datas[[i]])]))}))),ymin)
    for (iset in seq(nsets)) {
      data <- datas[[iset]]
      if (length(VarsFree) > 0) state[VarsFree] <- found[VarsFree]
      if (length(ParsFree) > 0) parms[ParsFree] <- found[ParsFree]
      
      if (!is.null(fixed)) 
        for (inum in seq(length(fixed))) {
          name <- names(fixed)[inum]
          if (isVar[name]) {
            state[match(name,names(state))] <- fixed[[inum]][iset] 
          } else {
            parms[match(name,names(parms))] <- fixed[[inum]][iset]
          }
        }
      if (!is.null(differ)) 
        for (i in seq(lendiff)) { # Copy found[iset] into state and parms
          value <- found[lenfree+i+(iset-1)*lendiff]
          if (isVar[differ[i]]) {
            state[match(differ[i],names(state))] <- value 
          } else {
            parms[match(differ[i],names(parms))] <- value
          }
        }
      if (initial) {
        if (data[1,1] > 0) stop("Data doesn't start at time=0")
        state[1:length(state)] <- unlist(data[1,2:ncol(data)])
      }
      tmaxi <- ifelse(is.null(tmax),ifelse(add,tmaxn,max(data[,1])),tmax)
      nsol <- do.call('run',c(list(tmax=tmaxi,state=state,parms=parms,odes=odes,table=TRUE,timeplot=FALSE),dots_run))
      ymaxi <- ifelse(is.null(ymax),ifelse(add,ymaxn,max(na.omit(data[2:ncol(data)]),nsol[2:ncol(nsol)])),ymax)
      ymini <- ifelse(is.null(ymin),ifelse(add,yminn,min(na.omit(data[2:ncol(data)]),nsol[2:ncol(nsol)])),ymin)
      solnames <- names(nsol)[2:ncol(nsol)]
      colnames <- names(data)[2:ncol(data)]
      imain <- main[min(length(main),iset)]
      isub <- sub[min(length(sub),iset)]    
      if (is.null(show)) {
        timePlot(nsol,tmin=tmin,tmax=tmaxi,ymin=ymini,ymax=ymaxi,log=log,main=imain,sub=isub,add=ifelse(iset>1,add,FALSE),xlab=xlab,ylab=ylab,font.main=font.main,font.sub=font.sub,legend=legend)
        timePlot(data,draw=points,add=TRUE,legend=FALSE,lwd=1.5,colMap=index(colnames,solnames),pchMap=pchMap)
      }else{
        for (i in show) {
          timePlot(nsol,tmin=tmin,tmax=tmaxi,ymin=ymini,ymax=ymaxi,log=log,show=i,main=imain,sub=isub,add=ifelse(i != show[1],add,FALSE),xlab=xlab,ylab=ylab,font.main=font.main,font.sub=font.sub,legend=legend)
          if (i %in% colnames)
            timePlot(data,draw=points,add=TRUE,legend=FALSE,lwd=1.5,show=i,colMap=index(colnames,solnames),pchMap=pchMap)
        }
      }
    }}
  if (bootstrap == 0) return(f)
  imat <- sapply(seq(bootstrap), function(i) {
    samples <- lapply(seq(nsets),function(j){datas[[j]][sample(nrow(datas[[j]]),replace=TRUE),]})
    ifit <- do.call('modFit',c(list(f=costfun,p=f$par,datas=samples,odes=odes,state=state,parms=parms,free=free,differ=differ,fixed=fixed,fun=fun,logpar=logpar,ParsFree=ParsFree,lower=lower,upper=upper,initial=initial,isVar=isVar),dots_fit))
    ifit$par
  })
  if (length(found) == 1) 
    imat <- matrix(imat,nrow=1,dimnames=list(names(found[1]),NULL))
  print(apply(imat, 1, function(i)c(mean=mean(i),sd=sd(i),median=median(i),quantile(i,c(.025, .975)))))
  f$bootstrap <- t(imat)
  return(f)
}

cost <- function(datas, odes, state, parms, guess, free, differ, fixed, fun, logpar, ParsFree, initial, isVar, ...) {
  dots <- list(...)
  if (!is.null(dots)) {
    dots_run <- dots[names(dots) %in% args_run]
    dots_fit <- dots[names(dots) %in% args_fit]
    if ("run_method" %in% names(dots)) dots_run[["method"]] <- dots$run_method
  }
  if (!is.null(fun)) fun <- match.fun(fun)
  VarsFree <- free[isVar[free]]
  ParsFree <- free[!isVar[free]]
  lenfree <- length(free)
  lendiff <- length(differ)
  if (length(VarsFree) > 0) state[VarsFree] <- guess[VarsFree]
  if (length(ParsFree) > 0) parms[ParsFree] <- guess[ParsFree]
  nsets <- length(datas)
  totcost <- NULL
  for (iset in seq(nsets)) {
    data <- datas[[iset]]
    if (initial) {
      if (data[1,1] > 0) stop("Data doesn't start at time=0: data[",iset,"]")
      state[1:length(state)] <- unlist(data[1,2:ncol(data)])
    }
    if (!is.null(fixed)) 
      for (inum in seq(length(fixed))) {
        name <- names(fixed)[inum]
        if (isVar[name]) {
          state[match(name,names(state))] <- fixed[[inum]][iset] 
        } else {
          parms[match(name,names(parms))] <- fixed[[inum]][iset]
        }
      }
    if (!is.null(differ)) for (i in seq(lendiff)) {
      value <- guess[lenfree+i+(iset-1)*lendiff]
      if (isVar[differ[i]]) {
        state[match(differ[i],names(state))] <- value 
      } else {
        parms[match(differ[i],names(parms))] <- value
      }
    }
    times <- sort(unique(data[,1]))
    if (!(0 %in% times)) times <- c(0,times)
    if (logpar) {
      if (iset == 1)  parms[ParsFree] <- exp(parms[ParsFree])
      parms[differ] <- exp(parms[differ])
    }
    nsol <- do.call('run',c(list(times=times,state=state,parms=parms,odes=odes,timeplot=FALSE,table=TRUE),dots_run))
    if (!is.null(fun)) {
      data[2:ncol(data)]=fun(data[2:ncol(data)])
      nsol[2:ncol(nsol)]=fun(nsol[2:ncol(nsol)])
    }
    totcost <- do.call('modCost',c(list(model=nsol,obs=data,cost=totcost),dots_fit)) 
  }
  return (totcost)
}

timePlot <- function(data, tmin=0, tmax=NULL, ymin=0, ymax=NULL, log="", xlab="Time", ylab="Density", show=NULL, legend=TRUE, draw=lines, lwd=2, add=FALSE, main=NULL, sub=NULL, colMap=NULL, pchMap=NULL, ...) {
  colnames <- names(data)[2:ncol(data)]
  ivar <- seq(ncol(data)-1)
  if (!is.null(show)) {
    ivar <- sort(index(show,colnames))
    data <- data[,c(1,ivar+1)]
  }
  if (!is.null(draw)) draw <- match.fun(draw)
  if (is.null(tmax)) tmax <- max(data[,1])
  if (is.null(ymax)) ymax <- max(na.omit(data[,2:ncol(data)]))
  logx <- ifelse(grepl('x',log), TRUE, FALSE)
  logy <- ifelse(grepl('y',log), TRUE, FALSE)
  if (tmin == 0 & logx) tmin <- min(data[,1])
  if (ymin == 0 & logy) {
    ymin <- min(na.omit(data[,2:ncol(data)]))
    if (ymin <= 0) ymin <- min(1,ymax/100)
  }
  if (!add)
    plot(1,1,type='n',xlim=c(tmin,tmax),ylim=c(ymin,ymax),log=log,xlab=xlab,ylab=ylab,main=main,sub=sub,font.main=font.main,font.sub=font.sub,...)
  
  for (i in seq(ncol(data)-1)) {
    #if (logy & min(data[,i+1]) <= 0) next
    j <- ifelse(is.null(colMap),ivar[i],colMap[ivar[i]])
    k <- ifelse(is.null(pchMap),j,pchMap[ivar[i]])
    draw(data[,1],data[,i+1],col=colors[min(j,ncolors)],lwd=lwd,pch=k)
  }
  if (legend) legend("topright", legend = colnames[ivar],
                     col=colors[ivar], lty=1, lwd=lwd, cex=sizeLegend,
                     pch = ifelse(identical(draw, lines), NA, ivar))
}

index <- function(strings, names, error=TRUE) {   # Return indices of strings in names
  hit <- strings %in% names
  if (error & length(strings[!hit] > 0)) stop("Unknown: ", paste(strings[!hit],collapse=", "))
  m <- match(strings[hit], names)
  if (length(m) > 0) return(m)
  return(NULL)
}

idrop <- function(strings, names) {            # Remove all strings from names
  hit <- strings %in% names
  m <- match(strings[hit], names)
  if (length(m) == length(names)) return(NULL)
  if (length(m) > 0) return(names[-m])
  return(names)
}

makelist <- function(strings,state=s,parms=p,nsets=1) { # Make a list with nsets default values
  all <- c(state,parms)
  nms <- names(all)
  hit <- strings %in% nms
  if (length(strings[!hit] > 0))
    stop("Unknown: ", paste(c(strings[!hit]),collapse=", "))
  lst <- as.list(all[match(strings[hit], nms)])
  return(lapply(lst,rep,lst,nsets))
}

args_plot <- unique(c(names(c(formals(plot.default),formals(axis),formals(axisTicks))),"xaxp","yaxp"))
args_fit <- unique(names(c(formals(modFit),formals(modCost))))
args_run <- unique(names(c(formals(run),formals(ode),formals(lsoda))))
args_run_dde <- unique(names(c(formals(run),formals(ode),formals(dede),formals(lsoda))))
args_steady <- unique(names(c(formals(steady),formals(stode))))
methods_run <- as.character(formals(ode)$method)
methods_fit <- as.character(formals(modFit)$method)

cat(paste("grind.R (",grind_version,") was sourced.\n\nNOTE! This version is no longer updated. For the latest version of grind, use install.packages().\n",sep=""))
