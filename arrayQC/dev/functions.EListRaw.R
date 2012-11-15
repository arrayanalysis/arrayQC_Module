assign("[.EListRaw",
	function (object, i, j, ...) 
	{
		if (nargs() != 3) 
			stop("Two subscripts required", call. = FALSE)
		other <- names(object$other)
		if (missing(i)) 
			if (missing(j)) 
				return(object)
			else {
				object$E <- object$E[, j, drop = FALSE]
				object$Eb <- object$Eb[, j, drop = FALSE]
				object$weights <- object$weights[, j, drop = FALSE]
				object$targets <- object$targets[j, , drop = FALSE]
				if (!is.null(object$design)) {
					object$design <- as.matrix(object$design)[j, 
					  , drop = FALSE]
					if (!is.fullrank(object$design)) 
					  warning("subsetted design matrix is singular", 
						call. = FALSE)
				}
				for (a in other) object$other[[a]] <- object$other[[a]][, 
					j, drop = FALSE]
			}
		else {
			if (is.character(i)) {
				i <- match(i, rownames(object))
				i <- i[!is.na(i)]
			}
			if (missing(j)) {
				object$E <- object$E[i, , drop = FALSE]
				object$Eb <- object$Eb[i, , drop = FALSE]
				object$weights <- object$weights[i, , drop = FALSE]
				object$genes <- object$genes[i, , drop = FALSE]
				for (a in other) object$other[[a]] <- object$other[[a]][i, 
					, drop = FALSE]
			}
			else {
				object$E <- object$E[i, j, drop = FALSE]
				object$Eb <- object$Eb[i, j, drop = FALSE]
				object$weights <- object$weights[i, j, drop = FALSE]
				object$genes <- object$genes[i, , drop = FALSE]
				object$targets <- object$targets[j, , drop = FALSE]
				if (!is.null(object$design)) {
					object$design <- as.matrix(object$design)[j, 
					  , drop = FALSE]
					if (!is.fullrank(object$design)) 
					  warning("subsetted design matrix is singular", 
						call. = FALSE)
				}
				for (a in other) object$other[[a]] <- object$other[[a]][i, 
					j, drop = FALSE]
			}
		}
		object
	}
)

rbind.EListRaw <- function (..., deparse.level = 1) 
{
    objects <- list(...)
    nobjects <- length(objects)
    out <- objects[[1]]
    other <- names(objects[[1]]$other)
    if (nobjects > 1) 
        for (i in 2:nobjects) {
            out$E <- rbind(out$E, objects[[i]]$E)
			out$Eb <- rbind(out$Eb, objects[[i]]$Eb)
            out$weights <- rbind(out$weights, objects[[i]]$weights)
            out$genes <- rbind(out$genes, objects[[i]]$genes)
            for (a in other) out$other[[a]] <- rbind(out$other[[a]], 
                objects[[i]]$other[[a]])
        }
    out
}

cbind.EListRaw <- function (..., deparse.level = 1) 
{
    objects <- list(...)
    nobjects <- length(objects)
    out <- objects[[1]]
    other <- names(objects[[1]]$other)
    if (nobjects > 1) 
        for (i in 2:nobjects) {
            out$E <- cbind(out$E, objects[[i]]$E)
			out$Eb <- cbind(out$Eb, objects[[i]]$Eb)
            out$weights <- cbind(out$weights, objects[[i]]$weights)
            out$targets <- rbind(out$targets, objects[[i]]$targets)
            for (a in other) out$other[[a]] <- cbind(out$other[[a]], 
                objects[[i]]$other[[a]])
        }
    out
}