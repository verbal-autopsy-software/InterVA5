#' Summarize population level cause-specific mortality fraction as InterVA5 suggested.
#'
#' The function takes input of a list of va object and calculates the
#' cause-specific mortality fraction. It only calculates CSMF5 as aggregation of up to the third largest causes.
#'
#' @param va The list of va object to summarize.
#' @return \item{dist.cod}{The cause-specific mortality fraction (including undetermined category).}
#' @author Jason Thomas, Zehang LI, Tyler McCormick, Sam Clark
#' @keywords interVA
#' @seealso \code{\link{CSMF5}}
#' @export CSMF.interVA5
#' @examples
#'
#' \dontrun{
#' data(RandomVA5)
#' # only fit first 20 observations for a quick illustration
#' RandomVA5 <- RandomVA5[1:20, ]
#' 
#' sample.output <- InterVA5(RandomVA5, HIV = "h", Malaria = "v", write=TRUE,
#'        directory = tempdir(), filename = "VA5_result", output = "extended", append = FALSE)
#' ## Get CSMF without plots
#' csmf <- CSMF.interVA5(sample.output$VA5)
#' }
#'
CSMF.interVA5 <- function(va) {
   # for future compatibility with non-standard input
    for (i in 1:length(va)) {
        if (!is.null(va[[i]]$wholeprob)) {
            causenames <- names(va[[i]]$wholeprob)
            causeindex <- 1:length(causenames)
            break
        }
    }

    include.probAC <- FALSE
    # fix for removing the first 3 preg related death in standard input
    if (causenames[ 1] == "Not pregnant or recently delivered" &&
        causenames[ 2] == "Pregnancy ended within 6 weeks of death" &&
        causenames[ 3] == "Pregnant at death"&&
        causenames[65] == "Culture" &&
        causenames[66] == "Emergency" &&
        causenames[67] == "Health systems" &&
        causenames[68] == "Inevitable" &&
        causenames[69] == "Knowledge" &&
        causenames[70] == "Resources") {
          causeindex <- causeindex[-c(1:3, 65:70)]
          causenames <- causenames[-c(1:3, 65:70)]
          include.probAC <- TRUE
    }

    ## Check if there is a valid va object
    if (length(va) < 1) {
        cat("No va object found")
        return()
    }
    ## Initialize the population distribution
    dist <- NULL
    for (i in 1:length(va)) {
        if (!is.null(va[[i]][15])) {
            dist <- rep(0, length(unlist(va[[i]][15])))
            break
        }
    }
    undeter <- 0

    ## pick not simply the top 3 causes, but the top 3 causes reported by InterVA5
    for (i in 1:length(va)) {
        if (is.null(va[[i]][15])) next
        this.dist <- unlist(va[[i]][15])
        if (include.probAC) this.dist[c(1:3, 65:70)] <- 0
        if (max(this.dist) < 0.4) {
            this.undeter <- ifelse(sum(this.dist) == 0, 1, sum(this.dist))
            undeter <- undeter + this.undeter
        } else {
            cutoff.3 <- this.dist[order(this.dist, decreasing = TRUE)[3]]
            cutoff.2 <- this.dist[order(this.dist, decreasing = TRUE)[2]]
            cutoff.1 <- this.dist[order(this.dist, decreasing = TRUE)[1]]
            cutoff <- min(max(cutoff.1 * 0.5 , cutoff.3), max(cutoff.1 * 0.5 , cutoff.2))

            undeter <- undeter + sum(this.dist[which(this.dist < cutoff)])
            this.dist[which(this.dist < cutoff)] <- 0
            if (!is.null(va[[i]][15])) dist <- dist + this.dist
        }
    }
    ## Normalize the probability for CODs
    if (undeter > 0) {
        dist.cod <- c(dist[causeindex], undeter)
        dist.cod <- dist.cod / sum(dist.cod)
        names(dist.cod) <- c(causenames, "Undetermined")
    } else {
        dist.cod <- dist[causeindex]/sum(dist[causeindex])
        names(dist.cod) <- causenames
    }
    if (sum(is.nan(dist.cod)) == length(dist.cod)) {
        dist.cod[is.nan(dist.cod)] <- 0
    }
    return(dist.cod)
}

#' Summarize population level mortality fraction by Circumstance of Mortality Category
#'
#' The function takes input of a list of va object and calculates the
#' mortality fraction by Circumstance of Mortality Category. 
#'
#' @param va The list of va object to summarize.
#' @return \item{dist.cod}{The cause-specific mortality fraction (including undetermined category).}
#' @author Jason Thomas, Zehang LI, Tyler McCormick, Sam Clark
#' @keywords interVA
#' @seealso \code{\link{CSMF5}}
#' @export COMCAT.interVA5
#' @examples
#'
#' \dontrun{
#' data(RandomVA5)
#' # only fit first 20 observations for a quick illustration
#' RandomVA5 <- RandomVA5[1:20, ]
#' 
#' sample.output <- InterVA5(RandomVA5, HIV = "h", Malaria = "v", 
#'        write=TRUE, directory = tempdir(),
#'        filename = "VA5_result", output = "extended", append = FALSE)
#' ## Get CSMF without plots
#' comcat <- COMCAT.interVA5(sample.output$VA5)
#' }
#'
COMCAT.interVA5 <- function(va){
   # for future compatibility with non-standard input
    for(i in 1:length(va)){
        if(!is.null(va[[i]]$wholeprob)){
            causenames <- names(va[[i]]$wholeprob)[65:70]
            causeindex <- 65:70
            break
        }
    }
    ## Check if there is a valid va object
    if(length(va) < 1){
        cat("No va object found")
        return()
    }
    ## Initialize the population distribution
    dist <- rep(0, 6)
    multi <- 0

    ## pick not simply the top 3 causes, but the top 3 causes reported by InterVA5
    for(i in 1:length(va)){
        if(is.null(va[[i]][15])) next
        this.dist <- unlist(va[[i]][15])[65:70]
        if(max(this.dist) < 0.5){
          multi <- multi + 1
        }else{
            dist[which.max(this.dist)] <- dist[which.max(this.dist)] + 1
        }
    }
    ## Normalize the probability for CODs
    if(multi > 0){
        dist.cod <- c(dist, multi)
        dist.cod <- dist.cod/sum(dist.cod)
        names(dist.cod)<-c(causenames, "Multiple")
    }else{
        dist.cod <- dist/sum(dist)
        names(dist.cod)<-causenames
    }
    if(sum(is.nan(dist.cod)) == length(dist.cod)){
        dist.cod[is.nan(dist.cod)] <- 0
    }

    return(dist.cod)
}

#' Summarize and plot a population level distribution of va probabilities.
#'
#' The function takes input of a list of va object and produces a summary plot
#' for the population distribution.
#'
#'
#' @param va The list of va object to summarize.
#' @param top.aggregate Integer indicating how many causes from the top need to go into
#' summary. The rest of the probabilities goes into an extra category
#' "Undetermined".  When set to NULL, default is all causes to be considered.
#' This is only used when \code{InterVA.rule} set to "FALSE".
#' @param InterVA.rule If it is set to "TRUE", only the top 3 causes reported by
#' InterVA5 is calculated into CSMF as in InterVA5. The rest of probabilities
#' goes into an extra category "Undetermined". Default set to "FALSE".
#' @param noplot A logical value indicating whether the plot will be shown. If
#' it is set to "TRUE", only the CSMF will be returned.
#' @param title A character string for the title of the CSMF plot.  
#' @param top.plot the maximum number of causes to plot in bar plot
#' @param min.prob The minimum probability that is to be plotted in bar chart,
#' or to be labeled in pie chart.
#' @param type An indicator of the type of chart to plot.  "pie" for pie chart;
#' "bar" for bar chart.
#' @param return.barplot A logical indicating if the (barplot) ggplot() object
#' should be returned (instead of printed).  Default value is FALSE.
#' @param ... Arguments to be passed to/from graphic function
#' \code{\link[graphics]{barplot}}, \code{\link[graphics]{pie}}, and more
#' graphical paramters (see \code{\link[graphics]{par}}). They will affect the
#' main title, size and font of labels, and the radius of the pie chart.
#' @return \item{dist.cod}{The population probability of CODs.}
#' @author Jason Thomas, Zehang LI, Tyler McCormick, Sam Clark
#' @seealso \code{\link{CSMF.interVA5}}
#' @keywords interVA
#' @export CSMF5
#' @examples
#'
#' \dontrun{
#' data(RandomVA5)
#' # only fit first 20 observations for a quick illustration
#' RandomVA5 <- RandomVA5[1:20, ]
#' 
#' sample.output <- InterVA5(RandomVA5, HIV = "h", Malaria = "v", write = FALSE, 
#'        directory = tempdir(), filename = "VA5_result", output = "extended", 
#'        append = FALSE)
#'
#' ## Get CSMF by considering only top 3 causes reported by InterVA5.
#' ## This is equivalent to using CSMF.interVA5() command Note that
#' ## it's different from using all top 3 causses, since they may not
#' ## all be reported
#' CSMF.summary <- CSMF5(sample.output, InterVA.rule = TRUE,
#'    noplot = TRUE)
#'
#' ## Population level summary using pie chart
#' CSMF.summary2 <- CSMF5(sample.output, type = "pie",
#'  min.prob = 0.01, title = "population COD distribution using pie chart",
#'  clockwise = FALSE, radius = 0.7, cex = 0.7, cex.main = 0.8)
#'
#' ## Population level summary using bar chart
#' CSMF.summary3 <- CSMF5(sample.output, type = "bar",
#'   min.prob = 0.01, title = "population COD distribution using bar chart",
#'   cex.main = 1)
## Population level summary specified by number of top causes
#' CSMF.summary4 <- CSMF5(sample.output, type = "bar",
#'   top.plot = 5, title = "Top 5 population COD distribution",
#'   cex.main = 1)
#' }
#'
CSMF5 <- function (va, top.aggregate = NULL, InterVA.rule = FALSE, noplot = FALSE,
                   title = "Top CSMF Distribution", type = "bar",  top.plot = 10,
                   return.barplot = FALSE, min.prob = 0, ... ) {

    ## Check if there is a valid va object
    if (class(va) == "interVA5") {
        va <- va$VA5
    }

    # for future compatibility with non-standard input
    for (i in 1:length(va)) {
        if (!is.null(va[[i]]$wholeprob)) {
            causenames <- names(va[[i]]$wholeprob)
            causeindex <- 1:length(causenames)
            break
        }
    }

    include.probAC <- FALSE
    # fix for removing the first 3 preg related death in standard input
    if (causenames[1] == "Not pregnant or recently delivered" &&
        causenames[2] == "Pregnancy ended within 6 weeks of death" &&
        causenames[3] == "Pregnant at death"&&
        causenames[65] == "Culture" &&
        causenames[66] == "Emergency" &&
        causenames[67] == "Health systems" &&
        causenames[68] == "Inevitable" &&
        causenames[69] == "Knowledge" &&
        causenames[70] == "Resources") {
            causeindex <- causeindex[-c(1:3, 65:70)]
            causenames <- causenames[-c(1:3, 65:70)]
            include.probAC <- TRUE
    }

    if (length(va) < 1) {
        cat("No va object found")
        return()
    }
    ## Initialize the population distribution
    dist <- NULL
    for (i in 1:length(va)) {
        if (!is.null(va[[i]][15])) {
            dist <- rep(0, length(unlist(va[[i]][15])))
            break
        }
    }
    ## determine how many causes from top need to be summarized
    if (is.null(top.aggregate)) top.aggregate <- length(causeindex)
    undeter <- 0

    if (is.null(dist)) {cat("No va probability found in input"); return()}
    ## Add the probabilities together
    if (!InterVA.rule) {
        for (i in 1:length(va)) {
            if (is.null(va[[i]][15])) {
                undeter <- undeter + 1
                next
            }
            this.dist <- unlist(va[[i]][15])
            if (include.probAC) this.dist[c(1:3, 65:70)] <- 0
            if (sum(this.dist) == 0) {
                undeter <- undeter + 1
                next
            }
            cutoff <- this.dist[order(this.dist, decreasing = TRUE)[top.aggregate]]
            undeter <- undeter + sum(this.dist[which(this.dist < cutoff)])
            this.dist[which(this.dist < cutoff)] <- 0
            if (!is.null(va[[i]][15])) dist <- dist + this.dist
        }
            ## Normalize the probability for CODs
        if (undeter > 0) {
            dist.cod <- c(dist[causeindex], undeter)
            dist.cod <- dist.cod / sum(dist.cod)
            names(dist.cod) <- c(causenames, "Undetermined")
        } else {
            dist.cod <- dist[causeindex] / sum(dist[causeindex])
            names(dist.cod) <- causenames
        }
    } else {
        dist.cod <- CSMF.interVA5(va)
    }

    ## Check if there is CODs above the minimum cut-off for prob
    if (max(dist.cod) < min.prob) {
        cat("No COD larger than the minimum probability cut off line")
        return()
    }
    if (noplot) {
        return(dist.cod)
    }

    if (!is.null(top.plot)) {
        if (top.plot < length(dist.cod)) {
            thre <- sort(dist.cod, decreasing=TRUE)[top.plot]
            min.prob <- max(min.prob, thre)
        }
    }

    ## Make pie plot upon request
    if (type == "pie") {
        dist.cod.sort <- sort(dist.cod, decreasing=TRUE)
        pie.color <- grey.colors(length(dist.cod.sort[dist.cod.sort >= min.prob]))
        pie.color.left <- rep(pie.color[length(pie.color)], length(dist.cod.sort[dist.cod.sort < min.prob]))
        pie.color <- c(pie.color, pie.color.left)
        pie(dist.cod.sort, main = title,
            col = pie.color, labels = names(dist.cod.sort)[dist.cod.sort >= min.prob],
            ...)
    }
    ## Make bar plot upon request
    if (type == "bar") {
        ## dist.cod.min <- dist.cod[dist.cod >= min.prob ]
        dist.cod.min <- dist.cod[dist.cod >= min.prob ]
        dist.cod.min <- sort(dist.cod.min, decreasing = FALSE)
        if (requireNamespace("ggplot2", quietly = TRUE)) {
            barplot.df <- data.frame(Probability = dist.cod.min,
                                     Causes = names(dist.cod.min))
            g <- ggplot2::ggplot(barplot.df,
                                 ggplot2::aes(x = stats::reorder(Causes,
                                                                 seq(1:length(Causes))),
                                              y = Probability,
                                              fill = stats::reorder(Causes,
                                                                    seq(1:length(Causes))))) +
                ggplot2::geom_bar(stat="identity") +
                ggplot2::xlab("") +
                ggplot2::ylab("") +
                ggplot2::coord_flip() +
                ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
                ggplot2::ggtitle(title) +
                ggplot2::theme(legend.position = "none")
            if (return.barplot) {
                return(g)
            } else {
                par(las = 2)
                par(mar = c(5,15,4,2))
                print(g)
            }

        } else {
            bar.color <- grey.colors(length(dist.cod.min))
            bar.color <- rev(bar.color)
            barplot(dist.cod.min, horiz = TRUE, names.arg = names(dist.cod.min),
                    main = title, col = bar.color, cex.names=0.8,
                    xlab = "Probability", ...)
        }
    }
    ## Save the population distribution
    dist.cod
}

#' Plot an individual-level distribution of va probabilities.
#'
#' The function takes an input of a single va object and produces a summary plot
#' for it.
#'
#'
#' @param va A va object
#' @param min.prob The minimum probability that is to be plotted in bar chart,
#' or to be labeled in pie chart.
#' @param type An indicator of the type of chart to plot.  "pie" for pie chart;
#' "bar" for bar chart.
#' @param title A character string for the title of the CSMF plot.  
#' @param ... Arguments to be passed to/from graphic function
#' \code{\link[graphics]{barplot}}, \code{\link[graphics]{pie}}, and more
#' graphical paramters (see \code{\link[graphics]{par}}). They will affect the
#' main title, size and font of labels, and the radius of the pie chart.
#' @seealso \code{\link{CSMF5}}
#' @export InterVA5.plot
#' @keywords InterVA
#' @examples
#'
#' \dontrun{
#' data(RandomVA5)
#' # only fit first 20 observations for a quick illustration
#' RandomVA5 <- RandomVA5[1:20, ]
#' #' sample.output <- InterVA5(RandomVA5, HIV = "h", Malaria = "v", write = FALSE, 
#'     directory = tempdir(), filename = "VA5_result", output = "extended", append = FALSE)
#'
#' ## Individual level summary using pie chart
#' InterVA5.plot(sample.output$VA5[[3]], type = "pie", min.prob = 0.01,
#'     main = "1st sample VA analysis using pie chart", clockwise = FALSE,
#'     radius = 0.6, cex = 0.6, cex.main = 0.8)
#'
#'
#' ## Individual level summary using bar chart
#' InterVA5.plot(sample.output$VA5[[3]], type = "bar", min.prob = 0.01,
#'     main = "2nd sample VA analysis using bar chart", cex.main = 0.8)
#' }
#'
InterVA5.plot <- function(va, type = "bar", title = "Top CSMF Distribution", min.prob = 0.01, ... ){

    # for future compatibility with non-standard input
    if (!is.null(va$wholeprob)) {
        causenames <- names(va$wholeprob)
        causeindex <- 1:length(causenames)
    } else {
        cat("Cause of death undetermined for this case\n")
        return()
    }

    # fix for removing the first 3 preg related death in standard input
    if(causenames[1] == "Not pregnant or recently delivered" &&
        causenames[2] == "Pregnancy ended within 6 weeks of death" &&
        causenames[3] == "Pregnant at death"){
            causeindex <- causeindex[-c(1:3)]
            causenames <- causenames[-c(1:3)]
    }


    ## Check if there is a valid va object
    if (length(va) < 1) {
        cat("No va object found")
        return()
    }
    ## Find the probability distribution
    dist <- unlist(va[15])
    dist.cod <- dist[causeindex]/sum(dist[causeindex])
    ## Check if there is CODs above the minimum cut-off for prob
    if (max(dist.cod) < min.prob) {
        cat("No COD larger than the minimum probability cut off line")
        return()
    }
    names(dist.cod) <- causenames
    ## Make pie plot upon request
    if (type == "pie") {
        dist.cod.sort <- sort(dist.cod, decreasing=TRUE)
        pie.color <- grey.colors(length(dist.cod.sort[dist.cod.sort >= min.prob]))
        pie.color.left <- rep(pie.color[length(pie.color)], length(dist.cod.sort[dist.cod.sort < min.prob]))
        pie.color <- c(pie.color, pie.color.left)
        pie(dist.cod.sort, main = title, col = pie.color,labels = names(dist.cod.sort)[dist.cod.sort >= min.prob], ...)
    }
    ## Make bar plot upon request
    if (type == "bar") {
        dist.cod.min <- dist.cod[dist.cod >= min.prob]
        dist.cod.min <- sort(dist.cod.min, decreasing = FALSE)
        par(las = 2)
        par(mar = c(5,15,4,2))
        if (requireNamespace("ggplot2", quietly = TRUE)) {
            barplot.df <- data.frame(Probability = dist.cod.min,
                                     Causes = names(dist.cod.min))
            g <- ggplot2::ggplot(barplot.df,
                                 ggplot2::aes(x = stats::reorder(Causes,
                                                                 seq(1:length(Causes))),
                                              y = Probability,
                                              fill = stats::reorder(Causes,
                                                                    seq(1:length(Causes))))) +
                ggplot2::geom_bar(stat="identity") +
                ggplot2::xlab("") + ggplot2::ylab("") +
                ggplot2::coord_flip() +
                ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
                ggplot2::ggtitle(title) +
                ggplot2::theme(legend.position = "none")
            print(g)
        } else {
            bar.color <- grey.colors(length(dist.cod.min))
            bar.color <- rev(bar.color)
            barplot(dist.cod.min, horiz = TRUE, names.arg = names(dist.cod.min),
                    main = title, col = bar.color, cex.names=0.8,
                    xlab = "Probability", ...)
        }
    }
}
