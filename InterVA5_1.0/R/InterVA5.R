InterVA5 <- function (Input, HIV, Malaria, directory = NULL, filename = "VA_result", 
                      output = "classic", append = FALSE, groupcode = FALSE,
                      write = TRUE, 
                      ...) 
    ## FOR DEBUGGING: 
    ## load('data/SampleInput.rda'); load('data/causetext.rda'); load('data/probbase.rda'); Input <- SampleInput
    ## HIV='l'; Malaria='l'; directory = NULL; filename = "VA_result";
    ## output = "classic"; append = FALSE; groupcode = FALSE; write=TRUE
{
    va <- function(ID, MALPREV, HIVPREV, PREGSTAT, PREGLIK, PRMAT, 
                   INDET, CAUSE1, LIK1, CAUSE2, LIK2, CAUSE3, LIK3, wholeprob, 
                   ...) {
        ID <- ID
        MALPREV <- as.character(MALPREV)
        HIVPREV <- as.character(HIVPREV)
        PREGSTAT <- paste(PREGSTAT, paste(rep(" ", 5 - nchar(PREGSTAT)), 
                                          collapse = ""), collapse = "")
        PREGLIK <- PREGLIK
        PRMAT <- PRMAT
        INDET <- as.character(INDET)
        wholeprob <- wholeprob
        va.out <- list(ID = ID, MALPREV = MALPREV, HIVPREV = HIVPREV, 
                       PREGSTAT = PREGSTAT, PREGLIK = PREGLIK, PRMAT = PRMAT, 
                       INDET = INDET, CAUSE1 = CAUSE1, LIK1 = LIK1, CAUSE2 = CAUSE2, 
                       LIK2 = LIK2, CAUSE3 = CAUSE3, LIK3 = LIK3, wholeprob = wholeprob)
        va.out
    }
    save.va <- function(x, filename, write) { ## HERE -- need to see if taking out 14th column still applies (should it be a different col?)
        if (!write) {
            return()
        }
        x <- x[-14]
        x <- as.matrix(x)
        filename <- paste(filename, ".csv", sep = "")
        write.table(t(x), file = filename, sep = ",", append = TRUE, 
                    row.names = FALSE, col.names = FALSE)
    }
    save.va.prob <- function(x, filename, write) {
        if (!write) {
            return()
        }
        prob <- unlist(x[14])  ## HERE -- need to see if this still applies (and next line)
        x <- x[-14]
        x <- unlist(c(as.matrix(x), as.matrix(prob)))
        filename <- paste(filename, ".csv", sep = "")
        write.table(t(x), file = filename, sep = ",", append = TRUE, 
                    row.names = FALSE, col.names = FALSE)
    }
    if (is.null(directory)) 
        directory = getwd()
    dir.create(directory, showWarnings = FALSE)
    globle.dir <- getwd()
    setwd(directory)
    data("probbase", envir = environment())
    probbase <- get("probbase", envir = environment())
    probbase <- as.matrix(probbase)
    data("causetext", envir = environment())
    causetext <- get("causetext", envir = environment())
    ## HERE -- note sure if we want groupcode
    ## groupcode: A logical value indicating whether or not the group code
    ##            will be included in the output causes.
    ## HERE -- you may need to copy the "extra" column over from the previous data
    ## if (groupcode) {
    ##     causetext <- causetext[, -2]
    ## }
    ## else {
    ##     causetext <- causetext[, -3]
    ## }
    if (write) {
        cat(paste("Error log built for InterVA", Sys.time(), 
            "\n"), file = "errorlog.txt", append = FALSE)
        cat(paste("Warning log built for InterVA", Sys.time(), 
            "\n"), file = "warnings.txt", append = FALSE)
    }

    Input <- as.matrix(Input)
    if (dim(Input)[1] < 1) {
        stop("error: no data input")
    }
    N <- dim(Input)[1]
    S <- dim(Input)[2]
    if (S != dim(probbase)[1]) {
        stop("error: invalid data input format. Number of values incorrect")
    }
    if (tolower(colnames(Input)[S]) != "w610459o") {
        stop("error: the last variable should be 'w610459o -- q costs'")
    }
    data("SampleInput", envir = environment())
    SampleInput <- get("SampleInput", envir = environment())
    valabels = colnames(SampleInput)
    count.changelabel = 0
    for (i in 1:S) {
        if (tolower(colnames(Input)[i]) != tolower(valabels)[i]) {
            warning(paste("Input column '", colnames(Input)[i], "' does not match InterVA standard: '", valabels[i], "'", sep = ""),
                    call. = FALSE, immediate. = TRUE)
            count.changelabel = count.changelabel + 1
        }
    }
    if (count.changelabel > 0) {
        warning(paste(count.changelabel, "column names changed in input. \n If the change in undesirable, please change in the input to match standard InterVA5 input format."), 
            call. = FALSE, immediate. = TRUE)
        colnames(Input) <- valabels
    }
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "I"  ] <- 1
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "A+" ] <- 0.8
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "A"  ] <- 0.5
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "A-" ] <- 0.2
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "B+" ] <- 0.1
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "B"  ] <- 0.05
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "B-" ] <- 0.02
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "B -"] <- 0.02
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "C+" ] <- 0.01
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "C"  ] <- 0.005
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "C-" ] <- 0.002
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "D+" ] <- 0.001
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "D"  ] <- 5e-04
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "D-" ] <- 1e-04
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "E"  ] <- 1e-05
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == "N"  ] <- 0
    probbase[,18:ncol(probbase)][probbase[,18:ncol(probbase)] == ""   ] <- 0
    probbase[1, 1:17] <- rep(0, 17)
    Sys_Prior <- as.numeric(probbase[1, ])
    D <- length(Sys_Prior)
    HIV <- tolower(HIV)
    Malaria <- tolower(Malaria)
    if (!(HIV %in% c("h", "l", "v")) || !(Malaria %in% c("h","l", "v"))) {
        stop("error: the HIV and Malaria indicator should be one of the three: 'h', 'l', and 'v'")
    }
    if (HIV == "h") 
        Sys_Prior[23] <- 0.05
    if (HIV == "l") 
        Sys_Prior[23] <- 0.005
    if (HIV == "v") 
        Sys_Prior[23] <- 1e-05
    if (Malaria == "h") {
        Sys_Prior[25] <- 0.05 
        Sys_Prior[45] <- 0.05 
    }
    if (Malaria == "l") {
        Sys_Prior[25] <- 0.005
        Sys_Prior[45] <- 1e-05
    }
    if (Malaria == "v") {
        Sys_Prior[25] <- 1e-05
        Sys_Prior[45] <- 1e-05
    }
    ID.list <- rep(NA, N)
    VAresult <- vector("list", N)
    if (write && append == FALSE) {
        header = c("ID", "MALPREV", "HIVPREV", "PREGSTAT", "PREGLIK", 
            "PRMAT", "INDET", "CAUSE1", "LIK1", "CAUSE2", "LIK2", 
            "CAUSE3", "LIK3", "INDET", "COMCAT", "COMNUM")
        if (output == "extended") 
            header = c(header, as.character(causetext[, 2])) ## HERE -- in previous version causetext had 3 columns
        write.table(t(header), file = paste(filename, ".csv", sep = ""),
                    row.names = FALSE, col.names = FALSE, sep = ",")
    }
    nd <- max(1, round(N/100))
    np <- max(1, round(N/10))
    
    for (i in 1:N) {
        if (i%%nd == 0) {
            cat(".")
        }
        if (i%%np == 0) {
            cat(paste(round(i/N * 100), "% completed\n", sep = "")) ## HERE -- need to fix this
        }
        index.current <- as.character(Input[i, 1])   
        Input[i, which(is.na(Input[i, ]))] <- "0"    
        Input[i, which(toupper(Input[i, ]) != "Y")] <- "0"
        Input[i, which(toupper(Input[i, ]) == "Y")] <- "1"
        input.current <- as.numeric(Input[i, ])
        input.current[1] <- 0                      
        if (sum(input.current[6:12]) < 1) {
            if (write) {
                cat(paste(index.current, " Error in age indicator: Not Specified ", 
                          "\n"), file = "errorlog.txt", append = TRUE)
            }
            next
        }
        if (sum(input.current[4:5]) < 1) {
            if (write) {
                cat(paste(index.current, " Error in sex indicator: Not Specified ", 
                          "\n"), file = "errorlog.txt", append = TRUE)
            }
            next
        }
        if (sum(input.current[20:354]) < 1) { ## HERE - not sure where to start/stop (not sure how InterVA4 determines "indicators" -- stops at 223, but there are 246 columns)
            if (write) {
                cat(paste(index.current, " Error in indicators: No symptoms specified ", 
                          "\n"), file = "errorlog.txt", append = TRUE)
            }
            next
        }
        
        ## HERE -- Byass also checks the "samb" column in data processing -- part of deleting incomplete data
        ## HERE -- need this error message "the following records are incomplete and excluded from further processing:" which excludes row -- nodata
        
        for (k in 1:2) {
            for (j in 2:S) {
                
                if(i==1 & k==1 & j==2 & write ){
                    cat("\n\n the following data discrepancies were identified and handled: \n\n", file="warnings.txt", append=TRUE)
                }
                pass.iteration <- ifelse(k==1, "(first pass)", "(second pass)")

                if (input.current[j] == 1) {     
                    
                    Dont.ask.tmp   <- substr(probbase[j, 8:15],1,5)
                    Dont.ask.indic <- Dont.ask.tmp[Dont.ask.tmp!=""]
                    Dont.ask.who   <- probbase[match(toupper(Dont.ask.indic), toupper(probbase[,1])),4]
                    Dont.ask.rows  <- match(toupper(Dont.ask.indic), toupper(probbase[,1]))
                    Dont.ask.list  <- input.current[match(toupper(Dont.ask.who), toupper(colnames(Input)))]
                    
                    Dont.ask.val.tmp <- substr(probbase[j, 8:15],6,6)
                    Dont.ask.val     <- Dont.ask.val.tmp[Dont.ask.val.tmp!=""]
                    ##Dont.ask.val <- probbase[match(toupper(Dont.ask.indic), toupper(probbase[,1])),6]
                    Dont.ask.val[Dont.ask.val=="N"] <- 0
                    Dont.ask.val[Dont.ask.val=="Y"] <- 1
                    
                    Dont.ask.error <- Dont.ask.list==as.numeric(Dont.ask.val)
                    
                    if (sum(Dont.ask.error)>0) {
                        input.current[j] <- 0
                        if (write) {
                            cat(index.current, "   ",
                                paste(probbase[j, 4], " (", probbase[j, 3], ") value inconsistent with ",
                                      Dont.ask.who[which(Dont.ask.error==TRUE)],
                                      " (", probbase[Dont.ask.rows[which(Dont.ask.error==TRUE)], 3], ")",  sep=""),
                                " - cleared in working file", pass.iteration, "\n",
                                file = "warnings.txt", append = TRUE)
                        }
                    }
                }
                if (input.current[j] == 1) {
                    Ask.if.tmp <- substr(probbase[j, 16],1,5)
                    Ask.if     <- probbase[match(toupper(Ask.if.tmp), toupper(probbase[,1])),4]
                    Ask.if.row <- match(toupper(Ask.if.tmp), toupper(probbase[,1]))
                    
                    ##Ask.if.val <- probbase[match(toupper(Ask.if.tmp), toupper(probbase[,1])),6]
                    Ask.if.val <- substr(probbase[j, 16],6,6)
                    Ask.if.val[Ask.if.val=="N"] <- 0
                    Ask.if.val[Ask.if.val=="Y"] <- 1
                    
                    if (!is.na(match(toupper(Ask.if), toupper(colnames(Input))))) {
                        if (input.current[match(toupper(Ask.if), 
                                                toupper(colnames(Input)))] != as.numeric(Ask.if.val)) {
                            input.current[match(toupper(Ask.if), toupper(colnames(Input)))] <- as.numeric(Ask.if.val)
                            if (write) {
                                cat(index.current, "   ",
                                    paste(probbase[j, 4], " (", probbase[j, 3], ")  not flagged in category ", 
                                          Ask.if, " (", probbase[Ask.if.row, 3], sep=""),
                                    ") - updated in working file", pass.iteration, "\n",
                                    file = "warnings.txt", append = TRUE)
                            }
                        }
                    }
                }
                if (input.current[j] == 1) {
                    NN.only.tmp <- substr(probbase[j, 17],1,5)
                    NN.only     <- probbase[match(toupper(NN.only.tmp), toupper(probbase[,1])),4]

                    NN.only.val <- substr(probbase[j, 17],6,6)
                    NN.only.val[NN.only.val=="N"] <- 0
                    NN.only.val[NN.only.val=="Y"] <- 1
                    
                    if (!is.na(match(toupper(NN.only), toupper(colnames(Input))))) { 
                        if (input.current[match(toupper(NN.only), 
                                                toupper(colnames(Input)))] != NN.only.val) {
                            input.current[j] <- 0
                            if (write) {
                                cat(index.current, "   ",
                                    paste(probbase[j, 4], " (", probbase[j, 3], sep=""),
                                    ")  only required for neonates - cleared in working file", pass.iteration, "\n",
                                    file = "warnings.txt", append = TRUE)
                            }
                        }
                    }
                }
            }
        }
        
        reproductiveAge <- 0
        preg_state      <- " "
        lik.preg        <- " "
        if (input.current[5] == 1 && (input.current[8] == 1 || input.current[9] == 1)) {
            reproductiveAge <- 1
        }
        prob <- Sys_Prior[18:D]
        temp <- which(input.current[2:length(input.current)] == 1)
        for (jj in 1:length(temp)) {
            temp_sub <- temp[jj]
            for (j in 18:D) {
                prob[j - 17] <- prob[j - 17] * as.numeric(probbase[temp_sub + 1, j])
            }
            if (sum(prob[1:3]) > 0)
                prob[1:3] <- prob[1:3]/sum(prob[1:3])
            if (sum(prob[4:64]) > 0)
                prob[4:64] <- prob[4:64]/sum(prob[4:64])
        }
        names(prob) <- causetext[, 2]
        prob_A <- prob[ 1: 3]
        prob_B <- prob[ 4:64]
        prob_C <- prob[65:70]
        if (sum(prob_A) == 0 || reproductiveAge == 0) {
            preg_state <- "Indet"
            lik.preg <- 0
        }
        if (which.max(prob_A) == 1 && prob_A[1] != 0 && reproductiveAge == 1) {
            preg_state <- "nrp"
            lik.preg <- round(prob_A[1]/sum(prob_A) * 100)
        }
        if (which.max(prob_A) == 2 && prob_A[2] != 0 && reproductiveAge == 1) {
            preg_state <- "pr6w"
            lik.preg <- round(prob_A[2]/sum(prob_A) * 100)
        }
        if (which.max(prob_A) == 3 && prob_A[3] != 0 && reproductiveAge == 1) {
            preg_state <- "preg"
            lik.preg <- round(prob_A[3]/sum(prob_A) * 100)
        }
        lik_mat <- " "
        if (reproductiveAge == 1 && sum(prob_A) != 0) 
            lik_mat <- round((prob_A[2] + prob_A[3])/sum(prob_A) * 100)
        if (sum(prob_B) != 0) 
            prob_B <- prob_B/sum(prob_B)
        prob.temp <- prob_B
        if (max(prob.temp) <= 0.4) {
            indet <- "Indet"
            cause1 <- lik1 <- cause2 <- lik2 <- cause3 <- lik3 <- " "
        }
        if (max(prob.temp) > 0.4) {
            indet <- " "
            lik1 <- round(max(prob.temp) * 100)
            cause1 <- names(prob.temp)[which.max(prob.temp)]
            prob.temp <- prob.temp[-which.max(prob.temp)]
            lik2 <- round(max(prob.temp) * 100)
            cause2 <- names(prob.temp)[which.max(prob.temp)]
            if (max(prob.temp) < 0.5 * max(prob_B)) 
                lik2 <- cause2 <- " "
            prob.temp <- prob.temp[-which.max(prob.temp)]
            lik3 <- round(max(prob.temp) * 100)
            cause3 <- names(prob.temp)[which.max(prob.temp)]
            if (max(prob.temp) < 0.5 * max(prob_B)) 
                lik3 <- cause3 <- " "
        }
        ID.list[i] <- index.current
        VAresult[[i]] <- va(ID = index.current, MALPREV = Malaria, 
                            HIVPREV = HIV, PREGSTAT = preg_state, PREGLIK = lik.preg, 
                            PRMAT = lik_mat, INDET = indet, CAUSE1 = cause1, 
                            LIK1 = lik1, CAUSE2 = cause2, LIK2 = lik2, CAUSE3 = cause3, 
                            LIK3 = lik3, wholeprob = c(prob_A, prob_B))
        if (output == "classic") 
            save.va(VAresult[[i]], filename = filename, write)
        if (output == "extended") 
            save.va.prob(VAresult[[i]], filename = filename, write)
    }
    setwd(globle.dir)
    out <- list(ID = ID.list[which(!is.na(ID.list))], VA = VAresult[which(!is.na(ID.list))], 
                Malaria = Malaria, HIV = HIV)
    class(out) <- "interVA"
    return(out)
}
#<environment: namespace:InterVA5>
 
