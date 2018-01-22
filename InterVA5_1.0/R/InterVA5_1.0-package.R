

#' Translation list of COD codes
#'
#' This is the translation of COD abbreviation codes into their corresponding
#' full names.
#'
#'
#' @name causetext
#' @docType data
#' @format A data frame with the translation of COD codes to their names on 68
#' CODs (both the version of COD only and COD with group code).
#' @keywords datasets
#' @examples
#'
#' data(causetext)
#'
NULL





#' Perform InterVA5 algorithm and provide graphical summarization of COD
#' distribution.
#'
#' Computes individual cause of death and population cause-specific mortality
#' fractions using the InterVA5 algorithm. Provides a simple graphical
#' representation of the result.
#'
#' To get the most up-to-date version of the package, as well as the past
#' versions, please check the github repository at:
#' \url{https://github.com/jarathomas/InterVA5-R-Replicate/}
#'
#' \tabular{ll}{ Package: \tab InterVA5\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2018-02-01\cr License: \tab GPL-3\cr }
#'
#' @name InterVA5-package
#' @docType package
#' @author Jason Thomas, Zehang Li, Tyler McCormick, Sam Clark
#'
#' Maintainer: Jason Thomas <jarathomas@@gmail.com>
#' @references http://www.interva.net/
#' @keywords InterVA
#' @examples
#'
#' data(SampleInput)
#' sample.output <- InterVA5(SampleInput, HIV = "h", Malaria = "v", directory = "VA test",
#'     filename = "VA_result", output = "extended", append = FALSE)
#'
NULL





#' Conditional probability of InterVA5
#'
#' This is the table of conditional probabilities of symptoms given CODs. The
#' values are from InterVA-5
#'
#'
#' @name probbase
#' @docType data
#' @format A data frame with 354 observations on 87 variables. Each observation
#' is the conditional probability.
#' @keywords datasets
#' @examples
#'
#' data(probbase)
#'
NULL

#' 10 records of Sample Input
#'
#' This is a dataset consisting of 10 arbitrary sample input deaths in the
#' acceptable format of InterVA5. Any data that needs to be analyzed by this
#' package should be in the same format. The orders of the input fields must
#' not be changed.
#'
#'
#' @name SampleInput
#' @docType data
#' @format 10 arbitrary input records.
#' @keywords datasets
#' @examples
#'
#' data(SampleInput)
#'
NULL



