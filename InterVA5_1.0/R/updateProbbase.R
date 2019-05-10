#' Update the Symptom-Cause-Information source (aka probbaseV5).
#'
#' The function takes input of a list of va object and calculates the
#' cause-specific mortality fraction. It only calculates CSMF5 as aggregation of up to the third largest causes.
#'
#' @return \item{newProbbase}{The Symptom-Cause-Information ("Probbase") used to assign causes of death.}
#' @author Jason Thomas, Zehang LI, Tyler McCormick, Sam Clark
#' @keywords interVA
#' @seealso \code{\link{probbaseV5}}
#' @export update.SCI
#' @examples
#'
#' \dontrun{
#' data(RandomVA5)
#' RandomVA5 <- RandomVA5[1:2, ]
#' newProbbase <- update.SCI()
#' out <- InterVA5(RandomVA5, sci = newProbbase, HIV = "h", Malaria = "l", write=FALSE, 
#'     directory = tempdir(), filename = "VA5_result", output = "extended", append = FALSE)
#' }
#'
update.SCI <- function(){

    zipFile = tempfile(fileext = ".zip")
    curl::curl_download("http://www.byass.uk/interva/InterVA_5_v5.0_release.zip", zipFile)
    utils::unzip(zipfile = zipFile, files = "probbase.xls", exdir = tempdir())
    xlsName <- paste(tempdir(), "probbase.xls", sep = "/")
    newProbbase <- readxl::read_xls(xlsName)
    newProbbase <- newProbbase[-1, ]

    return(newProbbase)
}
