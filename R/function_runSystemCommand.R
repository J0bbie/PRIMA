#' @title Perform a system command.
#' 
#' @description Performs a system command and captures STDERR and STDOUT in a list.
#' If the exitStatus is not zero (0), then the command failed.
#'
#' @param cmd (character): Command to execute.
#' 
#' @examples
#'  r <- runSystemCommand('ls -lh | head -n 5')
#' 
#' @return (list) List of exitcode, stdout and stderr.
#' 
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @export
runSystemCommand <- function (cmd){
    
    # Input validation --------------------------------------------------------
    
    checkmate::assertCharacter(cmd)
    
    # Start -------------------------------------------------------------------
    
    futile.logger::flog.trace('runSystemCommand - %s', cmd)
    
    # Temp. files to store STDERR/STDOUT.
    stderrFile <- base::tempfile(pattern ="R_robust.system_stderr", fileext = base::as.character(Sys.getpid()))
    stdoutFile <- base::tempfile(pattern ="R_robust.system_stdout", fileext = base::as.character(Sys.getpid()))
    
    # Perform command.
    cmd.output <- list()
    cmd.output$exitStatus <- base::system(paste0(cmd, " 2> ", base::shQuote(stderrFile), " > ", base::shQuote(stdoutFile)))
    cmd.output$stdout <- base::readLines(stdoutFile)
    cmd.output$stderr <- base::readLines(stderrFile)
    
    # Remove temp. files.
    base::unlink(c(stdoutFile, stderrFile))
    
    if(cmd.output$exitStatus != 0){
        warning(base::sprintf('%s returned error: %s', cmd, cmd.output$stderr))
    }
    
    return(cmd.output)
}