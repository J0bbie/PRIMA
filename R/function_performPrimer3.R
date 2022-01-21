#' @title Generate PCR primers (local primer3).
#' 
#' @description Performs a local primer3 analysis to generate PCR primers based on a given DNA or RNA sequence. 
#' Sequences which \strong{must} be included in the final PCR product can be given via \emph{target}, this will generate flanking primers accordingly.
#' 
#' Please see the primer3 instructions on more detailed information on the primer3-specific settings.
#' In addition, the path to the primer3_config folder need to specified and (if not detected) also the primer3 executable. 
#' 
#' @param input.sequence (\link[Biostrings]{XStringSet}): A DNAString or RNAString object containing the sequence.
#' @param target (\link[Biostrings]{XStringSet} or numeric): A DNAString or RNAString object containing a sequence that must be present in the PCR product.
#' Alternatively, a numeric vector can be given to include that starting position + width, e.g. c(10, 15) to include position 10 to 25.
#' @param exclude (\link[Biostrings]{XStringSet} or numeric): A DNAString or RNAString object containing a sequence that must not be present in the PCR product.
#' Alternatively, a numeric vector can be given to exclude that starting position + width, e.g. c(10, 15) to exlude position 10 to 25.
#' @param product.sizes (character): A vector of size ranges. Primer3 first tries to pick primers in the first range. If that is not possible, it goes to the next range and tries again.
#' @param primer.size (character): Min. optimal and max. primer sizes. E.g. c(18, 20, 23).
#' @param primer.tm (character): Min. optimal and max. primer Tm. E.g. c(57, 59, 62).
#' @param nPrimer (numeric): Max. number of primers pairs to return.
#' @param internalOligo (logical): Also generate an internal oligo?
#' @param mismatching.library (character): Path to a primer mismatching library (provided by primer3). Leave NULL to ignore.
#' @param primer3.args (list): Additional arguments for primer3 settings. E.g. list('PRIMER_PICK_INTERNAL_OLIGO' = 0, 'PRIMER_EXPLAIN_FLAG' = 1).
#' @param primer3_core.binary (character): Optional location of the binary of primer3_core (if not found by Sys.which())
#' @param primer3_config (character): Path to primer3_config folder.
#' @param primer3_defaultsettings (character): Path to primer3 default settings file.
#' 
#' @examples
#'  \dontrun{
#'  
#'  DNA.MDM2 <- Biostrings::DNAString('ATGTGCAATACCAACATGTCTGTACCTACTGATGGTGCTGTAACCACCTCACAGATT
#'  CCAGCTTCGGAACAAGAGACCCTGGTTCTTTTTTATCTTGGCCAGTATATTATGACTAAACGATTATATGATGAGAAGCAACAACATATTGTA')
#'  primer3.results <- performPrimer3(DNA.MDM2, target = Biostrings::DNAString('TGGTTCTTTTT'), nPrimer = 3)
#'  
#' 	}
#' @return (\link[data.table]{data.table}) Results of primer3 in a data.table.
#' 
#' @author Job van Riet \email{j.vanriet@erasmusmc.nl}
#' @export
performPrimer3 <- function(input.sequence, target = NULL, exclude = NULL, product.sizes = c('50-150', '150-250'), primer.size = c(18, 20, 23), primer.tm = c(57, 59, 62), nPrimer = 3, 
                           internalOligo = FALSE, mismatching.library = NULL, primer3.args = NULL, primer3_core.binary = NULL, primer3_config = NULL, primer3_defaultsettings = NULL){
  
  # Input validation. ----
  
  checkmate::assert(
    checkmate::checkClass(input.sequence, 'DNAString'),
    checkmate::checkClass(input.sequence, 'RNAString')
  )
  
  checkmate::assert(
    checkmate::checkClass(target, 'DNAString', null.ok = T),
    checkmate::checkClass(target, 'RNAString', null.ok = T),
    checkmate::checkNumeric(target, len = 2, null.ok = T)
  )
  
  checkmate::assert(
    checkmate::checkClass(exclude, 'DNAString', null.ok = T),
    checkmate::checkClass(exclude, 'RNAString', null.ok = T),
    checkmate::assertNumeric(exclude, len = 2, null.ok = T)
  )
  
  checkmate::assertCharacter(product.sizes)
  checkmate::assertNumeric(primer.size)
  checkmate::assertNumeric(primer.tm)
  checkmate::assertNumber(nPrimer)
  checkmate::assertLogical(internalOligo)
  checkmate::assertCharacter(mismatching.library, null.ok = T)
  checkmate::assertList(primer3.args, null.ok = T)
  checkmate::assertCharacter(primer3_core.binary, null.ok = T)
  checkmate::assertDirectory(primer3_config, access = 'r')
  
  # Check primer3 settings. If none given is given, use the default settings.
  if(base::is.null(primer3_defaultsettings)){
    primer3_defaultsettings <- '~/git/PRIMA/inst/extdata/primer3_v1_1_4_default_settings.txt'
  }
  checkmate::assertAccess(primer3_defaultsettings, 'r')
  
  # Check if primer3 can be executed.
  if(base::is.null(primer3_core.binary)) primer3_core.binary <- Sys.which('primer3_core')
  checkmate::assertAccess(base::Sys.which(primer3_core.binary), 'x')
  
  
  # Start ----
  
  futile.logger::flog.debug('PRIMA - Performing primer3: %s', primer3_core.binary)
  
  # Create temporary files for input/output results.
  temp.settings <- base::tempfile(pattern = '', fileext = '.primer3settings')
  
  
  # Discover target and exclusion sites. ----
  
  # Discover the position(s) of the targets within the given DNA sequence.
  target.setting <- if(!is.numeric(target) & !is.null(target)){
    target.found <- Biostrings::matchPattern(pattern = target, subject = input.sequence)
    base::paste(base::paste(GenomicRanges::start(target.found), GenomicRanges::width(target.found), sep = ','), collapse = ' ')
  }else{
    paste(target, collapse = ',')
  }
  
  exclusion.setting <- if(!is.numeric(exclude) & !is.null(exclude)){
    exclude.found <- Biostrings::matchPattern(pattern = exclude, subject = input.sequence)
    base::paste(base::paste(GenomicRanges::start(exclude.found), GenomicRanges::width(exclude.found), sep = ','), collapse = ' ')
  }else{
    paste(exclude, collapse = ',')
  }
  
  
  # Write settings file. ----------------------------------------------------
  
  settings <- base::list('SEQUENCE_ID' = paste0('Temp', stats::runif(1)),
                         'SEQUENCE_TEMPLATE' = as.character(input.sequence),
                         'PRIMER_TASK' = 'generic',
                         'PRIMER_PICK_LEFT_PRIMER' = 1,
                         'PRIMER_PICK_RIGHT_PRIMER' = 1,
                         'PRIMER_PICK_INTERNAL_OLIGO' = as.integer(internalOligo),
                         'PRIMER_PAIR_MAX_DIFF_TM' = 3,
                         'PRIMER_EXPLAIN_FLAG' = 0,
                         'SEQUENCE_TARGET' = target.setting,
                         'SEQUENCE_EXCLUDED_REGION' =  exclusion.setting,
                         'PRIMER_MIN_TM' = primer.tm[1],
                         'PRIMER_OPT_TM' = primer.tm[2],
                         'PRIMER_MAX_TM'= primer.tm[3],
                         'PRIMER_MIN_SIZE' = primer.size[1],
                         'PRIMER_OPT_SIZE' = primer.size[2],
                         'PRIMER_MAX_SIZE'= primer.size[3],
                         'PRIMER_PRODUCT_SIZE_RANGE' = paste(product.sizes, collapse = ' '),
                         'PRIMER_MISPRIMING_LIBRARY' = mismatching.library,
                         'PRIMER_NUM_RETURN' = nPrimer,
                         'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT' = 1,
                         'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT' = 0,
                         'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' = primer3_config
  )
  
  # Add additional user primer3 parameters.
  settings <- c(settings, primer3.args)
  
  # Remove empty parameters.
  settings <- settings[!base::sapply(settings, base::is.null)]
  
  # Create a string with \n between parameters.
  settings.concat <- base::paste(base::lapply(base::names(settings), function(x) base::sprintf('%s=%s', x, settings[[x]]) ), collapse = '\n')
  
  # Write the settings to the temporary file. (Also add final closing =)
  futile.logger::flog.debug('PRIMA - Writing primer3 setting to: %s', temp.settings)
  base::write(paste0(settings.concat, '\n='), file = temp.settings)
  
  # Make primer3 command.
  command.primer3 <- sprintf('%s %s -p3_settings_file %s',  base::Sys.which(primer3_core.binary), temp.settings, primer3_defaultsettings)
  
  # Run the command.
  command.output <- PRIMA::runSystemCommand(command.primer3)
  
  # Remove temporary files.
  unlink(temp.settings)
  
  # Check if command failed.
  if(command.output$exitStatus != 0){
    stop('Could not run primer3: ', command.output$stderr)
  }
  
  
  # Convert primer3 output to data.table ------------------------------------
  
  futile.logger::flog.trace('performPrimer3 - Converting primer3 output to data.table.')
  primer3.STDOUT <- data.table::data.table(utils::read.table(base::textConnection(command.output$stdout), sep = '=', strip.white = T, stringsAsFactors = F))
  
  primer3.results <- data.table::data.table(t(primer3.STDOUT$V2))
  base::colnames(primer3.results) <- t(primer3.STDOUT$V1)
  
  
  # Check for primer3 errors ------------------------------------------------
  
  if('PRIMER_ERROR' %in% colnames(primer3.results)){
    stop('Error in primer3: ', primer3.results$PRIMER_ERROR)
  }
  
  
  # Return results ----------------------------------------------------------
  
  return(primer3.results)
  
}