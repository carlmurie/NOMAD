#' Apply a highly scaleable linear model normalization on Mass Spectrometry data
#' 
#' Normalize MS data quickly and validly
#' @docType package
#' @name NOMAD
#' 
#' @import dplR
#' @import stringr
#' @import stats
#' @import graphics
NULL

#' Check whether logging the data is necessary. y and x are the same inputs
#' as nomadNormalization
#' Check whether logging the data is necessary
#'
#' @param y \code{vector} dependant variable (vector of real numbers)
#' @param x \code{matrix} independent or explanatory variables (matrix). x must include columns
#'     for Peptide, Protein, iTRAQ, and Run
#' @param meth \code{string} manner of calculating center of peptide abundances per protein
#'        - choices: TukeyBW, Median, or Mean
#' @param rawTrim \code{real} Quantile used to reduce scale for plotting raw data. Will plot
#'           raw means up to the specified quantile (default=1)
#' @return None
#' 
#' @examples 
#' data(BalfPeptides)
#' nomadCheckLogTransform(y=BalfPeptides$Abundance, x=BalfPeptides)
#' 
#' @export
nomadCheckLogTransform <- function(y, x, meth="TukeyBW", rawTrim=1) {
  
  ## remove itraqCorrectionIndex if it exists
  x$iTRAQCorrectionIndex <- NULL
  
  ## peptide level
  pMeans <- tapply(y, x$Peptide, mean, na.rm=TRUE)
  pSD  <- tapply(y, x$Peptide, sd, na.rm=TRUE)
  pMeanslog <- tapply(log2(y), x$Peptide, mean, na.rm=TRUE)
  pSDlog  <- tapply(log2(y), x$Peptide, sd, na.rm=TRUE)
  
  ## protein level
  rawProt <- nomadAssembleProteins(y,x, method=meth)$scores
  logProt <- nomadAssembleProteins(log2(y),x, method=meth)$scores
  
  protmeans <- apply(rawProt, 1, mean, na.rm=TRUE)
  protsds <-  apply(rawProt, 1, sd, na.rm=TRUE)
  protmeanslog <- apply(logProt, 1, mean, na.rm=TRUE)
  protsdslog <-  apply(logProt, 1, sd, na.rm=TRUE)
  
  ## retain only finite values - stops warning message about removing
  ## non-finfite values. sd may produce NA if there is only one value
  ## in protsdslog
  notInfInd <- is.finite(pMeanslog)
  notInfIndProt <-  is.finite(protmeanslog) & !is.na(protsdslog)
  
  ind1 <- pMeans <= quantile(pMeans, rawTrim, na.rm=TRUE)
  ind2 <- protmeans <= quantile(protmeans, rawTrim, na.rm=TRUE)
  
  ## plot results
  par(mfrow=c(2,2))
  smoothScatter(pMeans[ind1], pSD[ind1], xlab="peptide mean", ylab="peptide sd", main="Peptide mean vs sd: Raw")
  lines(supsmu(pMeans[ind1], pSD[ind1], bass=5), col=8)
  
  smoothScatter(pMeanslog[notInfInd], pSDlog[notInfInd], xlab="peptide mean", ylab="peptide sd", main="Peptide mean vs sd: Log")
  lines(supsmu(pMeanslog[notInfInd], pSDlog[notInfInd], bass=5), col=8)
  
  smoothScatter(protmeans[ind2], protsds[ind2], xlab="protein mean", ylab="protein sd", main="Protein mean vs sd: Raw")
  lines(supsmu(protmeans[ind2], protsds[ind2], bass=5), col=8)
  
  smoothScatter(protmeanslog[notInfIndProt], protsdslog[notInfIndProt], xlab="protein mean", ylab="protein sd", main="Protein mean vs sd: Log")
  lines(supsmu(protmeanslog[notInfIndProt], protsdslog[notInfIndProt], bass=5), col=8)
  
} # nomadCheckLog


#' Check level of bias based on batch effect.
#'
#' @param dat \code{matrix} Rows are protein, columns are samples (each iTRAQ sample 
#'       per day. The same format as the "scores" element of the output of
#'       assembleProteins with format=bySample. Missing data must be shown
#'       as NA.
#' @param fact \code{vector} numeric or character vector of labels identifying batches.
#' Generally a batch is a single MS run of 4 or 8 iTRAQ channels.
#' @param numNAs \code{integer} Number of missing data points allowed. Defaults to number of
#'          samples in data
#' @param label \code{string} Text for title of graph
#' @param yMax \code{real} Upper count range for plotting histograms.
#' @param showSampleSize \code{boolean} If TRUE add text box showing number of samples.
#' 
#' @return None
#' 
#' @examples 
#' data(BalfPeptides)
#' 
#' ## assemble the protein scores without normalization
#' scores <- nomadAssembleProteins(BalfPeptides$Abundance, BalfPeptides)
#' runInd <- rep(1:3, each=8)  ## each integer is a single batch (3 batches total)
#' nomadCheckBias(scores$scores, fact=runInd)
#' @export
nomadCheckBias <- function(dat, fact, numNAs=dim(dat)[[2]], label=NULL, yMax=NULL, showSampleSize=TRUE) {
  
  ## Apply one-way anova to data.
  theTest <- function(x, indy, nas=numNAs) {
    tmp <- data.frame(y=x, exp=indy)
    
    ## if too much missing data then return NA
    if(sum(is.na(x)) > nas) {
      return(NA)
    }
    
    res <- lm(y~indy, data=tmp)
    return(anova(res)$"Pr(>F)"[1])
  } ## end theTest
  
  allRes <- apply(dat, 1, theTest, indy=fact)
  
  maxCount <- max(hist(allRes, nclass=20, plot=FALSE)$counts)
  if(!is.null(yMax)) {
    maxCount <- yMax
  }
  hist(allRes, nclass=20, main=paste(label, "bias"), xlab="p-value",
       ylab="count", ylim=c(0,maxCount))
  
  if(showSampleSize) {
    n <- length(allRes)
    legend("topright", legend=c(paste("n=", n)), bty="n")
  }
  
} ## end nomadCheckBias





#' Apply anova model to mass spec data. Use residuals of anova model as
#' final scores.
#'
#' @param y \code{real} Dependant variable (vector of positive real numbers)
#' @param x \code{matrix} Independent or explanatory variables (matrix). In this context they
#'     will be "Peptide","Protein", "Run", and "iTRAQ" identifiers. An
#'     optional "iTRAQCorrectionIndex" column is necessary if iTRAQ isotope
#'     correction is desired.
#' @param factors \code{vector} Default is to use the 4 main factors (Peptide, Protein,
#'            Run, and iTRAQ)  and the three interactions between
#'            Run and the other three main effects. The user can
#'            define their own factors. The "factors" list must match the
#'            data.frame labels of x.
#' @param doRobust \code{logical} If TRUE use medians of levels instead of the default means.
#'            Default is to use median.
#' @param doLog \code{logical} If TRUE then y will be logged (log 2). Default is to log.
#' @param doiTRAQCorrection  \code{logical} If true then orrect peptide abundance scores based on
#'                     their iTRAQ labelling. Default is TRUE. iTRAQ
#'                    correction values are hardcoded for 8 iTRAQs. The
#'                    iTRAQCorrectionTable must be provided for iTRAQs that
#'                    do not number 8.
#' @param iTRAQCorrectionTable \code{matrix} Matrix of iTRAQ corrections. If iTRAQ equals 8
#'                        then an internal correction table is used. If a
#'                        table is provided when iTRAQ equals 8 then the
#'                        user supplied table will be used. If the number
#'                        of iTRAQs is not eight than a table must be
#'                        supplied by the user.
#'                        
#' @return \code{List} elements are \code{y}: normalized scores after anova normalization,
#' \code{x}: design matrix, \code{removedData}: index of rows removed from input data,
#' \code{call}: function call. 
#'
#' @examples 
#' data(BalfPeptides)
#' 
#' ret <- nomadNormalization(y=BalfPeptides$Abundance, x=BalfPeptides)
#'
#' @export
#'
nomadNormalization <- function(y, x, factors=NULL, doRobust=TRUE,
                               doLog=TRUE, doiTRAQCorrection=FALSE,
                               iTRAQCorrectionTable=NULL) {
  
  ## error check parameters
  if (missing(y))
    stop("Need to specify dependent variable: y.")
  if (missing(x))
    stop("Need to specify independent variable: x.")
  
  if(!is.vector(y) | !is.numeric(y)) {
    stop("y must be a numeric vector")      
  }
  
  if(sum(y < 0, na.rm=TRUE) > 0) {
    stop("y contains negative numbers")      
  }
  
  if(!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  
  if(length(y) != dim(x)[[1]]) {
    outChar <- paste("Length of y (", length(y),") does not equal length of x (", dim(x)[[1]], ")", sep="")
    stop(outChar)
  }
  
  if(is.null(x$Protein)) {
    stop(paste("x is missing Protein column"))
  }
  
  if(is.null(x$Peptide)) {
    stop(paste("x is missing Peptide column"))
  }
  
  if(is.null(x$Run)) {
    stop(paste("x is missing Run column"))
  }
  
  if(is.null(x$iTRAQ)) {
    stop(paste("x is missing iTRAQ column"))
  }
  
  if(doiTRAQCorrection==TRUE) {
    if(is.null(x$iTRAQCorrectionIndex)) {
      stop("x is missing iTRAQCorrectionIndex column")
    }
  }
  
  ## assemble factors for normalization including interactions
  if(is.null(factors)) {
    factorList <- list("Protein", "Peptide", "Run", "iTRAQ", 
                       c("Run", "Protein"),
                       c("Run", "Peptide"),
                       c("Run", "iTRAQ"))
  } else {
    factorList <- factors
  }
  
  ## check that factors exist in data
  for(i in 1:length(factorList)) {
    aFac <- factorList[[i]]
    for(j in 1:length(aFac)) {
      thisFac <- aFac[j]
      if(is.null(x[[thisFac]])) {
        stop(paste("x is missing", thisFac,  "column"))
      }
    }
  }
  
  # get number of iTRAQs used in assay
  numiTRAQs <- length(unique(x$iTRAQ))
  
  ## find missing data - full row will be removed for any missing data
  naInd <- is.na(y) | is.na(as.character(x$Peptide)) |
    is.na(as.character(x$Protein)) | is.na(as.character(x$Run)) |
    is.na(as.character(x$iTRAQ))
  
  ## remove rows with missing data
  #x <- x[!naInd,]
  iterResids <- y[!naInd]            ## initial set of data points
  newx <- x[!naInd,]
  len <- length(iterResids)          ## number of data points
  
  
  ## correct peptide abundances based on iTRAQ label
  if(doiTRAQCorrection) {
    
    #ordy <- order(as.factor(paste(newx$iTRAQCorrectionIndex, newx$iTRAQ, sep="_")))
    fac <- newx$iTRAQCorrectionIndex
    #newx <- newx[ordy,]
    #newy <- tmpy <- iterResids[ordy]
    newy <- tmpy <- iterResids
    itraqy <- data.frame(newy, newx$iTRAQ)
    
    if(is.null(iTRAQCorrectionTable)) {
      iTRAQTab <- NULL 
    } else {
      iTRAQTab <-  iTRAQCorrectionTable
    }
    
    cat("doing iTRAQCorrection with", numiTRAQs, "iTRAQs\n")
    if(numiTRAQs == 8) {
      split(tmpy, fac) <- lapply(split(itraqy, fac),
                                 solveiTRAQCorrection8,
                                 iTRAQTable=iTRAQTab)
    } else if(numiTRAQs == 4) {
      split(tmpy, fac) <- lapply(split(itraqy, fac),
                                 solveiTRAQCorrection4,
                                 iTRAQTable=iTRAQTab)
    } else {
      if(is.null(iTRAQTab)) {
        stop("User needs to input iTRAQCorrectionTable to apply iTRAQCorrection")
      }
      if(numiTRAQs != dim(iTRAQTab)[[1]] &
         numiTRAQs != dim(iTRAQTab)[[1]]) {
        stop("dimension of iTRAQCorrectionTable must be equal to number of iTRAQs")
      }
      
      split(tmpy, fac) <- lapply(split(itraqy, fac),
                                 solveiTRAQCorrection,
                                 iTRAQTable=iTRAQTab)
    }
    
    ## replace negative scores with NA
    tmpy[ tmpy < 0 ] <- NA
    iterResids <- tmpy
    
  } ## end doiTRAQCorrection
  
  ## log data if requested. Zero scores may produce a negative infinity
  ## so set them to minimum observed non-infinity score.
  if(doLog) {
    iterResids <- log2(iterResids)
    iterResids[iterResids==-Inf] <- min(iterResids[iterResids!=-Inf])
    
  } ## end doLog
  
  cat("Running normalization with ", len, " number of data points\n")
  
  ## iterate through single and interaction factors: subtract
  ## interaction level means
  for(i in 1:length(factorList)) {
    
    thisFactor <- factorList[[i]]
    cat("Normalizing for factor: ", thisFactor, "\n")
    
    ## get labels of levels for this factor(s)
    newInd <- as.character(newx[, thisFactor[1]])
    if(length(thisFactor) > 1) {
      newInd <- paste(newInd, newx[,thisFactor[2]], sep="-")
    }
    
    ## calculate mean of each level. simplify ensures that tapply returns
    ## a list so we can track level IDs
    temp <- iterResids
    if(doRobust) {
      split(temp, newInd) <- lapply(split(iterResids, newInd),
                                    median, na.rm=TRUE)
      
    } else {
      split(temp, newInd) <- lapply(split(iterResids, newInd),
                                    mean, na.rm=TRUE)
    }
    
    iterResids <- iterResids - temp
    
  } # end for i
  
  ## remove iTRAQCorrectionIndex column, we don't need it anymore
  newx$iTRAQCorrectionIndex <- NULL
  
  ## generate function call parameters
  newFac <- paste(factorList[[1]], collapse="_")
  for(i in 2:length(factorList)) {
    newFac <- paste(newFac, paste(factorList[[i]], collapse="_"), sep=", ")
  }
  
  itraqTable <- NULL
  if(doiTRAQCorrection & numiTRAQs==8) {
    itraqTable <- "hardcoded"
  } 
  
  model <- c(paste("factors =", newFac), paste("doRobust =", doRobust),
             paste("doLog = ", doLog),
             paste("doiTRAQCorrection =", doiTRAQCorrection),
             paste("iTRAQCorrectionTable =", itraqTable)
  )
  
  return(list(y=as.vector(iterResids, mode="numeric"), x=newx, removedData=which(naInd), call=model))
  
} ## end nomadNormalization




#'
#' Assemble Protein abundance measurements from Peptide abundance
#' measurements
#'
#' @param y \code{vector} Dependant variable (vector of real numbers)
#' @param x \code{vector} Independent or explanatory variables (matrix). In this context they
#'     will be Peptide, Protein, Run, and iTRAQ identifiers.
#'     Run and iTRAQ identifiers must be numeric
#' @param method \code{string} Method of calculating center of peptide abundances per protein
#'        - choices: TukeyBW, Median, or Mean. Default is TukeyBW
#' @param format \code{string} BySample - each column is an individual sample
#'          ByProtein - rows are proteins-columns are abundance,day and itraq
#' @param combineDupPeptides \code{string} If TRUE then average (based on method) duplicated
#'                  peptides for the same protein before calculating protein
#'                  abundance
#' @param standardizeAbundance \code{logical} If TRUE then transform each sample (experiment and
#'                       iTRAQ) with a robust Z adjustment. The default is
#'                       FALSE.
#' @param mySep \code{logical} Character string that must not exist in any factor in "x".
#'
#' @return \code{List} with elements \code{scores} matrix of normalized protein scores. 
#' each row is a protein and each column is a single sample (run and iTRAQ)
#'
#' @examples 
#' data(BalfPeptides)
#' 
#' ret <- nomadNormalization(y=BalfPeptides$Abundance, x=BalfPeptides)
#' scores <- nomadAssembleProteins(ret$y, ret$x)
#'
#' @export
nomadAssembleProteins <- function(y, x, method="TukeyBW", format="BySample", combineDupPeptides=FALSE, standardizeAbundance=FALSE, mySep="&_MYSEP_&") {
  
  ############### error check parameters ###########################
  
  if (missing(y))
    stop("Need to specify dependent variable: y.")
  if (missing(x))
    stop("Need to specify independent variable: x.")
  
  if(!is.vector(y) | !is.numeric(y)) {
    stop("y must be a numeric vector")      
  }
  
  if(!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  
  if(is.null(x$Peptide)) {
    stop("x is missing Peptide column")
  }
  
  if(is.null(x$Protein)) {
    stop("x is missing Protein column")
  }
  
  if(is.null(x$Run)) {
    stop("x is missing Run column")
  }
  
  if(is.null(x$iTRAQ)) {
    stop("x is missing iTRAQcolumn")
  }
  
  if(length(method) != 1) {
    stop("method must be one of TukeyBW, Median, or Mean")
  }
  
  methy <- intersect(method, c("TukeyBW", "Median", "Mean"))
  if(length(methy) != 1) {
    stop("method must be one of TukeyBW, Median, or Mean")
  }
  
  methy <- intersect(format, c("BySample", "ByProtein"))
  if(length(methy) != 1) {
    stop("format variable must be one of BySample or ByProtein" )
  }
  
  ## find missing data - full row will be removed for any missing data
  naInd <- is.na(y) | is.na(as.character(x$Peptide)) |
    is.na(as.character(x$Protein)) | is.na(as.character(x$Run)) |
    is.na(as.character(x$iTRAQ))
  
  ## check whether crazy biologist has managed to use my separator in
  ## their annotations.
  if(length(grep(mySep, as.character(x$Protein))) != 0 | 
     length(grep(mySep, as.character(x$Run))) |
     length(grep(mySep, as.character(x$iTRAQ)))) {
    stop("mySep variable " , mySep, " is contained in design matrix x. Please redefine variable mySep to a character string not contained in design matrix x")
  }
  
  ## calculate estimate of center of peptide abundance (x)
  centre <- function(x, type) {
    outie <- switch(type,
                    Mean = mean(x, na.rm=TRUE),
                    Median = median(x, na.rm=TRUE),
                    TukeyBW = tbrm(x))
    
    if(is.null(outie)) {outie <- NA}
    return(outie)
  }
  
  ############## now assemble proteins #############################
  
  allPepsTmp <- as.matrix(cbind(as.character(y), as.character(x$Peptide), as.character(x$Protein), as.character(x$Run), as.character(x$iTRAQ)))
  
  ## calculate centre of duplicated peptides for same protein before
  ## calculating protein abundance
  if(combineDupPeptides) {
    
    ## get unique label for each peptide by protein, experiment and itraq
    pepIDs <- paste(x$Peptide, x$Protein, x$Run, x$iTRAQ, sep=mySep)
    
    ## calculate center of scores for each potential duplicate peptide
    ## (by protein, experiment and itraq)
    pepScores <- tapply(y, pepIDs, centre, type=method, simplify=FALSE)
    
    ## unpack protein, experiment, and itraq labels - attach to protein
    ## scores
    nameys <- names(pepScores)
    allPepsTmp <- cbind(unlist(pepScores), matrix(unlist(str_split(nameys, mySep)), ncol=4, byrow=TRUE))
    
  } ## end sumDupPeptides
  
  dimnames(allPepsTmp) <- list(NULL, c("Abundance", "Peptide", "Protein", "Run", "iTRAQ"))
  
  ## get unique label for each protein by experiment and itraq
  protIDs <- paste(allPepsTmp[,"Protein"], allPepsTmp[,"Run"],
                   allPepsTmp[,"iTRAQ"], sep=mySep)
  
  ## calculate center of scores for each protein (by experiment and itraq)
  protScores <- tapply(as.numeric(allPepsTmp[,"Abundance"]), protIDs, centre, type=method, simplify=FALSE)
  
  ## unpack protein, experiment, and itraq labels - attach to protein
  ## scores
  nameys <- names(protScores)
  allProtsTmp <- cbind(unlist(protScores), matrix(unlist(str_split(nameys, mySep)), ncol=3, byrow=TRUE))
  
  dimnames(allProtsTmp) <- list(NULL, c("Abundance", "Protein", "Run", "iTRAQ"))
  
  ## order by protein, experiment, itraq
  ordInd <- order(allProtsTmp[,2], as.numeric(allProtsTmp[,3]), as.numeric(allProtsTmp[,4]))
  allProts <- allProtsTmp[ordInd,]
  
  if(standardizeAbundance) {
    myRobZ <- function(x) { return( (x - median(x, na.rm=TRUE))/mad(x, na.rm=TRUE) ) }
    
    standInd <- paste(allProts[,3], allProts[,4], sep="_")
    scores <- tmpy <- as.numeric(allProts[,1])
    
    split(tmpy, standInd) <- lapply(split(scores, standInd), myRobZ)
    
    allProts[,1] <- tmpy
    
    ## testing
    if(FALSE) {
      scores2 <- scores3 <- raw <- NULL
      for(i in unique(allProts[,3])) {
        for(j in unique(allProts[,4])) {
          indy <- allProts[,3] == i & allProts[,4] == j
          stag <- paste(i, j, sep="_")
          sind <- standInd == stag
          stmp <- as.numeric(allProts[sind,1])
          scores3 <- c(scores3, myRobZ(stmp))
          raw <- c(raw, stmp)
          cat(i, " ", j, " ", sum(indy), "\n")
          tmp <- myRobZ(as.numeric(allProts[indy,1]))
          scores2 <- c(scores2, tmp)
        }
      }
      cor(scores2, raw, use="complete.obs")
      cor(tmpy, scores, use="complete.obs")
      
    } ## emd if FALSE
    
  } ## end if standardizeAbundance
  
  ## construct function call
  model <- c(paste("method =", method), paste("format =", format),
             paste("combineDupPeptides =", combineDupPeptides),
             paste("standardizeAbundance =",  standardizeAbundance))
  
  
  ## choose which type of format to return 
  
  ## reformat data so gene IDs are rows and samples are columns
  if(identical(format, "BySample")) {
    genies <- sort(unique(allProts[,2]))
    days <- sort(unique(as.numeric(allProts[,3])))
    traqs <-  sort(unique(allProts[,4]))
    numSamps <- length(unique(days))*length(unique(traqs))
    newDat <- matrix(NA, ncol=numSamps, nrow=length(genies))
    row.names(newDat) <- genies
    
    count <- 0
    colNames <- NULL
    for(i in 1:length(days)) {
      for(j in 1:length(traqs)) {
        colNames <- c(colNames, paste("Run", i, "_iTRAQ", j, sep=""))
        ##print(i)
        count <- count + 1
        indy <- allProts[,3] == days[i] & allProts[,4] == traqs[j]
        gNames <- allProts[indy,2]
        gInd <- match(gNames, genies)
        tmp <- allProts[indy,1]
        newDat[gInd,count] <- as.numeric(tmp)
      } # end for j
    } # end for i
    
    colnames(newDat) <- colNames
    outie <- newDat
  } ## end if bySample
  
  if(identical(format, "ByProtein")) {
    outie <- allProts
  }
  
  return(list(scores=outie, call=model))
  
} ## end nomadAssembleProteins



#'
#' Correct peptide abundances based on iTRAQ specification. The iTRAQTable
#' below characterizes the components of the abundance scores.
#'
#' @param y \code{logical} A vector of 8 peptide abundances ordered by iTRAQ (1 to 8)
#' for a single protein
#' @param iTRAQTable \code{logical}  A 8x8 matrix of isotope corrections. Probably shouldn't modify
#'              unless you have a very good reason.
#'              
#' @return vector of peptide abundances after correcting for isotope bleed over.
#' 
solveiTRAQCorrection8 <- function(y, iTRAQTable=NULL) {
  
  itraqs <- 1:8
  nas <- setdiff(itraqs, y[,2])
  thisy <- rep(0, 8)
  thisy[ y[,2]] <- y[,1]
  
  if(length(thisy) != 8) {
    stop("solveiTRAQCorrection failed: length of y does not equal 8")
  }
  
  if(is.null(iTRAQTable)) {
    iTRAQTable <- rbind(
      c(.9287, 0.0689, .0024,  0,    0,     0,      0,     0),
      c(.0094,   0.93,   .059,  .0016, 0,     0,      0,     0),
      c(0,     .0188,  .9312, .049,  .001,  0,      0,     0),
      c(0,     0,      .0282, .9321, .039,  0.0007, 0,     0),
      c(0,     0,      .0006, .0377, .9329, .0288,  0,     0),
      c(0,     0,      0,     .0009, .0471, .9329, .0191, 0),
      c(0,     0,      0,     0,     .0014, .0566,  .9333, .0087),
      c(0,     0,      0,     0,     0,     .0027,  .0744, .9211)
    )
  }
  
  res <- solve(iTRAQTable, as.numeric(thisy))
  res[nas] <- NA
  return(res)
  
} # end solveiTRAQCorrection8


#'
#' Correct peptide abundances based on iTRAQ specification. The iTRAQTable
#' below characterizes the components of the abundance scores.
#'
#' @param y \code{vector} A vector of 4 peptide abundances ordered by iTRAQ (1 to 4)
#' @param iTRAQTable \code{matrix} An iTRAQTable - 4x4 matrix of isotope corrections. Probably shouldn't modify
#'              unless you have a very good reason.
#'              
#' @return vector of peptide abundances after correcting for isotope bleed over.
solveiTRAQCorrection4 <- function(y, iTRAQTable=NULL) {
  
  itraqs <- 1:4
  nas <- setdiff(itraqs, y[,2])
  thisy <- rep(0, 4)
  thisy[ y[,2]] <- y[,1]
  
  if(length(thisy) != 4) {
    stop("solveiTRAQCorrection failed: length of y does not equal 4")
  }
  
  if(is.null(iTRAQTable)) {
    iTRAQTable <- rbind(
      c(.9290, .0590, .0020, 0),
      c(.0200, .9230, .0560, .0010),
      c(0.00, .0300, .9240, .0450),
      c(0, .0010, .0400, .923))
  }
  
  res <- solve(iTRAQTable, as.numeric(thisy))
  res[nas] <- NA
  return(res)
} ## end solveiTRAQCorrection4


#'
#' Correct peptide abundances based on iTRAQ specification. The iTRAQTable,
#' below characterizes the components of the abundance scores.
#'
#' @param y \code{vector} A vector of n peptide abundances ordered by iTRAQ (n iTRAQ labels)
#' @param iTRAQTable \code{matrix} n by n matrix of isotope corrections
#' 
#' @return vector of peptide abundances after correcting for isotope bleed over.
solveiTRAQCorrection <- function(y, iTRAQTable=NULL) {
  
  len <- dim(iTRAQTable)[[1]]
  itraqs <- 1:len
  nas <- setdiff(itraqs, y[,2])
  thisy <- rep(0, len)
  thisy[ y[,2]] <- y[,1]
  
  if(length(thisy) != len) {
    stop("solveiTRAQCorrection failed: length of y does not equal dimension of iTRAQTable")
  }
  
  res <- solve(iTRAQTable, as.numeric(thisy))
  res[nas] <- NA
  return(res)
  
} # end solveiTRAQCorrection

