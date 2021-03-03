##############################################
# utility functions
##############################################

#' Progress function
#'
#' Non exported function, plots the prgress of an iterative procedure using "[=>   ]", "[==> ]", etc.
#'
#' @param p Value
#' @return NULL
loader <- function(p)
{
  if (p==0) cat("0%                       50%                     100%\n")
 str <- paste0(rep(c("\r[", "=", ">", " ", "]"), c(1, floor(p*50), 1, 50 - floor(p*50), 1)), collapse = "")
 cat(str)
 utils::flush.console()
 if (floor(p) == 1) cat("\n")
}

fGetValsScan <- function(Nums,RPs){
  
  #Nums - character vector with values or the name of refererence points to be substituted
  #RPs - named vector of reference points 

  if (length(Nums)==1){
    if (is.na(Nums)){
      return(NA)
    }
  }

  #extract the numeric ones
  ret <- as.numeric(Nums[!is.na(suppressWarnings(as.numeric(Nums)))])
  
  #now the specified reference points
  if (sum(is.na(suppressWarnings(as.numeric(Nums))))>0) {
    ret <- c(ret,as.numeric(RPs[Nums[is.na(suppressWarnings(as.numeric(Nums)))]]))
  }
  
  #add a warning if there are unmatched RPs
  
  sort(ret)
  
}

##' Write ICES/CEFAS data file from matrix 
##' @param x a matrix where rownames are taken as years and colnames are taken as ages
##' @param fileout file name or connection
##' @param ... Arguments to be passed to write 
##' @details
##' 
##' Takes the data and writes them in the ICES/CEFAS format. It is assumed that rows represent consecutive years and cols consecutive ages 
##' 
##' @export
write.ices <- function(x, fileout, ...){
  top <- paste0(fileout, " auto written\n1 2\n", paste0(range(as.integer(rownames(x))), collapse="  "),"\n",paste0(range(as.integer(colnames(x))), collapse="  "), "\n1")    
  write(top,fileout,...)
  write(t(x),fileout,ncolumns=ncol(x),append=TRUE,sep="  \t",...)
}

##' Write all data files from a list as created by 'setup.sam.data'  
##' @param dat A list as created by 'setup.sam.data'
##' @param dir Directory where the files are written  
##' @details
##' 
##' Write all data files from a list as created by 'setup.sam.data'
##' 
##' @export
write.data.files<-function(dat, dir="."){
  od <- setwd(dir)
  write.ices(dat$data$catchMeanWeight, "cw.dat")
  write.ices(dat$data$disMeanWeight, "dw.dat")
  write.ices(dat$data$landMeanWeight, "lw.dat")
  write.ices(dat$data$landFrac, "lf.dat")  
  write.ices(dat$data$propMat, "mo.dat")    
  write.ices(dat$data$stockMeanWeight, "sw.dat")
  write.ices(dat$data$propF, "pf.dat")
  write.ices(dat$data$propM, "pm.dat")
  write.ices(dat$data$natMor, "nm.dat")
  # fit <- list(data=dat)
  write.ices(getFleet(dat,1), "cn.dat")
  write.surveys(fit, "survey.dat")
  setwd(od)
}

##' Extract a fleet from a fitted object 
##' @param fit A fitted object as returned from sam.fit
##' @param fleet The number of the fleet 
##' @export
getFleet <- function(fit, fleet){
  fidx <- fit$data$aux[,"fleet"]==fleet
  aux <- fit$data$aux[fidx,]
  logobs <- fit$data$logobs[fidx ]
  .goget <- function(y, a) {
    ret <- exp(logobs[aux[, "year"] == y & aux[, "age"] == a])
    ifelse(length(ret) == 0, 0, ret)
  }
  yr <- min(aux[,"year"]):max(aux[,"year"])
  ar <- min(aux[,"age"]):max(aux[,"age"])
  tmp <- outer(yr, ar, Vectorize(.goget))
  dimnames(tmp)[[1]] <- yr
  dimnames(tmp)[[2]] <- ar
  return(tmp)
}

##' Write surveys in ICES/CEFAS data file from a model object  
##' @param fit A fitted object as returned from sam.fit
##' @param fileout file name or connection
##' @param ... Arguments to be passed to write 
##' @details
##' 
##' Takes the survey data from the fitted object and writes them in the ICES/CEFAS format. 
##' 
##' @export
write.surveys <- function(fit,fileout,...){
  sidx <- which(fit$data$fleetTypes%in%c(2,3,4))
  top <- paste0(fileout," auto written\n", 100+length(sidx))
  write(top,fileout,...)
  for(s in sidx){
    write(paste0(attr(fit$data, "fleetNames")[s]), fileout, append=TRUE, ...)
    S <- getFleet(fit,s)
    yr <- range(as.integer(rownames(S)))
    ar <- range(as.integer(colnames(S)))
    write(paste0(yr[1], " ", yr[2]), fileout, append=TRUE, ...)
    st <- fit$data$sampleTimes[s]
    write(paste0(1," ", 1, " ", st, " ", st), fileout, append=TRUE, ...)
    write(paste0(ar[1], " ", ar[2]), fileout, append=TRUE, ...)
    x <- cbind(1,S)
    write(t(x),fileout,ncolumns=ncol(x),append=TRUE,sep="  \t",...)
  }
}

# Load RData function
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# lowcase function
lowcase <- function(df) {
  names(df) <- tolower(names(df)) %>% gsub("\\?|\\s+|\\.+|_+|\\(|\\)","",.) 
  df
}

# Theme publication (for ggplot)
theme_publication <- function(base_size=14, base_family="Helvetica") {
  # library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title       = element_text(face = "bold",size = rel(0.8), hjust = 0.0),
            plot.margin      = unit(c(10,5,5,5),"mm"),
            plot.background  = element_rect(colour = NA),
            text             = element_text(),
            axis.title       = element_text(face = "bold",size = rel(1)),
            axis.title.y     = element_text(angle=90,vjust =2),
            axis.title.x     = element_text(vjust = -0.2),
            axis.text        = element_text(), 
            axis.line        = element_line(colour="black"),
            axis.ticks       = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(colour="black" , size=0.1),
            panel.background = element_rect(colour = NA),
            strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text       = element_text(face="bold"),
            legend.key       = element_rect(colour = NA),
            legend.position  = "bottom",
            legend.direction = "horizontal",
            legend.key.size  = unit(0.2, "cm"),
            legend.spacing   = unit(0, "cm"),  # updated from legend.margin which is deprecated
            legend.title     = element_text(face="italic", size=rel(0.8))
    ))
}

