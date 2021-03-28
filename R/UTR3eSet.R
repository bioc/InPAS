#' UTR3eSet-class
#'  
#' An object of class [UTR3eSet-class] represents the results of 3' 
#' UTR usage
#' @name UTR3eSet-class
#' @slot usage GRanges. 
#' @slot PDUI matrix. 
#' @slot PDUI.log2 matrix. 
#' @slot short matrix. 
#' @slot long matrix. 
#' @slot signals list. 
#' @slot testRes matrix. 
#' 
#' @return
#' @import methods
#' @export
#'

setClass("UTR3eSet",
         representation(usage="GRanges", PDUI="matrix",
                        PDUI.log2="matrix", 
                        short="matrix",
                        long="matrix",
                        signals="list",
                        testRes="matrix"))

#' Constructor 
#' 
#' @name new
#' @param ... data for constructing an object of [UTR3eSet-class]
#' @rdname UTR3eSet-class
#' @return an object of [UTR3eSet-class]
#' 
UTR3eSet <- function(...){
    new("UTR3eSet", ...)
}

#' Accessor function for UTR3eSet class
#' 
#' @name $
#' @param UTR3eSet an object of [UTR3eSet-class]
#' @param x an object of [UTR3eSet-class]
#' @param name slot name
#' @rdname UTR3eSet-class
#' @return a slot value of an object of UTR3eSet class
#' @exportMethod $
#' @aliases $,UTR3eSet-method
#'

setMethod("$", "UTR3eSet", function(x, name) slot(x, name))

#' Set slot information of UTR3eSet
#' 
#' @name $<-
#' @param x an object of [UTR3eSet-class]
#' @param name slot name
#' @param value a value to be assigned
#' @rdname UTR3eSet-class
#' @exportMethod $<-
#' @aliases $<-,UTR3eSet-method
#' 
setReplaceMethod("$", "UTR3eSet",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

#' convert an object of UTR3eSet class to an ExpressionSet object
#' 
#' @name as
#' @param from an object of [UTR3eSet-class]
#' @param to an ExpressionSet object
#' @rdname UTR3eSet-class
#' @export
#' 
setAs(from="UTR3eSet", to="ExpressionSet", function(from){
    ExpressionSet(assayData=from$PDUI.log2,
                  featureData=AnnotatedDataFrame(as.data.frame(from$usage)))
})

#' convert an object of UTR3eSet class to a GRanges object
#' 
#' @name as
#' @param from an object of [UTR3eSet-class]
#' @param to an object of GRanges
#' @rdname UTR3eSet-class
#' @export

setAs(from="UTR3eSet", to="GRanges", function(from){
    n <- from$usage$transcript
    short <- from$short
    long <- from$long
    colnames(short) <- paste(colnames(short), "short.form.usage", sep="_")
    colnames(long) <- paste(colnames(long), "long.form.usage", sep="_")
    PDUI <- from$PDUI
    colnames(PDUI) <- paste(colnames(PDUI), "PDUI", sep="_")
    
    r <- from$usage
    mcol_r <- as.data.frame(mcols(r))
    
    if(length(from$testRes)!=0){
        test <- from$testRes
        mcols(r) <- cbind(mcol_r, 
                          short[match(n, rownames(short)), , drop=FALSE], 
                          long[match(n, rownames(long)), , drop=FALSE],
                          PDUI[match(n, rownames(short)), , drop=FALSE],
                          test[match(n, rownames(test)), , drop=FALSE])
    }else{
        mcols(r) <- cbind(mcol_r, 
                          short[match(n, rownames(short)), , drop=FALSE], 
                          long[match(n, rownames(long)), , drop=FALSE],
                          PDUI[match(n, rownames(short)), , drop=FALSE])
    }
    r
})

#' print method for UTR3eSet-class 
#' 
#' @name show
#' @param object an object 
#' @rdname UTR3eSet-class
#' @exportMethod show
#' @aliases show,UTR3eSet-method
setMethod("show", "UTR3eSet", function(object){
    cat("This is an object of UTR3eSet\n", "slot usage:\n")
    show(object$usage)
    if(length(object$PDUI)>0){
        cat("slot PDUI: matrix object with ", nrow(object$PDUI), 
            " rows ", ncol(object$PDUI), "columns\n")
        show(head(object$PDUI))
        cat("...\n")
        cat("slot PDUI.log2: matrix object with ", nrow(object$PDUI.log2), 
            " rows ", ncol(object$PDUI.log2), "columns\n")
        show(head(object$PDUI.log2))
            cat("...\n")
        cat("slot short: matrix object with ", nrow(object$short), 
            " rows ", ncol(object$short), "columns\n")
        show(head(object$short))
            cat("...\n")
        cat("slot long: matrix object with ", nrow(object$long), 
            " rows ", ncol(object$long), "columns\n")
        show(head(object$long))
            cat("...\n")
    }
    if(length(object$signals)>0){
        cat("slot signals: length ", length(object$signals), "\n")
    }
    if(length(object$testRes)>0){
        cat("slot testRes: matrix object with ", nrow(object$testRes), 
            " rows ", ncol(object$testRes), "columns\n")
        show(head(object$testRes))
        cat("...\n")
    }
})

