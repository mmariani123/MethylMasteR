#!/usr/bin/env Rscript

##Michael Mariani PhD Dartmouth College 2022

##################### The MethylMasteR class #################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

##check out http://adv-r.had.co.nz/S4.html
##For more S4 class info

setClass("MethylMaster",
         representation(Sample_ID  = "character",
                        chrom      = "character",
                        loc.start  = "integer",
                        loc.end    = "integer",
                        num.mark   = "integer",
                        bstat      = "numeric",
                        seg.mean   = "numeric",
                        seg.median = "numeric",
                        pval       = "numeric",
                        state      = "integer",
                        treatment  = "character",
                        method     = "character",
                        sub.method = "character"),
         prototype(Sample_ID = NA_character_,
                   loc.start = NA_integer_),
         contains=NULL ##Inherits from
)

mmo <- new("MethylMaster",
           Sample_ID  = "test_sample",
           chrom      = "1",
           loc.start  = 1L,
           loc.end    = 1e6L,
           num.mark   = 1e3L,
           bstat      = 0.04,
           seg.mean   = 0.31,
           seg.median = 0.25,
           pval       = 0.04,
           state      = 3L,
           treatment  = "tumor",
           method     = "custom",
           sub.method = "")

##Note that Unlike S3, S4 checks that
##all of the slots have the correct type:

##To access slots of an S4 object you use @, not $:

mmo@Sample_ID

##Or if you have a character string giving
##a slot name, you use the slot function:

slot(mmo, "Sample_ID")

##getSlots will return a description of all the slots of a class:

getSlots("MethylMaster")
#        name         age
# "character"   "numeric"

##https://stackoverflow.com/questions/5411919/error-handling-with-s4-classes

##"if the assignment fails,
##I'd return an error instead of a warning.
##A warning tells you that the process went through,
##but might give an unexpected result. In your case,
##the process is aborted :

##setReplaceMethod(f="setInd",signature="foo",
##def=function(object,value){
##   if(!is.numeric(value))
##     stop("Foobar")
##
##   object@ind <- value
##   return(object)}
##)

##"using stop allows you to use tryCatch() or try() constructs.
##See the relevant help pages for more information. Eg :

##tryCatch(setInd(thisFoo)<-"A",error=function(e){print("Hello")})

##X <- try(setInd(thisFoo) <- "A")
##Error in `setInd<-`(`*tmp*`, value = "A") : Foobar
##if(is(X,"try-error")) setInd(thisFoo) <- 5
##thisFoo
##An object of class "foo"
##Slot "ind":
##[1] 5

##"If you really need to work with warnings,
##look at withCallingHandlers. Using your original code :"

##withCallingHandlers({setInd(thisFoo)<-"A"},
##     warning = function(w) {print("Hello")})
##[1] "Hello"
##Warning message:
##In `setInd<-`(`*tmp*`, value = "A") : Foobar

##"Mind you, this is less straightforward to use than
##the abovementioned options using errors.
