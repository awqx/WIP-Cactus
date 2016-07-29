read_cactus <- function(URL, destfile) {
        out <- tryCatch(
{
        # Just to highlight: if you want to use more than one 
        # R expression in the "try" part then you'll have to 
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression 
        # in case the "try" part was completed successfully
        
        message("Trying URL")
        #readLines(con = URL,warn = FALSE) # 1st approach
         Map(download.file, url = URL, destfile = destfile) #2nd approach
        # The return value of `Map()` is the actual value 
        # that will be returned in case there is no condition 
        # (e.g. warning or error). 
        # You don't need to state the return value via `return()` as code 
        # in the "try" part is not wrapped insided a function (unlike that
        # for the condition handlers for warnings and error below)
},
error=function(cond) {
        message(paste("URL does not seem to exist:", URL))
        message("Original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
},
warning=function(cond) {
        message(paste("URL caused a warning:", URL))
        message("Original Warning message:")
        message(cond)
        # Choose a return value in case of warning
        return(NULL)
},
finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you 
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>' 
        message(paste("Processed URL:", URL))
        message("Some other message at the end")
}
        )    
return(out)
}