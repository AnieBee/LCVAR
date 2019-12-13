createOutputList <- function(LagCombinations, Rand, Rational, Initialization, PreviousSol)
{
    OutputListLarge =  vector(mode = "list", LagCombinations)
    OutputListSmall = vector(mode = "list", Rand + as.numeric(Rational) + as.numeric(!is.null(Initialization))
                                            + as.numeric(PreviousSol)) # to store output of all starts
    
    for (smallList in 1:LagCombinations)
    {
        OutputListLarge[[smallList]] = OutputListSmall
    }
    
    invisible(OutputListLarge)

}