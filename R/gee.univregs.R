gee.univregs <- function(target, reps = NULL, id, dataset, targetID = -1, test, wei = NULL, 
                              correl = "echangeable", se = "jack", ncores = 1) {
    univariateScore.gee(target = target, reps = reps, group = id, dataset = dataset, test = test, wei = wei, 
                                              targetID = targetID, correl = correl, se = se, ncores = ncores)
    
}  
