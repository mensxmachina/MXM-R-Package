tc.plot = function(target, tp, id, type = "l", ylab = "Values", xlab = "Time points",
           col = 2, lwd = 1, lty = 2, pch = 1, main = "Spaghetti plot") {
  
  x <- as.numeric( unique(id) )
  tp <- unique(tp)
  n <- length(x)
  target <- as.vector(target)

  if ( type == "l" ) {
    plot(tp, target[ x == 1 ], type= 'l', col = col, lwd = lwd, lty = lty, 
    xaxt = "n", ylim = c( min(target), max(target) ), ylab = ylab, xlab = xlab, main = main )
    for (i in 2:n) {
      lines(tp, target[x == i], col = col, lwd = lwd, lty = lty)  
    }

  } else {
    plot(tp, target[ x == 1 ], type= 'b', col = col, lwd = lwd, lty = lty, 
    xaxt = "n", ylim = c( min(target), max(target) ), ylab = ylab, xlab = xlab, main = main )
    for (i in 2:n) {
      points(tp, target[x == i])  
      lines(tp, target[x == i], col = col, lwd = lwd, lty = lty)  
    
    }

  }

  axis(1, at = tp, labels = tp, col.axis = "red")
 
}