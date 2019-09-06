nchoosek = function(cs , k) {  #i can also pass the compFun arg for selecting
  if (length(cs) == 1) { #if not vector
    res = choose(cs, k); #or nchoosek
  } else  res = combn(cs, k)
  res;
}