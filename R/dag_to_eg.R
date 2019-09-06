  ## Converts a DAG into Essential Graph 
  ## Input: dagx, format (0: no edge, 1: directed edge)
  ## Output: eg, format (0: no edge, 1: directed edge, 2:
  ## undirected edge)
  ## 
  ## Base Code is from BNT

dag_to_eg <- function(dagx) {      # function [eg, order] = dag_to_eg(dagx)

  ## Converts a DAG into an Essential Graph
  ## where edges are coded by 2 and 3, 2 is
  ## directed edge and 3 is bidirected edge and is at one (the same as the original DAG) of the two
  ## symetrical places. 

  ## Is implemented by the algorithm of Max Chickering in D.M.Chickering (1995). 
  ## A transformational characterization of equivalent Bayesian network structures. 
  ## In Proceedings of Eleventh Conference on Uncertainty in Artificial Intelligence, Montreal, QU,
  ## pages 87-98. Morgan Kaufmann 
  ## http://research.microsoft.com/~dmax/publications/uai95.pdf 

  ## Implemented by Tomas Kocka, AAU.

  ##print_dag(dagx); ## Just checking input

  ord = topological_sort(dagx); ## get the topological order of nodes and their number

  nx = ny = ncol(dagx); ## gets the number of nodes, note that nx == ny

  if ( any( is.na(ord) ) ) {
    cat('Nodes are not completely ordered.\n');
    eg = matrix(0, nx, ny)
    ord = NaN
  }

  ## fprintf('the topological order is: ##d',order);
  ## fprintf('\n');

  eg = matrix(0, nx, ny)
  a = which(dagx > 0, arr.ind = TRUE)  ## finds all nonzero elements in the adjacency matrix, i.e. arcs in the DAG - however we will overwrite it in a special order
  I = a[, 1]    ;    J = a[,2]
  ## we will sort the arcs from lowest possible y and highest possible x, arcs are x->y
  e = 1;
  for ( y in 1:ny ) {
    for ( x in nx:1 ) {

        ##fprintf('x ##d ',order(x)); fprintf('y ##d ',order(y));

        if ( dagx[ ord[x], ord[y] ] == 1 ) {
            I[e] = ord[x]
            J[e] = ord[y]
            e = e + 1;
            ##fprintf('x order ##d',x);
            ##fprintf('y order ##d',y);
            ##fprintf('\n');
        }
    }
  }


  ## fprintf('the arcs are: ##d',I);
  ## fprintf('\n');
  ## fprintf('the arcs are: ##d',J);
  ## fprintf('\n');


  ## Now we have to decide which arcs are part of the essential graph and
  ## which are undirected edges in the essential graph.
  ## Undecided arc in the DAG are 1, directed in EG are 2 and undirected in EG
  ## are 3.


  for ( e in 1:length(I) ) {
    if ( dagx[ I[e], J[e] ] == 1 ) {
        cont = 1
        for ( w in 1:nx ) {

            if ( dagx[ w, I[e] ] == 2 ) {
                if ( dagx[w, J[e] ] != 0 ) {
                    dagx[ w, J[e] ] = 2

                } else {

                    for ( ww  in 1:nx ) {

                        if ( dagx[ ww, J[e] ] != 0 ) {
                           dagx[ ww, J[e] ] = 2
                        }

                    } ## and now skip the rest and start with another arc from the list

                    w = nx
                    cont = 0
                }
            }
        }

        if ( cont == 1 ) {
 
           exists = 0

           for ( z in 1:nx ) {
               ##fprintf('test ##d',dagx(z,J(e)));

               if ( dagx[ z, J[e] ] != 0  &  z != I[e]  &  dagx[ z, I[e] ] == 0 ) {
                  exists = 1; 
                  for ( ww in 1:nx ) {
                        if ( dagx[ ww, J[e] ] == 1 ) {
                           dagx[ ww, J[e] ]  = 2;
                        } 
                  }

               }
           }

           if ( exists == 0 ) {  
               for ( ww in 1:nx ) {
                   if ( dagx[ ww, J[e] ] == 1 ) {
                      dagx[ ww, J[e] ] = 3
                   } 
               }  
           }

        }
    }
            
  }

  ##print_dag(dagx); ## Just checking output
  ## dagx is now in BNT's eg format, convert to use 
  ## (0: no edge, 1: directed edge, 2: undirected edge format)
  ## Also, make sure that any undirected edge is listed in both
  ## places x1 -- x2 -> (x1, x2) & (x2, x1) = 2;
  eg[ which( dagx == 2 ) ] = 1
  eg[ which( dagx == 3 ) ] = 2

  for ( i in 1:nx ) {
    for ( j in 1:ny ) {
        if ( eg[i, j] == 2 )   eg[j, i] = 2
    }
 }

 list(eg = eg, ord = ord)
}



