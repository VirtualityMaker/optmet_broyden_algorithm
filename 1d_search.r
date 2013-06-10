svenn_algorithm <- function(func, start_point, step)
{
    start_point = c(start_point, func(start_point))
  
    cat("Svenn's Algorithm report:\n")
    cat("Start point: ")
    print(start_point)
    cat("Start step: ")
    print(step)
    cat("\n")
                    
    points <- list(start_point)
  
    a = c(start_point[1] - step, func(start_point[1] - step))
    b = c(start_point[1] + step, func(start_point[1] + step))
    
    if (a[2] >= start_point[2] && b[2] >= start_point[2]) {
      cat("a: ")
      print(a)
      cat("\n")
      cat("b: ")
      print(b)
      cat("\n")
      return(list(start_point, start_point))
    }
    
    step_mod = 2

    #     |
    # . | |
    if ((a[2] <= b[2]) &&
        (a[2] < points[[1]][2])) {
        step_mod <- step_mod * (-step)
        points[[length(points)+1]] = points[[length(points)]]
        points[[length(points)-1]] <- b
        points[[length(points)+1]] <- a
    }

    # |   
    # | | .
    else if ((b[2] <= a[2]) &&
             (b[2] < points[[1]][2])) {
        step_mod <- step_mod * step
        points[[length(points)+1]] = points[[length(points)]]
        points[[length(points)-1]] <- a
        points[[length(points)+1]] <- b
    }
  
    # |   |
    # | | |
    else if ((a[2] == b[2]) &&
             (a[2]  > points[[1]][2])) {
        cat("RESULT:\n\n")
        print(list(a, b))
        
        return(list(a, b))
    }
  
    while (TRUE)
    {
        prev_point = points[[length(points)]]
        next_point = c(prev_point[1] + step_mod, func(prev_point[1] + step_mod))
        
        if (next_point[2] < prev_point[2])
        {
            points[[length(points)+1]] <- next_point
            step_mod <- step_mod * 2
        }
        else if (next_point[2] >= prev_point[2])
        {
            step_mod <- step_mod / 2
            points[[length(points)+1]] <- c(next_point[1] - step_mod, func(next_point[1] - step_mod))
           
            print(points)
            cat("RESULT:\n\n")
            
            result = list()
            i = length(points)
            if (i == 2) {
                result = points
            } else {
                if (step_mod > 0) {
                    result = list(points[[i - 2]], points[[i]])
                } else {
                    result = list(points[[i]], points[[i - 2]])
                }
            }
            
            print(result)
            
            return(result)
        }
    }
}

method_golden_ratio <- function(func, a, b, accuracy) {
    L  = b[1] - a[1]
    
    x1 = c(a[1] + 0.382 * L, func(a[1] + 0.382 * L))
    x2 = c(a[1] + 0.618 * L, func(a[1] + 0.618 * L))
    
    #cat("Golden Ratio Method report:\n\n")
    #cat('a:  ')
    #print(a)
    #cat('x1: ')
    #print(x1)
    #cat('x2: ')
    #print(x2)
    #cat('b:  ')
    #print(b)

    while (L > accuracy) {
        if (x1[2] >= x2[2]) {
            a  = x1
            L  = b[1] - a[1]
            
            x1 = x2
            x2 = c(a[1] + 0.618 * L, func(a[1] + 0.618 * L))
        }
        else if (x1[2] <= x2[2]) {
            b  = x2
            L  = b[1] - a[1]
            
            x2 = x1
            x1 = c(a[1] + 0.382 * L, func(a[1] + 0.382 * L))
        }
        
        #cat('\n')
        #cat('a:  ')
        #print(a)
        #cat('x1: ')
        #print(x1)
        #cat('x2: ')
        #print(x2)
        #cat('b:  ')
        #print(b)
    }
    
    #cat("------// END OF REPORT //------\n")
    
    return(list(a,b))
}

method_dihotomy <- function(func, a, b, accuracy)
{
    cat("Method dihotomy\n")
    cat("--------------------------------\n")

    L = b - a
    
    xmiddle  = (a + b) / 2
    fxmiddle = func(xmiddle)
    
    while (L > accuracy) {
        x1 = a + L/4
        x2 = b - L/4
        
        fx1 = func(x1)
        fx2 = func(x2)
        
        cat("a          = ", a, "\n")
        cat("b          = ", b, "\n")
        cat("L          = ", L, "\n")
        cat("x1         = ", x1, "\n")
        cat("f(x1)      = ", fx1, "\n")
        cat("x2         = ", x2, "\n")
        cat("f(x2)      = ", fx2, "\n")
        cat("xmiddle    = ", xmiddle, "\n")
        cat("f(xmiddle) = ", fxmiddle, "\n\n")
        
        if (fx1 <= fxmiddle) {
            b        = xmiddle
            xmiddle  = x1
            fxmiddle = fx1
        } else {
            if (fx2 < fxmiddle) {
                a        = xmiddle
                xmiddle  = x2
                fxmiddle = fx2
            } else if (fx1 > fxmiddle) {
                a = x1
                b = x2
            }
        }

        L = L / 2
    }
    
    return(c(xmiddle, fxmiddle))
}

method_dsk <- function(func, start_point, svenn_step = 0.1, accuracy_var = 0.1, accuracy_func = accuracy_var)
{
    start_point = c(start_point, func(start_point))

    while (TRUE) {
        lst_bounds = svenn_algorithm(func, start_point[1], svenn_step)
        a = lst_bounds[[1]]
        b = lst_bounds[[2]]
        
        svenn_step = svenn_step/10
    
        central = c((a[1] + b[1]) / 2, func((a[1] + b[1]) / 2))
        deltax  = abs(b[1] - central[1])
        
        if (a[2] - 2*central[2] + b[2] == 0) return(a)
        
        xpoint = central[1] + deltax * (a[2] - b[2]) / (2 * (a[2] - 2*central[2] + b[2]))
        x = c(xpoint, func(xpoint))
        
        # find minimal point
        xmin = a
        for (itm in list(b, central)) {
          if (xmin[2] > itm[2]) {
            xmin = itm
          }
        } 

        cat("A:                      ")
        print(a)
        cat("Central point:          ")
        print(central)
        cat("B:                      ")
        print(b)
        cat("Point of minimum:       ")
        print(xmin)
        cat("\nResult point:           ")
        print(x)
        cat("Variables check error:  ")
        print(abs(xmin[1] - x[1]))
        cat("Functions check error:  ")
        print(abs(xmin[2] - x[2]))
        cat('\n')
        
        if ((abs(xmin[1] - x[1]) <= accuracy_var) &&
            (abs(xmin[2] - x[2]) <= accuracy_func)) {
            return(x)
        } else {
            if (xmin[2] < x[2]) {
              start_point = xmin  
            } else {
              start_point = x
            }            
        }
    }
}

method_dsk_pauell <- function(func, start_point, svenn_step = 0.1, accuracy_var = 0.1, accuracy_func = accuracy_var)
{
  start_point = c(start_point, func(start_point))
  
  lst_bounds = svenn_algorithm(func, start_point[1], svenn_step)
  a = lst_bounds[[1]]
  b = lst_bounds[[2]]
  
  central = c((a[1] + b[1]) / 2, func((a[1] + b[1]) / 2))
  deltax  = abs(b[1] - central[1])
  
  if (a[2] - 2*central[2] + b[2] == 0) return(a)
  
  xpoint = central[1] + deltax * (a[2] - b[2]) / (2 * (a[2] - 2*central[2] + b[2]))
  x = c(xpoint, func(xpoint))
  
  # find minimal point
  xmin = a
  for (itm in list(b, central)) {
    if (xmin[2] > itm[2]) {
      xmin = itm
    }
  } 
  
  cat("A:                      ")
  print(a)
  cat("Central point:          ")
  print(central)
  cat("B:                      ")
  print(b)
  cat("Point of minimum:       ")
  print(xmin)
  cat("\nResult point:           ")
  print(x)
  cat("Variables check error:  ")
  print(abs(xmin[1] - x[1]))
  cat("Functions check error:  ")
  print(abs(xmin[2] - x[2]))
  cat('\n')
  
  if ((abs(xmin[1] - x[1]) <= accuracy_var) &&
        (abs(xmin[2] - x[2]) <= accuracy_func)) {
    return(x)
  } else {
    if (xmin[2] < x[2]) {
      return(xmin)
    }    
  }
  
  x_prev = x
  
  cat("\nStart of Pauell's algorithm\n\n")
  while (TRUE) { 
    if (x[1] < a[1]) {
      b = central
      central = a
      a = x
    } else if (x[1] < central[1]) {
      b = central
      central = x
    } else if (x[1] < b[1]) {
      a = central
      central = x
    } else if (x[1] > b[1]) {
      a = central
      central = b
      b = x
    }
    
    if (central[1] - a[1] == 0) return(a)
    if (b[1] - central[1] == 0) return(b)
    if (b[1] - a[1] == 0) return(b)
    
    a1 = (central[2] - a[2])/(central[1] - a[1])
    a2 = ((b[2] - a[2]) / (b[1] - a[1]) - a1) / (b[1] - central[1])
    
    if (a2 == 0) return(b)
    
    x  = ((a[1] + central[1]) / 2) - (a1 / (2*a2))
    
    x  = c(x, func(x))
    
    if (x[2] == x_prev[2] && x[1] == x_prev[1]) return(x)
    x_prev = x
    
    xmin = a
    for (itm in list(b, central)) {
      if (xmin[2] > itm[2]) {
        xmin = itm
      }
    }
    
    cat("A:                      ")
    print(a)
    cat("Central point:          ")
    print(central)
    cat("B:                      ")
    print(b)
    cat("Point of minimum:       ")
    print(xmin)
    cat("\nResult point:           ")
    print(x)
    cat("Variables check error:  ")
    print(abs(xmin[1] - x[1]))
    cat("Functions check error:  ")
    print(abs(xmin[2] - x[2]))
    cat('\n')
    
    if ((abs(xmin[1] - x[1]) <= accuracy_var) &&
          (abs(xmin[2] - x[2]) <= accuracy_func)) {
      return(x)
    } else {
      if (x[2] > xmin[2]) {
        t = x
        x = xmin
        xmin = x
        
        if (xmin[1] < a[1]) {
          b = central
          central = a
          a = xmin
        } else if (xmin[1] < central[1]) {
          b = central
          central = xmin
        } else if (xmin[1] < b[1]) {
          a = central
          central = xmin
        } else if (xmin[1] > b[1]) {
          a = central
          central = b
          b = xmin
        }
      }    
    }
  }
}