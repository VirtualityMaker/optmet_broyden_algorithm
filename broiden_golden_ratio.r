source("1d_search.r");
source("data.r");

f_by_direction <- function(f, start_point, direction)
{
    return(function(lambda){ f(start_point + lambda * direction); });
}

get_next_hesse_approximation <- function(prev_A, delta_x, delta_g)
{
    delta_x = matrix(delta_x, nrow=2, ncol=1)
    delta_g = matrix(delta_g, nrow=2, ncol=1)
    
    #~ broidena-shafe
    
    #mult1 = diag(2) - (delta_x %*% t(delta_g)) / (t(delta_x) %*% delta_g)[1,1]
    #mult2 = (delta_x %*% t(delta_x)) / (t(delta_x) %*% delta_g)[1,1]

    #A = mult1 %*% prev_A %*% mult1 + mult2
    
    #~ broidena
    
    mult1 = delta_x - prev_A %*% delta_g
  
    A = prev_A + (mult1 %*% t(mult1)) / (t(mult1) %*% delta_g)[1,1]
    
    #~ devisona-fletchera-pauella
    
    #mult1 = (delta_x %*% t(delta_x)) / (t(delta_x) %*% delta_g)[1,1]
    #mult2 = (prev_A %*% delta_g %*% t(delta_g) %*% prev_A) / (t(delta_g) %*% prev_A %*% delta_g)[1,1]
    
    #A = prev_A + mult1 - mult2
    
    cat("\nApproximated matrix:\n")
    print(A)
    cat("\n")
    
    return(A)
}

get_svenn_accuracy <- function(i)
{
  0.6 ^ (i)
}

get_1d_accuracy <- function(i)
{
  10^(-2)
}

x_prev = c(-1.2, 0)
A0 = matrix(c(1,0,0,1),nrow=2,ncol=2)
#x_prev = c(-1, -1)
A_prev = A0

cat("Start point: \n")
print(x_prev)
cat("Start Hesse approximation: \n")
print(A_prev)
cat("\n")

x_points = c()
y_points = c()
z_points = c()

points = list()
points[[1]] = c(x_prev, f(x_prev))
x_points = c(x_points, x_prev[1])
y_points = c(y_points, x_prev[2])
z_points = c(z_points, f(x_prev))

restarts = list()

maximum_iterations_number = 100
epsilon_for_finish_criteria = 10^(-6)

iteration_number = 0
adjacent_restarts_number = 0
grinshtadts_recalculations = 0

abs_of_deltag = Inf
norm_of_deltax = Inf

for (i in 1:maximum_iterations_number) {
  if (adjacent_restarts_number == 2) break
  
  iteration_number = iteration_number + 1
  
  cat("--------- Iteration #", i, " --------\n\n")
  #eigenvals = eigen(A_prev, TRUE, TRUE)
  
  #for (j in 1:length(eigenvals$values)) {
  #  if (eigenvals$values[j] < 0) {
  #    cat("Not positive definite matrix!\n")
  #    grinshtadts_recalculations = grinshtadts_recalculations + 1
      
  #    if (abs(A_prev[1,1]) < 0.0001) A_prev[1,1] = 0.0001;
  #    if (abs(A_prev[2,2]) < 0.0001) A_prev[2,2] = 0.0001;
      
  #    C_ = matrix(c((1/abs(A_prev[1,1])^0.5), 0, 0, (abs(1/A_prev[2,2])^0.5)), nrow=2, ncol=2)
  #    P = C_ %*% A_prev %*% C_
      
  #    ev = eigen(P, TRUE)
      
  #    for (t in 1:length(ev$values)) {
  #      if (ev$values[t] < 10^(-4)) {
  #        ev$values[t] = 10^(-4) 
  #      }
  #    }
      
  #    A_prev = abs(ev$values[1]) * matrix(c(ev$vectors[1,1],ev$vectors[1,2],0,0), nrow=2, ncol=2) %*% 
  #      matrix(c(ev$vectors[1,1],0,ev$vectors[1,2],0), nrow=2, ncol=2) +
  #      abs(ev$values[2]) * matrix(c(ev$vectors[2,1],ev$vectors[2,2],0,0), nrow=2, ncol=2) %*% 
  #      matrix(c(ev$vectors[2,1],0,ev$vectors[2,2],0), nrow=2, ncol=2)
      
  #    print(A_prev)
  #    cat("\n")
  #    break
  #  }
  #}

  s = - A_prev %*% grad_f(x_prev);
  
  if (norm(s) == 0) {
    A_prev = A0
    restarts[[length(restarts)+1]] = c("norm_s", iteration_number)
    adjacent_restarts_number = adjacent_restarts_number + 1
    
    cat("RESTART!\n\n")
    next
  }
  
  s = s / norm(s, "f") 
  
  cat("Direction:        \n")
  print(s)
  cat("\n")
  
  #lambda=0.001
  
  #lambda = method_dsk(f_by_direction(f, x_prev, s), 0, get_svenn_accuracy(iteration_number), get_1d_accuracy(iteration_number))[1];

  svenn_step =  min(get_svenn_accuracy(i), 0.1 * norm(matrix(x_prev, nrow=1, ncol=2), "f"))
  #svenn_step = get_svenn_accuracy(i)
  
  svenn_bounds = svenn_algorithm(f_by_direction(f, x_prev, s), 0, svenn_step)
  lambda = method_golden_ratio(f_by_direction(f, x_prev, s), svenn_bounds[[1]], svenn_bounds[[2]], get_1d_accuracy(i))
  lambda = lambda[[1]][1]
  
  #svenn_bounds = svenn_algorithm(f_by_direction(f, x_prev, s), 0, svenn_step)
  #lambda = method_dihotomy(f_by_direction(f, x_prev, s), svenn_bounds[[1]][1], svenn_bounds[[2]][1], get_1d_accuracy(i))
  #lambda = lambda[1]
  
  #lambda = method_dsk_pauell(f_by_direction(f, x_prev, s), 0, get_svenn_accuracy(iteration_number), get_1d_accuracy(iteration_number))[1]
  
  cat("Lamda:            ", lambda, "\n")
  
  if (lambda <= 0) {
    A_prev = A0
    restarts[[length(restarts)+1]] = c(lambda[1], iteration_number)
    adjacent_restarts_number = adjacent_restarts_number + 1
    
    cat("RESTART!\n\n")
    next
  }
  
  adjacent_restarts_number = 0
  
  x_cur = x_prev + lambda * s
  x_cur = c(x_cur[1,1], x_cur[2,1])
  
  cat("\nx[", i,"]: ")
  print(x_cur)
  cat("\n")
  
  points[[length(points)+1]] = c(x_cur, f(x_cur))
  x_points = c(x_points, x_cur[1])
  y_points = c(y_points, x_cur[2])
  z_points = c(z_points, f(x_cur))

  #norm_of_deltax = norm(matrix(x_prev - x_cur, nrow=1, ncol=2), "f") / norm(matrix(x_cur, nrow=1, ncol=2), "f")
  #abs_of_deltag  = abs(f(x_prev) - f(x_cur))
  
  #cat("|f(x_k+1) - f(x_k)| = ")
  #print(abs_of_deltag)
  #cat("||x_k+1 - x_k||     = ")
  #print(norm_of_deltax)
  
  #if (norm_of_deltax < epsilon_for_finish_criteria &&
  #   abs_of_deltag < epsilon_for_finish_criteria) break
  
  cat("||grad_f(x1)||     = ")
  print(norm(matrix(grad_f(x_cur), nrow=1, ncol=2), "f"))
  
  if (norm(matrix(grad_f(x_cur), nrow=1, ncol=2), "f") < epsilon_for_finish_criteria) break
  
  A_cur = get_next_hesse_approximation(A_prev, x_cur - x_prev, grad_f(x_cur) - grad_f(x_prev))

  x_prev = x_cur
  A_prev = A_cur
}

cat("Number of iterations:  ", iteration_number, "\n")
cat("Grinshtadts recalculations:  ", grinshtadts_recalculations, "\n")
cat("Function calculations: ", length(unique(f_values)), "\n")
cat("Gradient calculations: ", length(unique(f_grad_values)), "\n")
cat("Points:\n")
print(points)
cat("Restarts:\n")
print(restarts)
cat("Last A:\n")
print(A_cur)

x <- seq(-2, 2, length.out = 20)  
y <- x

rotf <- Vectorize(function(x,y){f(c(x,y))})

z <- outer(x,y,rotf)
contour(x,y,z)
points(x_points, y_points, "b")
