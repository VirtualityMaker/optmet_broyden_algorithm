real_grad_f_calculation_forward_difference <- function(x, h)
{
  grad_f_1 = (f(x + h * c(1, 0)) - f(x)) / h
  grad_f_2 = (f(x + h * c(0, 1)) - f(x)) / h
  
  return(c(grad_f_1, grad_f_2))
}

real_grad_f_calculation_central_difference <- function(x, h)
{
  grad_f_1 = (f(x + h * c(1, 0)) - f(x - h * c(1, 0))) / (2*h)
  grad_f_2 = (f(x + h * c(0, 1)) - f(x - h * c(0, 1))) / (2*h)
  
  return(c(grad_f_1, grad_f_2))
}

real_grad_f_calculation_backward_difference <- function(x, h)
{
  grad_f_1 = (f(x) - f(x - h * c(1, 0))) / h
  grad_f_2 = (f(x) - f(x - h * c(0, 1))) / h
  
  return(c(grad_f_1, grad_f_2))
}

real_f_calculation <- function(x)
{
  100*(x[1]^2- x[2])^2 + (x[1] - 1)^2
  #(10*(x[1] - x[2])^2 + (x[1] - 1)^2)^4
  #(10*(x[1] - x[2])^2 + (x[1] - 1)^2)^(1/4)
  #(x[1]^2 + x[2]^2)
}

#real_grad_f_calculation_analytic <- function(x)
#{
#  return(c(
#    400*x[1]*(x[1]^2 - x[2]) + 2*(x[1]-1),
#    -200*(x[1]^2 - x[2])))
#}

check_for_vector <- function(x)
{
  return(function(a){ a[1] == x[1] && a[2] == x[2] })
}

f <- function(x)
{
  f_values <- get('f_values', envir=.GlobalEnv)
  f_values[[length(f_values)+1]] = x
  assign("f_values", f_values, envir=.GlobalEnv)

  value = real_f_calculation(x)
  
  return(value)
}

grad_f <- function(x)
{
  f_grad_values <- get('f_grad_values', envir=.GlobalEnv)
  f_grad_values[[length(f_grad_values)+1]] = x
  assign("f_grad_values", f_grad_values, envir=.GlobalEnv)
  
  #cat("Grad analytic calculation: ", real_grad_f_calculation_analytic(x), "\n")
  #cat("Grad difference scheme calculation: ", real_grad_f_calculation_backward_difference(x, 0.000001), "\n")
  #value = real_grad_f_calculation_central_difference(x, 10^(-12))
  value = real_grad_f_calculation_central_difference(x, 10^(-12))
  #value = real_grad_f_calculation_analytic(x)
  
  return(value)
}

f_values <- list()
f_grad_values <- list()