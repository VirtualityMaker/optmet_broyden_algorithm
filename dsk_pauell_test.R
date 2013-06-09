source ("1d_search.r")

f <- function(x)
{
  10*(x^2 - 1)^2
}

print(method_dsk_pauell(f, 2, 0.1, 0.001))