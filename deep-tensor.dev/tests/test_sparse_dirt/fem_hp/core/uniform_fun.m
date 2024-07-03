function f = uniform_fun(x, left, right)

tol = 1E-10;
ind = x < left-tol | x > right+tol;
f   = ones(size(x));
f(ind) = 0;

end