function y = Fourier_sqrconv(coef, poly1, poly2, FG)
%{
coef is a cell array, containing the two vectors of coeffcients
poly1 is the Fourier bases for the input functions
poly2 is the Fourier bases for their square
FG is the modifier for Gaussian convolution
out put y contains the coeffcient of the product convoluted with Gaussian 
%}

c1 = coef{1};
c2 = coef{2};

fun = @(x) eval_radon(poly1, c1, x) .*eval_radon(poly1, c2, x);

nodes = poly2.nodes;
y = poly2.node2basis * fun(nodes);

y = y.*FG;

%
end