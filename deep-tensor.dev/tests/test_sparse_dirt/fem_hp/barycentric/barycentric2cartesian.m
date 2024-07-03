function x = barycentric2cartesian(vertices, lambdas)
%barycentric2cartesian
%
% transform from barycentric to cartesian
%
% Tiangang Cui, 29/Oct/2016
%
% inputs: 
%   lambdas:    n x nd array
%   vertices:   nd+1 x nd
%
% coordinates of all the inputs and outputs points are aligned columnwise
% along the row is the indices of points

x = [lambdas, 1-sum(lambdas,2)]*vertices;

end
