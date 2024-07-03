function M = align_cartesian_pts(pts)
%align_cartesian_pts
%
% adding a one column vector to the catesian points, and then the output
% can be used in coordinate transformation 
%
% Tiangang Cui, 29/Oct/2016
%
% inputs: 
%   pts: n x d array
%

M = [pts, ones(size(pts,1),1)];

end