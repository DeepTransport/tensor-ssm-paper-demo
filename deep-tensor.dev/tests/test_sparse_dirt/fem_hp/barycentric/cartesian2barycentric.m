function lambdas = cartesian2barycentric(vertices_invT, pts)
%cartesian2barycentric
%
% transform from cartesian barycentric
%
% Tiangang Cui, 29/Oct/2016
%
% inputs: 
%   pts:            n x d array
%   vertices_invT:  nd+1 x nd+1, precomputed for inverse transformation 
%
% coordinates of all the inputs and outputs points are aligned columnwise
% along the row is the indices of points

tmp = align_cartesian_pts(pts)*vertices_invT;
lambdas = tmp(:,1:end-1);

if mean(abs(sum(tmp,2))) < 1E-12
    disp('Error in cartesian 2 barycentric')
end

end
