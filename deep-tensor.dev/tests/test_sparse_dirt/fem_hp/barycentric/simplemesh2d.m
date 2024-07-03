function [p, t] = simplemesh2d(h, a, b)
%simplemesh2d 
%
% T Cui, 17 Oct 2016
%
%
% The geometry of the mesh looks like: (for N_side=4)
%
%     *----*----*----*----*
%     | \  | \  | \  | \  |
%     |  \ |  \ |  \ |  \ |
%     *----*----*----*----*
%     | \  | \  | \  | \  |
%     |  \ |  \ |  \ |  \ |
%     *----*----*----*----*
%     | \  | \  | \  | \  |
%     |  \ |  \ |  \ |  \ |
%     *----*----*----*----*
%     | \  | \  | \  | \  |
%     |  \ |  \ |  \ |  \ |
%     *----*----*----*----*
%
% where * are the FE nodes.

m_side = a/h;
n_side = b/h;
x = linspace(0,a,m_side+1);
y = linspace(0,b,n_side+1);

n_elem  = n_side*m_side*2;       % number of elements

[xx,yy] = meshgrid(x,y);  % Generate mesh of node coordinates.

% Now number the nodes (from top left downwards)
p       = [xx(:), yy(:)];

% Generate local-to-global node number mapping
% Local node numbering is:
%
%   4---3
%   |   |
%   1---2
% 1 2 4 and 2 3 4 form two triangle mesh
n_ref   = 3;             % number of nodes on reference element
t       = zeros(n_elem, n_ref);

blnode = 1; %assemble the node map from the bottom left node
for i = 1:n_elem/2
    
    t(i*2-1,:)  = blnode + [0, n_side+1, 1];
    t(i*2,:)    = blnode + [n_side+1, n_side+2, 1];
    blnode = blnode + 1;
    if rem(blnode,n_side+1) == 0 % reach the top of the column (y)
        blnode = blnode+1; % skip to next column
    end
end

end