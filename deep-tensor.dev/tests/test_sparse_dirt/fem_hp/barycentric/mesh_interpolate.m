function map = mesh_interpolate(mesh, local_elem, pts, h, debug)
%mesh_interpolate
%
% evaluating interpolation matrix for a given grid and a given point set
%
% Tiangang Cui, 8/Dec/2016

patch_size = 2;

n   = size(pts, 1);
tol = h/20;

% find element for each pts
lambdas = zeros(n, mesh.num_dim);
elem_indices = zeros(n, 1);

for i = 1:n
    ind = true(mesh.num_elem, 1);
    for d = 1:(mesh.num_dim+1)
        ind = ind & ( sqrt(sum( (repmat(pts(i,:), mesh.num_elem, 1) - mesh.grid_p(mesh.grid_t(:, d), :)).^2 , 2)) < h*patch_size );
    end
    is = find(ind);
    ls = zeros(length(is), mesh.num_dim);
    for j = 1:length(is)
        %pts = mesh.grid_p(mesh.grid_t(is(j), :), :);
        cartesian2barycentric(mesh.vertices_invT{is(j)}, pts(i,:));
        ls(j,:) = cartesian2barycentric(mesh.vertices_invT{is(j)}, pts(i,:));
    end
    jnd = sum(abs(ls), 2) < 1 & sum(ls > 0, 2) == mesh.num_dim;
    if sum(jnd) == 1 % interior 
        elem_indices(i) = is(jnd);
        lambdas(i,:) = ls( jnd, : );
    elseif sum(jnd) > 1
        js  = find(jnd);
        [~, col] = min( abs(sum(ls(js,:), 2) - 1) );
        elem_indices(i) = is( js(col) );
        lambdas(i,:)  = ls( js(col), :);
    else
        knd = sum(abs(ls), 2) < (1+tol) & sum(ls > -tol, 2) == mesh.num_dim;
        ks  = find(knd);
        [~, col] = min( abs(sum(ls(ks,:), 2) - 1) );
        elem_indices(i) = is( ks(col) );
        lambdas(i,:)  = ls( ks(col), :);
    end
end

fs = zeros(n, local_elem.num_node);
ii = zeros(n, local_elem.num_node);
jj = zeros(n, local_elem.num_node);
% interpolate 
for i = 1:n
    fs(i,:) = eval_barycentric_bases(local_elem.alpha, local_elem.A, lambdas(i,:));
    ii(i,:) = i;
    jj(i,:) = mesh.node_map(elem_indices(i),:);
end

map = sparse(ii, jj, fs, n, mesh.dof);

if debug 
    figure;
    m = ceil(sqrt(n));
    for i = 1:n
        subplot(m, m, i)
        ind = mesh.grid_t(elem_indices(i), :);
        plot(mesh.grid_p(ind,1), mesh.grid_p(ind,2), 'bo')
        hold on
        plot(pts(i,1), pts(i,2), 'rx')
    end
    
end

end