function [local_elem, local_elem_bnd] = local_elem_simplex(num_dim, poly_order, quad_order)
%local_elem_simplex
%
% build local element in a num_dim dimensional space, with polynomial order
% poly_order and quadrature order quad_order
%
% Tiangang Cui, 29/Oct/2016

% reference triangle (1,0) (0,1) (0,0)
% coordinates: 
%   x:      physical space
%   lambda: barycentric
%
% bases.f_quad_pts, bases.grad_f_quad_pts:
%   evaluation of basis functions and their gradient on quadrature points 
%   num_quad x num_node

tol = 1e-10;

% cartesian interpolation points by (Matteo Briani, Alvise Sommariva and Marco Vianello)
if num_dim == 2
    if (poly_order < 1 || poly_order > 18)
        disp('Polynomial order not supported');
    end
    x_pts = set_simplex_leb_gll(poly_order);
end

% convert to barycentric
vertices = [eye(num_dim); zeros(1,num_dim)];
vertices_invT = eye(num_dim+1)/align_cartesian_pts(vertices);
lambda_pts = cartesian2barycentric(vertices_invT, x_pts);

% assemble local element
local_elem = assemble_barycenteric_bases(lambda_pts, poly_order, quad_order);
local_elem.vertices = vertices;
local_elem.x_pts = x_pts;
local_elem.lambda_pts = lambda_pts;
local_elem.poly_order = poly_order; 
local_elem.quad_order = quad_order;

% bounadry points
bnd_lambda_pts = x_pts(abs(x_pts(:,end)) < tol, 1:end-1);

% evaluate boundary terms
local_elem_bnd = assemble_barycenteric_bases(bnd_lambda_pts, poly_order, quad_order);
local_elem_bnd.lambda_pts = bnd_lambda_pts;
local_elem_bnd.poly_order = poly_order; 
local_elem_bnd.quad_order = quad_order;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bases = assemble_barycenteric_bases(lambda_pts, poly_order, quad_order)

num_dim = size(lambda_pts, 2);

% set multi index
[bases.alpha, ~]   = multi_index(num_dim, poly_order);
if size(lambda_pts, 1) ~= size(bases.alpha, 1)
    disp('Error: Number of points does not equal to dimension of polynomials')
    return 
end

% evaluate polynomial constants in barycentric space 
bases.A = barycentric_bases(bases.alpha, lambda_pts);

% define gradient of barycentric bases
[bases.grad_alpha, bases.grad_A] = grad_barycentric_bases(bases.alpha, bases.A);

% define quadrature points (in barycentric) and weights
switch num_dim
    case 1
        [bases.quad_lambdas, bases.quad_weights] = gauss_rule(quad_order, 0, 1);
    case 2
        [bases.quad_lambdas, bases.quad_weights] = dunavant_rule(quad_order, 1);
        bases.quad_lambdas = bases.quad_lambdas';
        bases.quad_weights = bases.quad_weights';
        bases.quad_weights = 0.5*bases.quad_weights;
    case 3
        disp('not implemented')
end

bases.num_dim = num_dim;
bases.num_quad = size(bases.quad_lambdas, 1);
bases.num_node = size(lambda_pts, 1);

% evaluate barycentric bases and its gradient at quadrature points
bases.f_quad_pts = eval_barycentric_bases(bases.alpha, bases.A, bases.quad_lambdas);
bases.grad_f_quad_pts = cell(num_dim,1);
for d = 1:num_dim
    bases.grad_f_quad_pts{d} = eval_barycentric_bases(...
        bases.grad_alpha{d}, bases.grad_A{d}, bases.quad_lambdas);
end

% evaluate gradient of barycentric bases functions at nodal points
% bases.grad_f_lambda_pts = cell(num_dim,1);
% for d = 1:num_dim
%     bases.grad_f_lambda_pts{d} = eval_barycentric_bases(...
%         bases.grad_alpha{d}, bases.grad_A{d}, lambda_pts);
% end

% local mass and stiffness matrix, without Jacobian scaling 
bases.mass = bases.f_quad_pts'*diag(bases.quad_weights)*bases.f_quad_pts;
stiff_partial = cell(num_dim);
for di = 1:num_dim
    for dj = 1:num_dim
        stiff_partial{di,dj} = ...
            bases.grad_f_quad_pts{di}'*diag(bases.quad_weights)*bases.grad_f_quad_pts{dj};
    end
end


% local mass and stiffness matrix at quadrature points
bases.mass_at_q = zeros(bases.num_node^2, bases.num_quad);
for q = 1:bases.num_quad
    tmp = bases.f_quad_pts(q,:)'*bases.f_quad_pts(q,:);
    bases.mass_at_q(:,q) = tmp(:);
end
bases.stiff_partial_at_q = cell(num_dim);
for di = 1:num_dim
    for dj = 1:num_dim
        bases.stiff_partial_at_q{di,dj} = zeros(bases.num_node^2, bases.num_quad);
        for q = 1:bases.num_quad
            tmp = bases.grad_f_quad_pts{di}(q,:)'*bases.grad_f_quad_pts{dj}(q,:);
            bases.stiff_partial_at_q{di,dj}(:,q) = tmp(:);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = barycentric_bases(alpha, lambda_pts)

if size(lambda_pts,1) ~= size(alpha,1)
    disp('Error: Number of points does not equal to dimension of polynomials')
end

n = size(lambda_pts,1);

P = zeros(n);
for i = 1:n
    P(:,i) = prod(repmat(lambda_pts(i,:), n, 1).^alpha, 2);
end
disp(['det(P) = ' num2str(det(P))])

A = eye(n)/P;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grad_alpha, grad_A] = grad_barycentric_bases(alpha, A)

[n, nd] = size(alpha);

grad_A = cell(nd, 1);
grad_alpha = cell(nd, 1);

% gradient
for d = 1:nd
    % order of the polynomial 
    order = alpha(:,d);
    % indices where order is 0
    ind = order == 0;
    % constant matrix
    grad_A{d} = A.*repmat(order', n, 1);
    % reduce the order of the polynomial along dim d
    grad_alpha{d} = alpha;
    grad_alpha{d}(:,d) = order-1;
    % set order to be zero
    grad_alpha{d}(ind,d) = 0;
end

end