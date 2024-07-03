function true_x = load_test_image(mesh, prior, func, type, base, range)
%INIT_TEST_IMAGE
%
% initializes the test image for distributed parameters
%
% Tiangang Cui, 11/May/2014

switch type
    case {'CF'}
        true_x = conduct_cf(mesh, base, range);
    case {'Prior'}
        true_x = conduct_prior(prior, func);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cond = conduct_prior(prior, func)

r     = randn(prior.dof, 1);
u     = matvec_prior_L(prior, r) + prior.mean_u;
cond  = u2x(func, u);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cond = conduct_cf(mesh, base, range)

ref = [0.75, 0.25];
r = ref(1)-ref(2);
c = 0.06+1e-10;

alphas = linspace(0,pi/2,20);

ind = false(mesh.dof,1);
for i = 1:length(alphas)
    x = ref(1) - r*cos(alphas(i));
    y = ref(2) + r*sin(alphas(i));
    
    j = sqrt( (mesh.nodes(:,1)-x).^2 + (mesh.nodes(:,2)-y).^2 ) < c;
    ind = ind | j(:);
end

c = 0.1+1e-10;
ref = [0.7 0.3];
j = mesh.nodes(:,1)>(ref(1)-c) & mesh.nodes(:,1)<(ref(1)+c) & mesh.nodes(:,2)>(ref(2)-c) & mesh.nodes(:,2)<(ref(2)+c);
ind = ind | j(:);

cond = ones(mesh.dof,1)*(base+range);
cond(ind) = base;

end
