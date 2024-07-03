%%
model = ShockAbsorber();
f = @(theta) potential(model, theta);

%%
% Ranges of domain
a = -3*ones(1,model.D+2);
b = 3*ones(1,model.D+2);

%%
p = 32; % order
r = 16; % rank

% Vectors of grid points2
polys = cell(model.D+2,1);
domains = cell(model.D+2,1);
for i=1:model.D+2
    polys{i} = Lagrange1(p);
    domains{i} = BoundedDomain([a(i),b(i)]);
end
base = ApproxBases(polys, domains);
opt = TTOption('max_als', 4, 'als_tol', 1E-8, 'local_tol', 1E-20, 'kick_rank', 2, 'init_rank', r, 'max_rank', r);

ttirt = TTIRT(f, base, opt);

%test_irt(f, ttirt, 1E4)

%%
debug_size = 1E4;
debug_x = random(ttirt,debug_size);
sample_x = random(ttirt,debug_size);
deb = InputData(sample_x, debug_x);

% Vectors of grid points
polys = cell(model.D+2,1);
domains = cell(model.D+2,1);
for i=1:model.D+2
    polys{i} = Chebyshev1st(p);
    domains{i} = AlgebraicMapping(1);
end
base = ApproxBases(polys, domains);
opt = TTOption('max_als', 4, 'als_tol', 1E-8, 'local_tol', 1E-20, 'kick_rank', 2, 'init_rank', r, 'max_rank', r);

ttsirt_t = TTSIRT(f, base, opt, 'var', deb);
%test_irt(f, ttsirt_t, 1E4)

%%
nsteps = 2^16;
init = eval_rt(ttsirt_t, deb.sample_x(:,1));
init = invert_cdf(model.ref, init);

tic; % negatively correlated pCN
out_pcn = pCN(@(z)pullback_potential_pcn(model,ttsirt_t,z), init, nsteps, log(10));
u = eval_cdf(model.ref, out_pcn.samples);
test_x = eval_irt(ttsirt_t, u);
toc
