%%
deb = InputData(test_x(:,1:debug_size), test_x(:,(debug_size+1):debug_size*2));

%%
% Ranges of domain
a = -3*ones(1,model.D+2);
b = 3*ones(1,model.D+2);

%%
p = 16; % order
r = 8; % rank

%%
% Vectors of grid points
polys = cell(model.D+2,1);
domains = cell(model.D+2,1);
for i=1:model.D+2
    polys{i} = Chebyshev2nd(p);
    domains{i} = BoundedDomain([a(i),b(i)]);
end
base = ApproxBases(polys, domains);

ttsirt_b = TTSIRT(f, base, opt, 'var', deb);

test_irt(f, ttsirt_b, test_x)

%%
% Vectors of grid points
polys = cell(model.D+2,1);
domains = cell(model.D+2,1);
for i=1:model.D+2
    polys{i} = Chebyshev1st(p);
    domains{i} = AlgebraicMapping(1);
end
base = ApproxBases(polys, domains);

ttsirt_t = TTSIRT(f, base, opt, 'var', deb);
test_irt(f, ttsirt_t, test_x)
