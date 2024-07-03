function z = ttrankest_amen(bases, rho, N)
% Solve the cubic equation to estimate the #sweeps matching N evaluations
% rho is the kick rank in amen, equal to the initial rank
n = mean(cellfun(@(b)numel(b.nodes), bases.oneds));
d = numel(bases.oneds);

a = (d-2)*n*rho^2;
b = n*rho + 5*(d-2)*n*rho^2;
c = (9/2)*n*rho + (43/6)*(d-2)*n*rho^2;
q = 5*n*rho*(1+(d-2)*rho) - N;

% Crude initial guess
z = (N/(d*n))^(1/3)/rho+1;
% Newton iteration
for iter=1:20
    R = (a/3)*z^3 + (b/2)*z^2 + c*z + q;
    J = a*z^2 + b*z + c;
    z = z - R/J;
end

z = round(max(z+1,1));
end
