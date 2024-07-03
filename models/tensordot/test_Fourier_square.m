dom = BoundedDomain([-4,4]);
poly1 = ApproxBases(Fourier(10), dom, 1);
poly2 = ApproxBases(Fourier(20), dom, 1);

c1 = zeros(21, 1);
c2 = zeros(21, 1);
c1(1) = 1;
c1(2) = 1/sqrt(2);
c2(1) = 1;
c2(2) = 1/sqrt(2);
coef = {c1,c2};


% Fourier_square(coef, poly1.oneds{1}, poly2.oneds{1})';

%%% test passed