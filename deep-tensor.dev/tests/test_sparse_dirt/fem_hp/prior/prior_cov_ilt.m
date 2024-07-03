function v = prior_cov_ilt(prior, u)
%COV_INVLU 
%
% whitening transformation
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case {'MRF'}
        t = zeros(size(u));
        t(prior.pm,:) = prior.Lm'\u;
        v = prior.K*t;
    case {'GP'}
        v = prior.L'\u;
end

end

% Given a permutation matrix P and the corresponding vector p, we have
%         P'u = u(p) and PP' = I or inv(P) = P'
%
% prior.K => K, prior.M => M, both are symmetric
% Pk: permutation matrix s.t. Pk' K Pk = K(pk, pk) = Lk Lk'
% Pm: permutation matrix s.t. Pm' M Pm = M(pm, pm) = Lm Lm'
%
% C = inv(Q) = inv(K*inv(M)*K) = inv(K) M inv(K)
%   => define Lc = inv(K) L and Lc' = L' inv(K) with LL' = M
%
% Pm' M Pm = M(pm, pm) = Lm Lm' <=> PmPm' M PmPm' = Pm Lm Lm' Pm'
%   <=> M = Pm Lm Lm' Pm' => define L = Pm Lm
%
% u = L v = Pm Lm v <=> Pm' u =  Pm'Pm Lm v <=> u(pm) = Lm v
% u = L'v = Lm'Pm'v <=> u = Lm' v(pm)
%
% inv(K) = inv(PkPk' K PkPk') = inv(Pk')inv(Pk' K Pk)inv(Pk)
%        = Pk inv(Pk' K Pk) Pk' = Pk inv(Lk Lk') Pk'
%        = Pk inv(Lk') inv(Lk) Pk'
%
% u = inv(K) v = Pk inv(Lk') inv(Lk) Pk' v
%  <=> Pk' u = inv(Lk') inv(Lk) Pk' v
%  <=> u(pk) = Lk'\(Lk\v(pk))
%
% want to evaluate inv(Lc') v 
%
% Lc' = L' inv(K) <=> inv(Lc') = K inv(L')
%
% inv(L') = inv(Lm'Pm') = inv(Pm') inv(Lm') = Pm inv(Lm')
%
% Lc' v = K Pm inv(Lm') v
%    1. u(pm) = Lm'\ v
%    2. w = K u
