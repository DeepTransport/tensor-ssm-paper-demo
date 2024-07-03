function u = matvec_prior_L(prior, v)
%MATVEC_PRIOR_L
% 
% L*v
%
% v ~ N(0, I)
% u ~ N(0, C)
%
% Tiangang Cui, 17/Jan/2014

switch prior.type
    case {'Field'}
        u = prior_cov_l(prior, v);
    case {'Basis'}
        u = prior.basis*v;
end

end
