function [mllkd, mlp] = minus_log_post_uni(model, obs, prior, v, masks)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

if nargin < 5
    masks = [];
end

% prior, assuming the input v is wsoltten parameters
mlp = 0;
% map to correlated Gaussian
u = matvec_prior_L(prior, v) + prior.mean_u;
% forward solve
n = size(u,2);
if isscalar(masks) && masks == -1
    mllkd = zeros(obs.n_data,n);
    for i = 1:n
        sol = forward_solve(model, u(:,i));
        mllkd(:,i) = minus_log_like_dirt(obs, sol.d);
    end
else
    mllkd = zeros(1,n);
    for i = 1:n
        sol = forward_solve(model, u(:,i));
        mllkd(i) = minus_log_like(obs, sol.d, masks);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
