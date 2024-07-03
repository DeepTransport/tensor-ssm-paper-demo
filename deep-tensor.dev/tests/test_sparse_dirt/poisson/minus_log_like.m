function [mllkd, gmllkd, I] = minus_log_like(obs, d, mask)
%
% minus log likelihood for a given observable model output 
% gradient w.r.t. observable model output
% Fisher information matrix at observable model output

if isempty(mask)
    n_data = obs.n_data;
    data = obs.data;
else
    n_data = sum(mask);
    d = reshape(d(mask),[],1);
    data = reshape(obs.data(mask),[],1);
end

switch obs.like
    case {'normal'}
        if isscalar(obs.std)
            std = obs.std;
        else
            std = obs.std(mask);
        end
        misfit = (d - data)./std;
        mllkd  = 0.5*sum(misfit.^2); % minus log-likelihood
        gmllkd = misfit./std;
        if isscalar(std)
            I = speye(n_data, n_data)./std^2;
        else
            I = spdiags(1./std.^2, 0, n_data, n_data);
        end
    case {'poisson'}
        mllkd = - data'*log(d) + sum(d, 1); %- sum(log(factorial(obs.data(:))));
        gmllkd = 1 - data./d;
        d(d<eps) = eps;
        I = spdiags(1./d, 0, numel(d), numel(d));
        misfit = [];
    otherwise
        error('likelihood not implemented')
end

end