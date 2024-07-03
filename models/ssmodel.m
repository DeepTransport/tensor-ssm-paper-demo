classdef ssmodel
    properties
        name string % the model name, to import the corresponding funs
        type = 0 % characterise if this model has state process error; 1 indicates no error and 2 indicates state error
        d (1,1)  % dimension for theta
        m (1,1)  % dimension for state
        n (1,1)  % dimension for observations
        T (1,1)  % time steps
        theta (:, 1) double % true value for theta
        X % states
        Y % observations
        pre % structure containing all prescribed paras
%         prior % the prior pdf.
%         F % push forward x_t, errors not included yet
%         st_process % x_{t+1} = F([theta;x_t]) % for Y_model, P = F
%         ob_process % y_t = G([theta;x_t])
%         like % likelihood pdf, input: theta, x_t, y_t
    end


    methods
        function model = ssmodel(name, d, m, n, T, varargin)
            model.name = name;
            model.d = d;
            model.m = m;
            model.n = n;
            model.T = T;
            if ~isempty(varargin)
                model.pre = varargin{1};
            end
%             model.theta = theta;
        end

        function model = complete(model)
            sams = priorsam(model, 1);
            model.X = zeros(model.m, model.T+1);
            model.Y = zeros(model.n, model.T);
            model.X(:,1) = sams(model.d+1 : end, 1);
            for t = 1: model.T
                model.X(:, t+1) = st_process(model, [model.theta; model.X(:,t)], t);
                model.Y(:, t) = ob_process(model, [model.theta; model.X(:,t+1)], t);
            end
        end
        
        function [sams, w] = push_samples(model, sams, w, t)
            % sams = [theta;X]
            % given sams at t-1 of size [d+m, :] and weights w, give weighted
            % samples at t and update weights using Y(:,t)
            % samsold = sams;

            X_sams = st_process(model, sams, t);
            sams = [sams(1:model.d,:); X_sams];
            % update weights in log scale and normalize
            w = log(w) + log(like(model, sams, t));
            w = w - max(w);
            w = exp(w);
            w(isnan(w)) = 0;
            w = w/sum(w);
        end

        function plot(model)
            for k = 1:model.m
                figure
                plot(0:model.T, model.X(k, :));
            end
        end
    end
end