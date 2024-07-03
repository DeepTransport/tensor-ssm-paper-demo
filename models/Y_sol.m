classdef Y_sol
    properties
        model ssmodel % the underlying model
        poly % poly for ftt
        opt  % opt for ftt
        N % # of particles
        epd % expand the approximation domain
        SIRTs
        resample_time (1,:) = 0 % time steps of reapproximate and resample
        weight
        ESS_all (1,:)  % effective sample size for all time
        tau (1,1) = 1 % reapproximate threshold
        samples % contain all the samples of thetas, states of all time
        L % contain the matrix for linear transformation in all steps
        mu % conts. to fit in the FTT interval
        FTT_time (1, :)
        logmarginal_likelihood
    end


    methods
        function sol = Y_sol(model, poly, opt, N, epd, tau)
            sol.model = model;
            sol.poly = poly;
            sol.opt = opt;
            sol.N = N;
            sol.epd = epd;
            sol.weight = zeros(sol.N, model.T);
            sol.tau = tau;
            sol.samples = zeros(model.d+model.m, N, model.T+1);
            sol.L = zeros(model.d + model.m, model.d + model.m, model.T);
            sol.mu = zeros(model.d + model.m, model.T);
            sol.logmarginal_likelihood = [];
        end

        function sol = solve(sol)
            sol.samples(:,:,1) = priorsam(sol.model, sol.N);
            w = ones(1, sol.N)/sol.N;
            for t = 1:sol.model.T
                [sol.samples(:,:,t+1), w] = push_samples(sol.model, sol.samples(:,:,t), w, t);
                if ESS(w) < sol.tau*sol.N || t == sol.model.T
                    sprintf("FTT approximation at time %d because ESS =  %.2f / %d", t, ESS(w), sol.N)
                    covflag = 0; % when the ESS is too small, the approximate cov is not reliable, use previous one
                    if ESS(w)<30, covflag = 1; end
                    [sol, fun_post] = reapprox(sol, w, t, covflag);


                    % generate samples
                    z = rand(sol.model.d+sol.model.m, sol.N);
                    [r, f] = eval_irt(sol.SIRTs{t}, z);
                    sol.samples(:, :, t+1) = sol.L(:,:,t) * r + sol.mu(:, t);
                    w = exp(-fun_post(r))./exp(-f);
                end
                % compute ESS
                sol.ESS_all(t) = ESS(w);
                sol.weight(:, t) = w';
            end
        end

        function [sol, logfun_post] = reapprox(sol, w, t, covflag)
            T_r = sol.resample_time(end);

            if T_r == 0
                logfun_prior = @(x) -log(priorpdf(sol.model, x));
                L_old = sol.epd * eye(sol.model.d + sol.model.m);
            else
                logfun_prior = @(x) -log(eval_pdf(sol.SIRTs{T_r}, sol.L(:, :, T_r)\(x - sol.mu(:, T_r)))); 
                L_old = sol.L(:,:, T_r);
                % input true para states at the target time
                % output is potential function
            end
            [mu_temp, L_temp, ~] = computeL(sol.samples(:,:,t+1), w);
            
            L_temp = sol.epd * L_temp;
            if covflag == 1, L_temp = L_old; end

            sol.L(:, :, t) = L_temp;
            sol.mu(:, t) = mu_temp;
            const = 0;
            bound = sol.poly.oned_domains{1}.bound(2);
            z = rand(sol.model.d + sol.model.m, 1e3) *2 * bound - bound;
            fun_const = @(x) fun_into_sirt(sol, logfun_prior, t, T_r, L_temp*x + mu_temp, const);
            const = min(fun_const(z));
            logfun_post = @(x) fun_into_sirt(sol, logfun_prior, t, T_r, L_temp*x + mu_temp, const);

            % resample unweighted samples to serve as sample points and
            % debugger in TTSIRT
            samples_unweighted = L_temp\(datasample(sol.samples(:,:,t+1),...
                sol.N, 2, 'Weights', w) - mu_temp);
            samples_debug = samples_unweighted(:, 1:sol.N/2);
            samples_init = samples_unweighted(:, sol.N/2+1:end);
            deb = InputData(samples_init, samples_debug);
            tic
            sirt = TTSIRT(logfun_post, sol.poly, sol.opt, 'var', deb);
            sol.FTT_time = [sol.FTT_time ,toc]; 
            sol.SIRTs{t} = sirt;
            sol.resample_time = [sol.resample_time, t];
        end

        function logpdf = fun_into_sirt(sol, logf, t, T_r, x, const)
            logpdf = zeros(1, size(x, 2));

            for k =1: (t - T_r)
                % this is the x_values pull back to previous time
                if k == 1
                    x_cond = x;
                else
                    x_cond = Finv(sol.model, x_cond, t+1-k);
                end
                logpdf = logpdf - log(like(sol.model, x_cond, t+1-k)) - log(Finvgrad(sol.model, x_cond));
            end
            logpdf = logpdf + logf(Finv(sol.model ,x_cond, T_r));
            logpdf = logpdf - const;
        end
        

        function output = fullESS(sol, t)
            % generate samples of theta, X_T and hist each dimension
            % give the ESS of (theta, X_T|y_1:T}
            ftt = sol.SIRTs{t};
            z = rand(sol.model.d+sol.model.m, sol.N);
            [r, f] = eval_irt(ftt, z);
            end_samples = sol.L(:,:,t) * r + sol.mu(:, t);

            % compute the ESS
            logpdf = zeros(1, size(end_samples, 2));
            for k = 1: t
                y = sol.model.Y(:, t-k+1);
                % this is the x_values pull back to previous time
                if k == 1
                    x_cond = end_samples;
                else
                    x_cond = Finv(sol.model, x_cond, 1);
                end
                logpdf = logpdf - log(like(sol.model, x_cond, y)) - log(Finvgrad(sol.model, x_cond));
            end

            % prior
            logpdf = logpdf - log(priorpdf(sol.model, Finv(sol.model ,x_cond, 1)));
            logpdf = logpdf - min(logpdf);
            w = exp(-logpdf)./exp(-f);
            output = ESS(w);
        end
    end
end