function output = est_kf(sol, ind)
% input is the sol structure and the time ind indicating when to estimate
% p(\theta)

N = 40;
index = linspace(0.4, 1 ,N);
index2 = index + 1e-2*eye(1,N) - 1e-2* flip(eye(1,N));
% [X, Y] = meshgrid(index,index);
[X_index2, Y_index2] = meshgrid(index2, index2);
% [X_index, Y_index] = meshgrid(index, index);
% vgrid = [X(:),Y(:)]';
vgrid2 = [X_index2(:),Y_index2(:)]';
vgrid_tranform = norminv((vgrid2-0.4)/0.6);
output = zeros(N, N, length(ind));

for k = 1:length(ind)
    t = ind(k);
    if sol.sqr == 1
        sirt = sol.SIRTs{t};
        f4 = @(x) eval_pdf(sirt, sol.L(1:2,1:2,t)\(x - sol.mu(1:2,t)));
    else
        ftt = sol.FTTs{t};
        ftt_int = int_block(ftt, sol.model.d+1 : sol.model.d+sol.model.m);
        ftt_int.oneds = ftt_int.oneds(1:sol.model.d);
        ftt_int.oned_domains = ftt_int.oned_domains(1:sol.model.d);
        f4 = @(x) eval_reference(ftt_int, sol.L(1:2,1:2,t)\(x - sol.mu(1:2,t)));
    end
    
    f = f4(vgrid_tranform)./normpdf(vgrid_tranform(1,:))./normpdf(vgrid_tranform(2,:)); % 
    f = reshape(f, N, N);
    output(:, :, k) = f./(nansum(nansum(f))* (index(2)-index(1))^2);
end

end