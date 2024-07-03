function output = theta_pdf_all(model)
% for all input theta, calculate the theoritical pdf of joint theta 1 2
% theta is in true domain 0.4 to 1, with uniform prior
% output is a N*N*T array, which contains the true pdf for all time steps

N = 40;
index = linspace(0.4, 1 ,N);
index2 = index + 1e-2*eye(1,N) - 1e-2* flip(eye(1,N));
% [X, Y] = meshgrid(index,index);
[X_index2, Y_index2] = meshgrid(index2, index2);
% [X_index, Y_index] = meshgrid(index, index);
% vgrid = [X(:),Y(:)]';
vgrid2 = [X_index2(:),Y_index2(:)]';
output = zeros(N, N, model.T);

for t = 1:model.T
    p_thm = theta_pdf(model, vgrid2, t);
    p_thm = p_thm./(nansum(nansum(p_thm))* (index(2)-index(1))^2);
    output(:,:,t)  = reshape(p_thm, N, N);
end

end



