function L1 = computeL1(f1, f2)
N = 40;
index = linspace(0.4, 1 ,N);
index2 = index + 1e-2*eye(1,N) - 1e-2* flip(eye(1,N));
[X, Y] = meshgrid(index,index);
[X_index2, Y_index2] = meshgrid(index2, index2);
vgrid = [X(:),Y(:)]';
vgrid2 = [X_index2(:),Y_index2(:)]';

p_thm(:, :) = reshape(f1(vgrid2), N, N);
p_thm(:, :) = p_thm(:, :)./(nansum(nansum(p_thm(:, :)))* (index(2)-index(1))^2);

p_ftt(:, :) = max(reshape(f2(vgrid2), N, N), 0);
p_ftt(:, :) = p_ftt(:, :)./(nansum(nansum(p_ftt(:, :)))* (index(2)-index(1))^2);

L1 = .5 * nansum(nansum(abs(p_thm-p_ftt))) * (index(2)-index(1))^2;


end