function HL = computeHL2(f1, f2)
N = 40;
index = linspace(0.4, 1 ,N);
index2 = index + 1e-2*eye(1,N) - 1e-2* flip(eye(1,N));
[X, Y] = meshgrid(index,index);
[X_index2, Y_index2] = meshgrid(index2, index2);
vgrid = [X(:),Y(:)]';
vgrid2 = [X_index2(:),Y_index2(:)]';

p_thm(:, :) = reshape(f1(vgrid2), N, N);
p_thm(:, :) = p_thm(:, :)./(nansum(nansum(p_thm(:, :)))* (index(2)-index(1))^2);

p_ftt(:, :) = reshape(f2(vgrid2), N, N);
p_ftt(:, :) = p_ftt(:, :)./(nansum(nansum(p_ftt(:, :)))* (index(2)-index(1))^2);

HL = sqrt(0.5 * nansum(nansum((sqrt(p_thm(:,:))-sqrt(p_ftt)).^2))*(index(2)-index(1))^2);

end