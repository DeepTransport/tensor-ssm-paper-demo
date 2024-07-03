function value = theta_pdf(model, theta, t)
% for all input theta, calculate the theoritical pdf of joint theta 1 2
% theta is in 0.4 to 1

Y = model.Y(:, 1:t);
m = model.m;
n = model.n;
Y = reshape(Y,[],1);
value = zeros(1, length(theta));
C = model.pre.C;


for k = 1:size(theta, 2)
    
    theta1 = theta(1, k);
    theta2 = theta(2, k);

    A = spdiags(-sqrt(1-theta1^2)*ones(m*t,1),-m,m*(t+1),m*(t+1));
    A = A+eye(m*t+m);
    A(1:m, 1:m) = theta1*eye(m);
    H = sparse([zeros(t*m,m),kron(eye(t),C)]);

    sig_inv = theta1^(-2)*(A'*A) + theta2^(-2)*(H'*H);
    mu = sig_inv\(H'*Y*(theta2)^(-2));
    h = 0.5*(mu'*sig_inv*mu-theta2^(-2)*(Y'*Y));
    value(k) = theta1^(-m*(t+1)) * det(A) *theta2^(-n*t)*det(sig_inv)^(-0.5)*exp(h) /0.36;

end
end



