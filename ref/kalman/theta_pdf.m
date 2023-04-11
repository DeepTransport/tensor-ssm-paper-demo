function value=theta_pdf(theta_all,C,Y)
% for all input theta, calculate the theoritical pdf of joint theta 1 2
% theta is already transformed to infinite domain
T=size(Y,2);
d=size(Y,1);
Y=reshape(Y,[],1);
value = zeros(1, size(theta_all, 2));
for k = 1:size(theta_all, 2)
    theta = theta_all(:, k);
    theta1 = 0.6 * normcdf(theta(1)) + 0.4;
    theta2 = 0.6 * normcdf(theta(2)) + 0.4;

    B=spdiags(-sqrt(1-theta1^2)*ones(d*T,1),-d,d*(T+1),d*(T+1));
    B=sparse(B+eye(d*T+d));
    C1=sparse([zeros(T*d,d),kron(eye(T),C)]);

    sig_inv=theta1^(-2)*(B'*B)+theta2^(-2)*(C1'*C1);
    mu=sig_inv\(C1'*Y*(theta2)^(-2));
    h=0.5*(mu'*sig_inv*mu-theta2^(-2)*(Y'*Y));
    value(k)=theta1^(-d*(T+1)) * det(B) *theta2^(-d*T)*det(sig_inv)^(-0.5)*exp(h) .* mvnpdf(theta);

end
end


