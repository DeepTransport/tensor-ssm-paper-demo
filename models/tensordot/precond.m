function [C, Sigmak] = precond(Sigma1,Sigma2)
% sigma1 is of d+m, sigma2 is of m
% C is the linear preconditioner, Sigmak is CSigma2C'

% m = size(Sigma2, 1);
% d = size(Sigma1, 1) - size(Sigma2, 1);
% C1 = inv(chol(Sigma1(1:d, 1:d))');
% [C2,~] = eig(Sigma1(d+1:d+m,d+1:d+m)\Sigma2);
% C2 = C2';
% % C2 = real(diag(sqrt(diag(C2*Sigma1(d+1:d+m, d+1:d+m)*C2')))\C2);
% C2 = real(diag(sqrt(diag(C2*Sigma1(d+1:d+m, d+1:d+m)*C2')))\C2); 
% C = [C1, zeros(d,m);zeros(m,d),C2];
% 
% Sigmak = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m = size(Sigma2, 1);
% d = size(Sigma1, 1) - size(Sigma2, 1);
% Sigma2 = [zeros(d), zeros(d,m);zeros(m,d), Sigma2];
% [U1, D1] = eig(Sigma1);
% [U2, D2] = eig(Sigma2);
% [U3, ~] = eig(sqrt(D1)\(U1'*Sigma2*U1/sqrt(D1)));
% C = real((U3'/sqrt(D1))*U1');
% 
% Sigmak = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m = size(Sigma2, 1);
% d = size(Sigma1, 1) - size(Sigma2, 1);
% C1 = inv(chol(Sigma1(1:d, 1:d))');
% [C2,~] = eig(Sigma1(d+1:d+m,d+1:d+m)\Sigma2);
% C2 = C2';
% % C2 = real(diag(sqrt(diag(C2*Sigma1(d+1:d+m, d+1:d+m)*C2')))\C2);
% C2 = real(diag(sqrt(diag(C2*Sigma1(d+1:d+m, d+1:d+m)*C2')))\C2); 
% C = [C1, zeros(d,m);zeros(m,d),C2];
% 
% Sigmak = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size(Sigma2, 1);
d = size(Sigma1, 1) - size(Sigma2, 1);
Sigma2full = [zeros(d), zeros(d,m);zeros(m,d),Sigma2];
Sigma3 = Sigma1\Sigma2full;
[C, ~] = eig(Sigma3);
C = C';
C(1:d, 1:d) = inv(chol(Sigma1(1:d,1:d)))';
% C2 = [eye(d), zeros(d,m); zeros(m,d), C2'];

% C2 = real(diag(sqrt(diag(C2*Sigma1(d+1:d+m, d+1:d+m)*C2')))\C2);
C = real(diag(sqrt(diag(C*Sigma1*C')))\C); 

Sigmak = C(d+1:d+m,d+1:d+m)*Sigma2*C(d+1:d+m,d+1:d+m)';
Sigmak = diag([zeros(1, d), diag(Sigmak)']);

end

