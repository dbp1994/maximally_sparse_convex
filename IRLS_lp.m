function [L2Err,L1Err,SE,u_n1,supp_n,avg_FP,avg_FN] = IRLS_lp(A,y,x,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x_0 = ones(size(A,2),1);
X = diag(x_0);
iter = 1000;
tol = 1e-6;
stop_thres = 1e-4;
q = 1-p/2;
lambda_n = .5;
epsl = .1;
supp_xk = find(x_0~=0);
X_k = X(supp_xk,supp_xk);
A_k = A(:,supp_xk);
X_temp = (abs(diag(X_k)).^q);

for i = 1:iter
    
    B = (A_k'*A_k +p*lambda_n.*diag(1./X_temp));
    X_temp = [];
    b = A_k'*y;
    [x_n] = conjugate_grad(b,B,tol);
    
    x_n(find(abs(x_n)<.1)) = 0;
    u_n1 = zeros(size(A,2),1);
    u_n1(supp_xk,1) = x_n;
    
    X = diag(abs(u_n1));
    
    supp_xk = find(u_n1~=0);
    X_k = X(supp_xk,supp_xk);
    A_k = A(:,supp_xk);
    x_i = diag(X_k);
    
    X_temp = (((x_i.*x_i)+epsl).^q);
    if norm(u_n1 - x_0)<stop_thres        
        break;
    end
    x_0 = u_n1;    
end

u_n1(find(abs(u_n1)<.001)) = 0;
supp_n = find(u_n1~=0);


[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(u_n1,x);
end




