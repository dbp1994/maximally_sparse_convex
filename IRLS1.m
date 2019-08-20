function [L2Err,L1Err,SE,x_n1,supp_n,avg_FP,avg_FN] = IRLS1(A,y,x,lambda_n,epsl)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x_0 = ones(size(A,2),1);
X = diag(x_0);
lambda_n = .5;
epsl = .1;
iter = 100;
tol = 1e-6;
stop_thres = 1e-4;
for i = 1:iter
    
    supp_xk = find(x_0~=0);
    X_k = X(supp_xk,supp_xk);
    A_k = A(:,supp_xk);
%     phi_xk = (2/(sqrt(3)*a_n))*((atan((1+2*a_n*abs(diag(X_k)))/sqrt(3)))-pi/6);
    
    B = (A_k'*A_k +2*lambda_n.*diag(1./diag(X_k)));
    b = A_k'*y;
    [x_n] = conjugate_grad(b,B,tol);
    
    x_n(find(abs(x_n)<.1)) = 0;
    x_n1 = zeros(size(A,2),1);
    x_n1(supp_xk,1) = x_n;
    
    X = diag(abs(x_n1)+epsl);
    if norm(x_n1 - x_0)<stop_thres        
        break;
    end
    x_0 = x_n1;    
end

% x_n1(find(abs(x_n1)<.001)) = 0;
supp_n = find(x_n1~=0);


[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(x_n1,x);
end




