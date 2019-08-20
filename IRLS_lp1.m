function [L2Err,L1Err,SE,x_n,supp_n,avg_FP,avg_FN] = IRLS_lp1(A,b,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x_0 = ones(size(A,2),1);
X = diag(x_0);
iter = 100;
tol = 1e-4;
supp_xk = find(x_0~=0);
for i = 1:iter
    
    supp_xk = find(x_0~=0);
    X_k = X(supp_xk,supp_xk);
    A_k = A(:,supp_xk);
    x_n = ((X_k.*X_k)*A_k')*pinv(A_k*(X_k.*X_k)*A_k')*b;
%     x_n = ((X_k.*X_k)*A_k')*(A_k*(X_k.*X_k)*A_k')\b;
    
    x_n(find(abs(x_n)<.01)) = 0;
    x_n1 = zeros(size(A,2),1);
    x_n1(supp_xk,1) = x_n;  

    
    X = diag(abs(x_n1).^(.6));
    if norm(x_n1 - x_0)<tol        
        break;
    end
    x_0 = x_n1;    
end

x_n1(find(abs(x_n1)<.05)) = 0;
supp_n = find(x_n1~=0);
supp_x = find(x~=0);

[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(x_n1,x);
end

