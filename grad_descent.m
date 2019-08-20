function [L2Err,L1Err,SE,new_x,supp_n,avg_FP,avg_FN] = grad_descent(A,y,x,lambda_n,a_n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
prev_x = zeros(size(A,2),1);
AtA = A'*A;
[~,eign_val] = eig(AtA);
% eigen_min = max(diag(eign_val));
stop_thres = 1e-16;
for i = 1:length(y)
    grad_phi_x = sign(prev_x)./((a_n.*a_n.*prev_x.*prev_x) + a_n.*abs(prev_x)+1);
    grad_f = (AtA*prev_x - A'*y) + lambda_n.*grad_phi_x;
    new_x = prev_x - .01*grad_f;
    if norm(prev_x - new_x)<stop_thres   
        
        break;
        
    else
        prev_x = new_x;  
%         phi_xk = (2./(sqrt(3).*a_n)).*((atan((1+2.*a_n.*abs(new_x))/sqrt(3)))-pi/6);
%         f(i) = .5*(norm(y - A*new_x)^2) + sum(lambda_n.*phi_xk);
    end
    
    
    
end
new_x(find(abs(new_x)<.1)) = 0;
supp_n = find(new_x~=0);

[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(new_x,x);

end

