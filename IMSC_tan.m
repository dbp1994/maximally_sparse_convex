function [L2Err,L1Err,SE,x_new,supp_n,avg_FP,avg_FN] = IMSC_tan(A,y,x,a_n,lambda_n,epsl,K)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x_0 = ones(size(A,2),1);
x_new = zeros(length(x),1);

X = diag(x_0);
iter = 100;
tol = 1e-6;
stop_thres = 1e-4;
supp_xk = find(x_0~=0);
X_k = X(supp_xk,supp_xk);
A_k = A(:,supp_xk);
phi_xk = (2./(sqrt(3).*a_n)).*((atan((1+2.*a_n.*abs(diag(X_k)))/sqrt(3)))-pi/6);
X_temp = (diag(X_k).*diag(X_k)./phi_xk);

for i = 1:iter   

%     B = (A_k'*A_k +2*lambda_n(supp_xk).*diag(X_temp));
    B = (A_k'*A_k +2*lambda_n.*diag(X_temp));
    X_temp = [];
    b = A_k'*y;
    [x_n] = conjugate_grad(b,B,tol);
    
    x_n(find(abs(x_n)<.001)) = 0;
    u_n1 = zeros(size(A,2),1);
    u_n1(supp_xk,1) = x_n;
    
    X = diag(abs(u_n1));
    
    supp_xk = find(u_n1~=0);
    X_k = X(supp_xk,supp_xk);
    A_k = A(:,supp_xk);

%     phi_xk = (2./(sqrt(3).*a_n(supp_xk))).*((atan((1+2.*a_n(supp_xk).*abs(diag(X_k)))/sqrt(3)))-pi/6);    
    phi_xk = (2./(sqrt(3).*a_n)).*((atan((1+2.*a_n.*abs(diag(X_k)))/sqrt(3)))-pi/6);    
%     X_temp = (diag(1./X_k).*diag(1./X_k).*phi_xk) + 1/(100*epsl);  

%     denom_x = (diag(X_k).*diag(X_k))./phi_xk;
%     X_temp = denom_x + .01;
    X_temp = (diag(1./X_k).*diag(1./X_k).*phi_xk) + 1/(100*epsl);  
%     X_temp = (diag(1./X_k).*diag(1./X_k).*phi_xk) + epsl;    
    if norm(u_n1 - x_0)<stop_thres
        break;
    end
    x_0 = u_n1;
    
end


% u_n1(find(abs(u_n1)<.001)) = 0;

figure(100)
hold on
satter_supp = find(u_n1~=0);
x_phix = u_n1(satter_supp);
A_k = A(:,satter_supp);
s = ((A_k'*(y - A*u_n1))./lambda_n);
% s = (A_k'*(y - A*u_n1))./lambda_n(satter_supp);
% grad_phi_x = sign(x_phix)./((a_n(satter_supp).*a_n(satter_supp).*x_phix.*x_phix) + a_n(satter_supp).*abs(x_phix)+1);
grad_phi_x = sign(x_phix)./((a_n.*a_n.*x_phix.*x_phix) + a_n.*abs(x_phix)+1);


plot(u_n1(satter_supp),s,'*')
hold on;
plot(x_phix,grad_phi_x,'o')
grid on

title('Scatter Plot of derivative of \phi(x) for \lambda = 2')
xlabel('x_n')
ylabel('H^T(y-Hx)/\lambda')
legend('observed','actual')
hold off

x_new(K) = u_n1;
supp_n = find(x_new~=0);

[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(x_new,x);



end




