function [R] = SDP_rn(F0,alp_min)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

p = size(F0,1);
c = -ones(p,1);
b = alp_min*ones(p,1);
A = int8(eye(p));
At = zeros(p*p,p);

btt = -c;
ctt = [-b; vec(F0)];
Ai = [-1, zeros(1,p-1)];
for i = 1:p    
    Fi = diag(circshift(Ai,i-1));
    At(:,i) = int8(-vec(Fi));
end
Att = [-A; At];
K.l = size(A,1);
K.s = size(F0,1);
[x, y, info] = sedumi(double(Att),btt,ctt,K);
R = diag(y);
end

