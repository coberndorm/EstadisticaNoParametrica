function Reg = RobustReg(X,Y,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin <3
    method = 'fastMCD'
end


Cov = RobustCov(X,Y,method);
XX = Cov.XX; XY = Cov.XY;

Betas = inv(XX)*XY;
meansX = mean(X,1);
Beta0 = mean(Y)-meansX*Betas;

y = Beta0 + X*Betas;

OrdinaryR2 = abs(1- sum((Y-y).^2)/sum((Y-mean(Y)).^2));
AjustedR2 =  abs(OrdinaryR2*((length(Y)-1)/(length(Y)-2)));

Reg.y = y;
Reg.OrdinaryR2 = OrdinaryR2;
Reg.AjustedR2 = AjustedR2;
Reg.Betas = Betas;
end