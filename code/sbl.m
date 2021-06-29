function [x, telap] = sbl(y,Phi, sigma, eps, Nmax)
%SBL Function to apply the sbl (EM) algorithm
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
tic
gamma = ones(size(Phi, 2),1);
x = zeros(size(Phi,2),1);
if Nmax == 0
    while norm(y-Phi*x) > eps
        s = inv(Phi'*Phi/sigma^2 + inv(diag(gamma)));
        x = s*(Phi'*y)/sigma^2;
        gamma = x.^2 + diag(s);
    end
else
   for k = 1:Nmax  
        s = inv(Phi'*Phi/sigma^2 + inv(diag(gamma)));
        x = s*(Phi'*y)/sigma^2;
        gamma = x.^2 + diag(s);
   end
end
telap = toc;
end