function B = Bhelm(r, z0, a, coil, z);
%
% Gustavo Telles
% 
NN = length(coil); R = 0; in = 2.54;
N = length(z); drz = a*in + 0.1;
Bz_u = zeros(N,3); Bz_d = zeros(N,3);
By_u = zeros(N,3); By_d = zeros(N,3);
x = zeros(N,1); y = x;
for k=1:NN
    rk = r + (k-1)*drz;
    for j=1:coil(k)
        zj = z0 + (j-1)*drz;
        Bz_u = Bz_u + B_loopNt([x y z], zj, rk, 1);
%        By_u = By_u + B_loopNt([x z' y], zj, rk, 1);
        Bz_d = Bz_d + B_loopNt([x y z], -zj, rk, 1);
%        By_d = By_d + B_loopNt([x z' y], -zj, rk, 1);
    end
end
%
Bz = Bz_u(:,3) + Bz_d(:,3);
B = [z Bz]; %Br = By_u + By_d;