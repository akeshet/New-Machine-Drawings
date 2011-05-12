function B = helm (r, z0, I, a, coil, Z);
%
% Gustavo Telles
% 
BrzH(1.5, 1.8, 250, 1/8, [1 2 3 4 5 6], [-5 5]);
BrzH(6.8, 1.8, 330, 1/8, [6 6 6], [-5 5]);

x = -5:0.01:5; y = x; z = x;
N = length(z); N2 = round(N/2);
Bz_u = zeros(N,3); Bz_d = zeros(N,3);
By_u = zeros(N,3); By_d = zeros(N,3);
qc = [1 2 3 4 5 6]; fc = [6 6 6];
NN = length(qc); R = 0; drz = 2.54/8 + 0.1;
for k=1:NN
    rk = r + (k-1)*drz;
    for j=1:coil(k)
        zj = z0 + (j-1)*drz;
        for i=1:N
            Bz_u(i,:) = Bz_u(i,:) + [B_loopN([0 0 z(i)], zj, rk, I)];
            By_u(i,:) = By_u(i,:) + [B_loopN([0 z(i) 0], zj, rk, I)];
            Bz_d(i,:) = Bz_d(i,:) + [B_loopN([0 0 z(i)], -zj, rk, -I)];
            By_d(i,:) = By_u(i,:) + [B_loopN([0 z(i) 0], -zj, rk, -I)];
        end
    end
    R = R + coil(k)*rk;
end
Bz = zeros(N2,1); Bz = abs( Bz_u(N2:N,3) + Bz_d(N2:N,3) );
By = zeros(N2,1); By = abs( By_u(N2:N,2) + By_d(N2:N,2) );
Gz = mean( gradient(Bz,0.1) );     % average gradient in Z direction
Gr = mean( gradient(By,0.1) );     % average gradient in X,Y directions
a = a*in;                          % square size of the copper tubes in cm