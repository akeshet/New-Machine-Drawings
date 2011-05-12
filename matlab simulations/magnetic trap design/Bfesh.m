function y = Bfesh( Z, d, Iq, If)
Qc = [1 2 3 4 5 6];
Fc = [6 6 6];
zz = Z(1):d:Z(2); N = length(zz); 
Bq = Bhelm(1.5, 1.8, 1/8, Qc, zz');
Bf = Bhelm(6.96, 1.8, 1/8, Fc, zz');
B = abs( Iq.*Bq(:,2) + If.*Bf(:,2) );
N0 = find( zz==0 );
mean( gradient( B(N0:N),d ) )  % MEAN gradient
mean( del2( B(N0:N),d ) )      % MEAN curvature
Bz = [Bq(:,1) B];
plot( Bz(:,1), Bz(:,2), '.-r');