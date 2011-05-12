function B = Bz4h(a, z0, I, Z);
%
% Calculates the magnetic field of a square loop coil 
%
in = 2.54; % an inch in cm
mu = 0.4*pi; % mu0 in CGS units (G*cm/A)
z = Z(1):0.1:Z(2); N = length(z); Bz_u = zeros(N); Bz_d = Bz_u;
s = a/2; r_u = sqrt( s.^2 + (z-z0).^2); r_d = sqrt( s.^2 + (z+z0).^2);
D_u = sqrt( s.^2 + r_u.^2); D_d = sqrt( s.^2 + r_d.^2);
Bz_u = 2*mu*s*I ./ ( pi*in .* r_u .* D_u );
Bz_d = 2*mu*s*I ./ ( pi*in .* r_d .* D_d );
Bz = abs( Bz_u + Bz_d );
%
plot( z, Bz, '.-b');
ylabel('B (gauss)'); xlabel('z (in)'); 
title('B field of a square loop coil along the z direction');