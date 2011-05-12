function B = Bz4L(a, z0, I, Z);
%
% Calculates the magnetic field of a square loop coil 
%
in = 2.54; % an inch in cm
mu = 0.4*pi; % mu0 in CGS units (G*cm/A)
z = Z(1):0.05:Z(2); N = length(z); Bz = zeros(N);
s = a/2; r = sqrt( s.^2 + (z-z0).^2); D = sqrt( s.^2 + r.^2);
Bz = 2*mu*s*I ./ ( pi*in .* r .* D );
B = [z' Bz'];
%
%plot( z, Bz, '.-b');
%ylabel('B (gauss)'); xlabel('z (in)'); 
%title('B field of a square loop coil along the z direction');