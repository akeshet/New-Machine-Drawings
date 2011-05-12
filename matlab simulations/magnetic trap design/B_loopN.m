function B = B_loopN(r,z0,R,I);
%
% Calculates the magnetic field from an arbitrary current loop calculated from
% eqns (1) and (2) in Phys Rev A Vol. 35, N 4, pp. 1535-1546; 1987.
%
% arguments:
%
% * r is a position vector where the Bfield is evaluated: [x y z]
%   r is in cm
% * z0 is the location of the center of the loop in cm
% * R is the radius of the loop (cm)
% * I is in A
%
% returns:
% * B is a vector for the B field at the position r in gauss (CGS)
%
% Gustavo Telles (11/2004)
% 
x = r(1); y = r(2); z = r(3);
%
rho = sqrt( x.^2 + y.^2 );
%
k = sqrt( (4 * R * rho)/( (R + rho)^2 + (z-z0)^2) );
[K,E] = ellipke(k^2);
%
Bz = ( 1/sqrt((R+rho)^2 + (z-z0)^2) )*( K + E*(R^2-rho^2-(z-z0)^2)/((R-rho)^2+(z-z0)^2) );
%
pi=3.1421;
mu = 0.4*pi; % mu0 in CGS units (G*cm/A)
if rho ~= 0
    Brho = ( (z-z0)/(rho*sqrt((R+rho)^2 + (z-z0)^2)) )*( -K + E*(R^2+rho^2+(z-z0)^2)/((R-rho)^2+(z-z0)^2) );
    co= x/rho; si= y/rho;
    B = I*mu/(2*pi).*[co*Brho si*Brho Bz];
else
    B = I*mu/(2*pi).*[0 0 Bz];
end