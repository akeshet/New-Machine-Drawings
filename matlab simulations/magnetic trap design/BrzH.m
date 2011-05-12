function B = BrzH(r, z0, I, a, coil, Z);
%
% Calculates the magnetic field from an arbitrary coil geometry 
%
% arguments:
%
% * r is the starting (smaller) radius of the coil loop given in cm
% * z0 is the location of the center of the loop in cm
% * I is in A
% ...
% returns:
% * lots of information and plots...
%
% Gustavo Telles
% 
z = Z(1):0.05:Z(2); N = length(z); N2 = round(N/2);
Bz_u = zeros(N,3); Bz_d = zeros(N,3);
By_u = zeros(N,3); By_d = zeros(N,3);
NN = length(coil); R = 0; in = 2.54;
drz = (a+0.01)*in;
for k=1:NN
    rk = r + (k-1)*drz;
    for j=1:coil(k)
        zj = z0 + (j-1)*drz;
        for i=1:N
            Bz_u(i,:) = Bz_u(i,:) + [B_loopN([0 0 z(i)], zj, rk, I)];
            By_u(i,:) = By_u(i,:) + [B_loopN([0 z(i) 0], zj, rk, I)];
            Bz_d(i,:) = Bz_d(i,:) + [B_loopN([0 0 z(i)], -zj, rk, I)];
            By_d(i,:) = By_d(i,:) + [B_loopN([0 z(i) 0], -zj, rk, I)];
        end
    end
    R = R + coil(k)*rk;
end
%
% final assignments
%
if a>=3/16
    f=3/15;                        % 12 l/m, assumed flux for the cooling agent for the 3/16" tube
    id=in/8;                       % square inner size of the copper tube for the 3/16" tube
else
    f=1/20;                        % 3 l/m, assumed flux for the cooling agent for the 1/8" tube
    id=in/16;                      % square inner size of the copper tube for the 1/8" tube
end
%
Bz = zeros(N2,1); Bz = abs( Bz_u(N2:N,3) + Bz_d(N2:N,3) );
By = zeros(N2,1); By = abs( By_u(N2:N,2) + By_d(N2:N,2) );
Gz = mean( gradient(Bz,0.1) );     % average gradient in Z direction
Gr = mean( gradient(By,0.1) );     % average gradient in X,Y directions
a = a*in;                          % square size of the copper tubes in cm
d = 1.0;                           % cooling agent mass density, 1.0 kg/l (H2O here)
C = 4186;                          % cooling agent heat capacity, H2O here
dRdT = 0.00393;                    % resistance temperature coefficient
r0 = 1.7e-6;                       % copper wire resistance at room temperature (ohms*cm)
Lc = 4*pi*R+400;                   % coils full length (winding + water circuit connections)
R0 = r0*Lc/(a^2-id^2);             % coils total resistance
P0 = R0*I^2;                       % voltage drop on the coils at room temperature
dT = 1/( (d*C*f)/P0 - dRdT );      % temperature change in cooling agent due to heat transfer
rT = 1.7e-6*(1+0.004*dT);          % copper wire resistance w/ temperature correction (ohms*cm)
Rc = rT*Lc/(a^2 - id^2);           % coils' resistance (T corrected)
Vc = Rc*I;                         % voltage drop per coil (T corrected)
Pc = Vc*I;                         % power to be dissipated (T corrected)
Vigbt=1.8; Vlds=1.3; Vxtra=0.4;    % extra voltage drops in the system
Vt = Vc + Vigbt + Vlds + Vxtra;    % total voltage drop in the system
Pt = Vt*I;                         % total power to be dissipated (T corrected)
Dz = max(coil)*drz + a/2;           % coil height
Di = 2*(r - a/2);                  % inner coil diameter
Do = 2*(r + (NN-1/2)*a);           % outer coil diameter
%
% outputs & plotting
%
fprintf('\rCOILS fact sheet:\r\t\tdB/dz = %0.1f G/cm\r\t\tdB/dr = %0.1f G/cm\r\t\tRc = %0.4f ohms\r\t\tVc = %0.1f V\r\t\tPc = %1.3f kW',Gz,Gr,Rc,Vc,Pc/1000);
fprintf('\r\t\tinner D: %1.2fcm\r\t\touter D: %1.2fcm\r\t\tcoil height: %1.2fcm\r\r',Di,Do,Dz);
fprintf('\r\rTOTAL fact sheet:\r\t\tLt = %0.1f m\r\t\tVt = %0.1f V\r\t\tPt = %1.3f kW\r\t\tdT = %1.1fºC\r\r',Lc/100,Vt,Pt/1000,dT);
%
% Plots B as function of z and y (plane)
%
plot( z(:), abs( Bz_u(:,3) + Bz_d(:,3) ), '.-b', z(:), abs( By_u(:,2) + By_d(:,2) ), '-r');
ylabel('B (gauss)'); xlabel('r,z (cm)'); 
title('Coil B field the along the r (red) and z (blue) direction');