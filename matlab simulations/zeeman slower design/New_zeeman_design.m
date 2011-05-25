%
% Zeeman slower design calculator
% for Bec2 New Machine.
% Aviv Keshet, adapted from Peyman Ahmadi
%

%% Initialization and constants


clear
%----------------length segment-------------------------------
Delta = 0.1;
Epsilon = 0.01;
%--------------------------------------------------------------
%--------------------Slower Segments--------------------------
% 1 = Decreasing field coil
% 2 = Increasing field coil
% 3 = second part of increasing field coil
% 4 = compensation coil
%----------------Constants-------------------------------------
cm = 1;
mm = 0.1*cm;
inch = 2.54*cm;
mu0 = 4*pi/10;                    %In Gauss-cm/Amp
hbar=1.055*10^(-34);              %J.s

%----------------Coil properties-------------------------------
wire_thickness = (0.1424)*inch;         %new coil dimension in cm
wire_id = wire_thickness - 2*.032*inch;  %Inner diameter of the copper wire
wire_od = wire_thickness;                %Outer diameter of the copper wire
tube_OD = 2.54*cm;                     %Outer diameter
resistivity101 = 1.7*10^-6;       %square tubing
resistivity122 = 2.03*10^-6;      %round tubing
%--------------------------------------------------------------

%------------Atomic Na properties---------------------------------
Lambda = 589*10^(-9);            %Laser light wavelength (m)
mass = 3.82*10^(-26);            %Na mass (Kg)
Gamma = (2*pi)*9.7*10^6;          %Na linewidth (Hz)
vmax=950;                        %Max velocity for oven T of 553 K, m/s
S = 1.5; %Isat = 6.26 mW/cm^2, I = 9.6...
%--------------------------------------------------------------

%------------Atomic Rb properties---------------------------------
%Lambda = 780*10^(-9);            %Laser light wavelength (m)
%mass = (87/23)*3.82*10^(-26);            %Na mass (Kg)
%Gamma = (2*pi)*6*10^6;          %Na linewidth (Hz)
%--------------------------------------------------------------

% H= Mu. B where Mu= e g S /(2 m_e). For Na S for the last electron is
% hbar/2. g for an electron is 2. so Mu becomes, Mu= e hbar/(2 m_e)= Bohr
% Magneton.

Mu = 1.4*10^6;                   %Borh Magneton/h 

%In all of the equations Bohr magneton/hbar appears. So I think Mu should
%be multiplied by a factor of 2*pi. Check it regorously:
%A factor of 1/(2*pi) goes into delta to convert it to a frequency
%A factor of 1/(2*pi) goes into k to make kbar = h/Lambda


%% Desired field

% The zeeman slower is split into several pieces
% 1) a decreasing field slower which starts at some initial field and
% comes down to zero field
% 2) a zero-field "spin flip" region, where the bellows are located
% 3) a increasing field region with the same initial acceleration
% 4) a increasing field region at a lower acceleration


% 1
% * Decreasing field section*
% Theoretical curve (design for a = f*amax)
f = 0.6;
amax = (hbar*(2*pi)/Lambda)/mass*(Gamma/2);    % Maximum acceleration in the slower
detuning = -1000*10^6;                         %Cooling light detuning (negetiv sign is absorbed in B-field Eq)   
Bzero = (detuning+vmax/Lambda)/Mu; %=720.0     %Initial magnetic field at the entrance of the slower     
Bfinal = detuning/Mu;                    %B-field at the end of the slower

% Slower length (not including 0-field gap region (bellows)
len = vmax^2/(2*f*amax);                 %Slower length (m).
len = 100*len;                           %Slower length in cm

%Defining the length variable
z1=0:Delta:len;
Bideal1 =(Bzero - Bfinal)*sqrt(1 - z1./len) + Bfinal;   %B-field Gauss

iDecreasingFieldSection = find(Bideal1>0);
DecreasingFieldSection = Bideal1(iDecreasingFieldSection);

%length(zDecreasingFieldSection)*Delta

% 2
% * Spin flip section *
% (field = zero for some length)

spinFlipLength = 13; %cm
SpinFlipFieldSection=zeros(spinFlipLength/Delta,1);


% 3
% fast increasing section
fastIncreasingLength=20; %cm

iFastIncreasingField = find(Bideal1<0);
iFastIncreasingFieldSection=iFastIncreasingField(1:fastIncreasingLength/Delta);
FastIncreasingFieldSection=Bideal1(iFastIncreasingFieldSection);

% 4 
% slow increasing field section


% Theoretical curve, but designed for a lower acceleration. 
f = 0.4;
amax = (hbar*(2*pi)/Lambda)/mass*(Gamma/2);    % Maximum acceleration in the slower
detuning = -1000*10^6;                         %Cooling light detuning (negetiv sign is absorbed in B-field Eq)   
Bzero = (detuning+vmax/Lambda)/Mu; %=720.0     %Initial magnetic field at the entrance of the slower     
Bfinal = detuning/Mu;                    %B-field at the end of the slower
% Slower length (not including 0-field gap region (bellows)
len = vmax^2/(2*f*amax);                 %Slower length (m).
len = 100*len;                           %Slower length in cm

%Defining the length variable
z1=0:Delta:len;
Bideal2 =(Bzero - Bfinal)*sqrt(1 - z1./len) + Bfinal;   %B-field Gauss


lastvalue = FastIncreasingFieldSection(length(FastIncreasingFieldSection));
iSlowIncreasingFieldSection = find(Bideal2<lastvalue);

SlowIncreasingFieldSection=Bideal2(iSlowIncreasingFieldSection);



DesignFieldProfile=[DecreasingFieldSection(:); SpinFlipFieldSection(:); FastIncreasingFieldSection(:); SlowIncreasingFieldSection(:)];
z=0:length(DesignFieldProfile)-1;
z=Delta*z;
plot(z, DesignFieldProfile);
title('Design Field');
xlabel('Position (cm)');
ylabel('Field (G)');



%% Achieved field
AchievedFieldProfile=0;

%Coil A
TurnsA = [132,119, 104, 89, 73, 56, 39, 21, 11, 11, 0, 0]; %Peyman's optimized-by-hand slower design
currentA=15;

for i=1:length(TurnsA)
    AchievedFieldProfile=AchievedFieldProfile+bfield1(z, currentA, tube_OD + (1+2*i)*wire_thickness, 0, TurnsA(i), wire_thickness);
end

%Coil B
%This design uses a few turns of "double-spaced coils".
TurnsB=[76, 54, 40, 22, 7, 6];
WeakTurnsB=[7,5];
currentB=-33.9;
positionB = 100;

for i=1:length(TurnsB)
    AchievedFieldProfile=AchievedFieldProfile + bfield1(z, currentB, tube_OD + (2*i+1) * wire_thickness, positionB, TurnsB(i), -wire_thickness);
end

for i=1:length(WeakTurnsB)
    AchievedFieldProfile=AchievedFieldProfile + bfield1(z, currentB, tube_OD + (2*i + 1) * wire_thickness, positionB - TurnsB(i)*wire_thickness, WeakTurnsB(i), -2*wire_thickness);
end

figure(1)
plot (z, AchievedFieldProfile, z, DesignFieldProfile);
figure(2)
plot (z, AchievedFieldProfile-transpose(DesignFieldProfile));
