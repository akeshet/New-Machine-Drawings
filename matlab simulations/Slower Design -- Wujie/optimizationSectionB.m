%
% Zeeman slower design calculator
% for Bec2 New Machine.
% Aviv Keshet, adapted from Peyman Ahmadi
%

%% Initialization and constants


clear
%----------------length segment-------------------------------
Delta = 0.1;                 % step
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
% Magnetic Momentum of Atom

% H= Mu. B where Mu= e g S /(2 m_e). For Na S for the last electron is
% hbar/2. g for an electron is 2. so Mu becomes, Mu= e hbar/(2 m_e)= Bohr
% Magneton.

Mu = 1.4*10^6;                   %Borh Magneton/h 

%In all of the equations Bohr magneton/hbar appears. So I think Mu should
%be multiplied by a factor of 2*pi. Check it regorously:
%A factor of 1/(2*pi) goes into delta to convert it to a frequency
%A factor of 1/(2*pi) goes into k to make kbar = h/Lambda


%% Desired field 
% ref: Thesis D.Durfee on Zeeman slower design 

% The zeeman slower is split into several pieces
% 1) a decreasing field slower which starts at some initial field and
% comes down to zero field
% 2) a zero-field "spin flip" region, where the bellows are located
% 3) a increasing field region with the same initial acceleration
% 4) a increasing field region at a lower acceleration


% 1
% * Decreasing field section *
% Theoretical curve (design for decelertion a = f*amax)
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

spinFlipLength = 13; % in unit of cm
SpinFlipFieldSection=zeros(spinFlipLength/Delta,1);


% 3
% * fast increasing section *
fastIncreasingLength=20; %cm

iFastIncreasingField = find(Bideal1<0);
iFastIncreasingFieldSection=iFastIncreasingField(1:fastIncreasingLength/Delta);
FastIncreasingFieldSection=Bideal1(iFastIncreasingFieldSection);

% 4 
% * slow increasing field section *

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



%DesignFieldProfile=[DecreasingFieldSection(:); SpinFlipFieldSection(:); FastIncreasingFieldSection(:); SlowIncreasingFieldSection(:)];
%z=0:length(DesignFieldProfile)-1;
%z=Delta*z;
%plot(z, DesignFieldProfile);
%title('Design Field');
%xlabel('Position (cm)');
%ylabel('Field (G)');



%% Achieved field
AchievedFieldProfile=0;

%select the regime of interest for optimization
IncreasingFieldSectionOP=transpose([FastIncreasingFieldSection(:); SlowIncreasingFieldSection(:);zeros(200,1)]);

z=-100:length(IncreasingFieldSectionOP)-1; 
z=Delta*(z+length(DecreasingFieldSection)+length(SpinFlipFieldSection));
MinDeviation=10000; % ensure the MinDeviation always decreases

IncreasingFieldSectionOP=transpose([zeros(100,1);IncreasingFieldSectionOP(:)]);

%Optimization Coil B

%Peymann's Design TurnsB=[76, 54, 40, 22, 7, 6]; WeakTurnsB=[7,5]; TurnsC=[7, 6, 0];
%[76, 54, 40, 22, 7, 6, 7,5,7, 6, 0]

%anand's thesis TurnsA = [78,53,37,17,7,0,6,5,7,6,0]; 

%Optimize design around
N=[76, 54, 40, 22, 7, 6, 7,5,7, 6, 0];


OptimizationRange=0;

for N1=N(1)-OptimizationRange:N(1)+OptimizationRange;
for N2=N(2)-OptimizationRange:N(2)+OptimizationRange;
for N3=N(3)-OptimizationRange:N(3)+OptimizationRange;
for N4=N(4)-OptimizationRange:N(4)+OptimizationRange;
for N5=N(5)-OptimizationRange:N(5)+OptimizationRange;
for N6=N(6)-OptimizationRange:N(6)+OptimizationRange;
for N7=N(7)-OptimizationRange:N(7)+OptimizationRange;
for N8=N(8)-OptimizationRange:N(8)+OptimizationRange;
for N9=N(9)-OptimizationRange:N(9)+OptimizationRange;
for N10=N(10)-OptimizationRange:N(10)+OptimizationRange;
for N11=max(0,N(11)-OptimizationRange):N(11)+OptimizationRange;
    
TurnsB = [N1,N2,N3,N4,N5,N6]; 
WeakTurnsB = [N7,N8];
TurnsC=[N9,N10,N11];
currentB=-33.9;
currentC=-115;
positionB = 100;

AchievedFieldProfileOptimal=0;

for i=1:length(TurnsB)
    AchievedFieldProfileOptimal=AchievedFieldProfileOptimal + bfield1(z, currentB, tube_OD + (2*i+1) * wire_thickness, positionB, TurnsB(i), -wire_thickness);
end

for i=1:length(WeakTurnsB)
    AchievedFieldProfileOptimal=AchievedFieldProfileOptimal + bfield1(z, currentB, tube_OD + (2*i+1) * wire_thickness, positionB - TurnsB(i)*wire_thickness, WeakTurnsB(i), -2*wire_thickness);
end

for i=1:length(TurnsC)
    AchievedFieldProfileOptimal= AchievedFieldProfileOptimal + bfield1(z, currentC, tube_OD + (1+2*i) * wire_thickness, positionB, TurnsC(i), wire_thickness);
end


AchievedFieldProfileOptimalSelect=AchievedFieldProfileOptimal;
IncreasingFieldSection1Select=IncreasingFieldSectionOP;

%for i=1:100;
%AchievedFieldProfileOptimal1(i)=0;
%DecreasingFieldSection1(i)=0;
%end

Deviation=0;
Derivative=0;
for i=100:300;
Derivative(i)=(AchievedFieldProfileOptimal(i+1)-IncreasingFieldSection1Select(i+1))-(AchievedFieldProfileOptimal(i)-IncreasingFieldSection1Select(i));
Deviation(i)=AchievedFieldProfileOptimalSelect(i)-IncreasingFieldSection1Select(i);
end

ToTDeviation=sqrt(std(abs(Deviation))^2+std(abs(Derivative))^2);

    if (ToTDeviation<MinDeviation)
        OptimalCoil = [N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11]
        MinDeviation = ToTDeviation
        AchievedFieldProfile=AchievedFieldProfileOptimal;
    end

end
end
end
end
end
end
end
end
end
end
end

subplot(2,1,1)
plot (z, AchievedFieldProfile,'r'); hold on;
plot (z, IncreasingFieldSectionOP,'b');

subplot(2,1,2)
plot (z, AchievedFieldProfile-IncreasingFieldSectionOP); hold on;

%plot (z, AchievedFieldProfile,'r'); hold on;
%plot (z, IncreasingFieldSectionOP,'b');
