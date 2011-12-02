%This program calculates the ideal B-field in a Zeeman slower for 
%Na atoms. Then a current and winding pattern is calculated which 
%closely produces the predicted B-field pattern.

%% Initialization

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


%% Ideal field curve

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
startposition=len*(1-(-Bfinal/(Bzero-Bfinal))^2);      %creating an offset in z such that at z=0 B-filed=0
z=z1 - startposition;
Bideal =(Bzero - Bfinal)*sqrt(1 - z1./len) + Bfinal;   %B-field Gauss


%% Decreasing-field section of the slower

fprintf('\r\t\tlength of the decreasing slower: %1.2fcm',startposition);

current = 15;

%Turns = [135, 121, 107, 93, 77, 61, 43, 25, 15, 8, 0, 0]; %Number of turns anand's thesis
Turns = [132,119, 104, 89, 73, 56, 39, 21, 11, 11, 0, 0]; %Peyman's optimized-by-hand slower design

firstposition=-startposition;

%bfield1 function format (z-position, current through the coil, coil
%diameter, first loop of the coil position, numturns, wire diameter); 
Solenoidfield=0;

numLayers=length(Turns);
for i=1:numLayers
    Solenoidfield= Solenoidfield + ...
        bfield1(z, current, tube_OD + (1+2*i)*wire_thickness, firstposition, Turns(i), wire_thickness);
end


subplot(4,2,1)
plot(z,Bideal)
hold on             
plot(z,Solenoidfield,'r')
hold off
ylabel('Field (G)'); xlabel('z (cm)'); 
title('Ideal(blue) and created (red) B-field for decreasing slower');

subplot(4,2,2)
plot(z,Solenoidfield-Bideal)
xlim([-45 0]);
ylabel('Field (G)'); xlabel('z (cm)');
title('Error of decreasing slower');

%Shift the ideal b-field 7 cm according to the anand's thesis
z6=z(503:753);
b1=Bideal(503:753);


%% Increasing-field section of the slower
f = 0.4; 
len2 = vmax^2/(2*f*amax);
len2 = 100*len2;                  %Convert to cm 
Bfinal = detuning/Mu;             %Note this is larger than the max field due to finite atom velocity at end *)
z2=0:Delta:len2+20*cm;            %An offset of 20 cm is given to calculate the net B-fiel at the position of the MOT 

%B-field calculation in Gauss
Bideal2=zeros(1,size(z2,2));
for i=1:size(z2,2)
    if (z2(i)<len2)
        Bideal2(i)=(Bzero - Bfinal)*sqrt(1 - z2(i)/len2) + Bfinal;
    else
        Bideal2(i)=0;
    end
end     

startposition2=len2*(1-(-Bfinal/(Bzero-Bfinal))^2);    %creating an offset in z such that at z=0 B-filed=0

z3=z2-startposition2;
biggining=find(z3>0,1)

z0 = 35;
current1 = 33.9;
current2 = 115;       %109;                    %up to 100 A

little = 0;

%bfield1 function format (z-position, current through the coil, coil
%diameter, first loop of the coil position, numturns, diameter of the coil);
%
%Turns2 = [78,53,37,17,7,0,6,5]; %Ananth's thesis

Turns2=[76, 54, 40, 22, 7, 6];
WeakTurns2=[7,5];
% the beginning section of the increasing field slower starts with windings
% that are spaced by twice the normal wire thickness
% (to have a better gradual field turn-on)
% this is the classic technique
% numerically, this is simulated by the "WeakTurns" which are spaced at
% double the usual spacing
Incfield1=0;

numLayers=length(Turns2);
for i=1:numLayers
    Incfield1= Incfield1 + bfield1(z3, current1, tube_OD + (2*i+1) * wire_thickness, z0 - little, Turns2(i), -wire_thickness);
end

numLayers=length(WeakTurns2);
for i=1:numLayers
    Incfield1=Incfield1 + bfield1(z3, current1, tube_OD + (2*i + 1) * wire_thickness, z0 - little - Turns2(i)*wire_thickness, WeakTurns2(i), -2*wire_thickness);
end

%% Compensation coils for slower
Turns3=[7 6 0];
numLayers=length(Turns3);

shim = 0.025*inch
Incfield2=0;

for i=1:numLayers
    Incfield2= Incfield2 + bfield1(z3, current2, tube_OD + (1+2*i) * wire_thickness + shim, z0, Turns3(i), wire_thickness);
end

Incfield = Incfield1 + Incfield2;

subplot(4,2,3)
plot(z6,b1)
hold on
plot(z3(biggining:size(z3,2)),-Incfield(biggining:size(z3,2)),'r') 
hold on
ylabel('G'); xlabel('z (cm)');
title('Ideal  B-field of increasing slower');

Bshift=Bideal2(biggining+size(z6,2):size(z3,2));
z33=z3(biggining+size(z6,2):size(z3,2));

%Only plotting the part with f=0.4 after 25 cm 
plot(z33(102:size(z33,2))-11.75,Bshift(102:size(Bshift,2)),'g')  
hold on
plot(z3(biggining+size(z6,2):size(z3,2)),-Incfield(biggining+size(z6,2):size(z3,2)),'r') 
hold off


size(Incfield(biggining+size(z6,2):size(z3,2)));


subplot(4,2,4)
plot(z6(16:size(z6,2)),-Incfield(biggining:biggining+size(z6,2)-1-15)-b1(16:size(b1,2)),'r')
%plot(z3(biggining+102size(z6,2):size(z3,2)),-IncfieldPartII-Bshift(102:size(Bshift,2)-1),'r') 
hold on
diff2 = -Bshift(102:size(Bshift,2)) - Incfield(1014:1014+size(Bshift,2)-102);
zdiff2 = z33(102:size(Bshift,2))-11.75;
plot(zdiff2,diff2)
hold off
ylim([-20 20])

ylabel('G'); xlabel('z (cm)');
title('(Ideal-created)  B-field of increasing slower');

%-------------Compensation Coil-------------------------------------
current3 = 0;
cs2 = 0.1*cm;
cs2 = wire_thickness;
extra = 1*mm;
comppos = z0 + 2*wire_thickness + extra;
ztrap = comppos + 6*cs2 + wire_thickness/2 + 5*cm;                          %MOT position?
compod = tube_OD + 7*wire_thickness;

%Turns3= [6, 4]; 
Turns3 = [0, 0];            %Anand's thesis
Compfield= bfield1(z3, current3, compod, comppos, Turns3(1), cs2) +...
           bfield1(z3, current3, compod + 2*cs2, comppos + 2*wire_thickness, Turns3(2), cs2);

%plot(z3(biggining:size(z3,2)),Compfield(biggining:size(z3,2)),'r') 


%Total field of the second slower
Slower2= -Incfield + 0*Compfield;

subplot(4,2,5)
plot(z6,b1)
hold on
%plot(z3(biggining:size(z3,2)),-Incfield(biggining:size(z3,2)),'r')
plot(z3(biggining:size(z3,2)),Slower2(biggining:size(z3,2)),'r')
hold on
plot(z33(108:size(z33,2))-12.4,Bshift(108:size(Bshift,2)),'g')  
hold off

ylabel('B-field'); xlabel('z (cm)'); 
title('Ideal(blue) and created (blue) B-field the second part of the slower with compensation coil');

subplot(4,2,6)
plot(z3(biggining:size(z3,2)), Bideal2(biggining:size(z3,2))-Slower2(biggining:size(z3,2)))
xlim([z3(biggining) 16])
ylim([-2 4])
ylabel('G'); xlabel('z (cm)');
title('(Ideal-created)  B-field of increasing slower + compensation');




spacing = 7.5*cm;                   %offsetting the second slower

fprintf('\r\t\tend of the second part of the slower: %1.2fcm',z0+spacing);
%Plotting every thing
subplot(4,2,[7,8])
plot(z(1:503),Solenoidfield(1:503),'r')
hold on
plot(z(1:503),Bideal(1:503))
hold on
plot(z6+7,b1)
hold on
plot(z3+spacing-12,Bideal2,'k') 
hold on
plot(z3(biggining:size(z3,2))+spacing,Slower2(biggining:size(z3,2)),'r')
hold off
xlim([-50 50]);
ylabel('B-field'); xlabel('z (cm)'); 
title('Ideal(blue) and created (blue) B-field the slower');
%Coil Power and Voltage

%Lengths of wire in cm
%Resistances of each wire

%First section: Decreasing field coil
i=1:length(Turns);
CoilLengths1 =pi*(tube_OD + (-1+2*i)*wire_thickness).*Turns;
TotalLength1 = 200 + sum(CoilLengths1);
Resistance1 = resistivity101*TotalLength1/(wire_od^2-wire_id^2);
TotalVolts1 = Resistance1*current;
TotalPower1 = current*TotalVolts1;
fprintf('\r\t\tcoil length of the decreasing slower in cm \r\t\t')
for i=1:length(Turns)
    fprintf('%1.2f \r\t\t',CoilLengths1(i))
end
fprintf('\r\t\ttotal voltage and power requiered for decreasing slower: %1.2f and %1.2f \r',TotalVolts1,TotalPower1)

%Second section: Increasing field coil, part 1
N = length(Turns2);
j=1:N;
CoilLengths2(1:(N-1)) =pi*(tube_OD + (-1+2*j(1:(N-1)))*wire_thickness).*Turns2(1:(N-1));
CoilLengths2(N) = 2*pi*(tube_OD + wire_thickness).*Turns2(N);
TotalLength2 = 200 + sum(CoilLengths2);
Resistance2 = resistivity101*TotalLength2/(wire_od^2-wire_id^2);
TotalVolts2 = Resistance2*current1;
TotalPower2 = current1*TotalVolts2;

fprintf('\r\t\tcoil length of the increasing coils in cm (part I)\r\t\t')
for i=1:length(Turns2)
    fprintf('%1.2f \r\t\t',CoilLengths2(i))
end
fprintf('\r\t\ttotal voltage and power requiered for increasing field coil (part I): %1.2f and %1.2f \r',TotalVolts2,TotalPower2)

%Part 2
k=1:length(Turns3);
CoilLengths3 =pi*(tube_OD + (-1+2*k)*wire_thickness).*Turns3;
TotalLength3 = 200 + sum(CoilLengths3);
Resistance3 = resistivity101*TotalLength3/(wire_od^2-wire_id^2);
TotalVolts3 = Resistance3*current2;
TotalPower3 = current2*TotalVolts3;

fprintf('\r\t\tcoil length of the increasing coils in cm (part II)\r\t\t')
for i=1:length(Turns3)
    fprintf('%1.2f \r\t\t',CoilLengths3(i))
end
fprintf('\r\t\ttotal voltage and power requiered for increasing field coil (part II): %1.2f and %1.2f \r',TotalVolts3,TotalPower3)

%Compensation Coil
TotalLength4 = pi*(compod+compod+2*cs2);
fprintf('\r\t\tcoil length of compensation coil in cm: %1.2f',TotalLength4)
Resistance4 = resistivity101*TotalLength4/(wire_od^2-wire_id^2);
TotalVolts4 = Resistance4*current3;
TotalPower4 = current3*TotalVolts4;
fprintf('\r\t\ttotal voltage and power requiered for compensation coil (part II): %1.2f and %1.2f \r',TotalVolts4,TotalPower4)


