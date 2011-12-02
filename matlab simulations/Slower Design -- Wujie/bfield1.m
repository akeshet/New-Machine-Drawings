
function [B1]=bfield1(z, current, diameter, firstposition, numturns, spacing)

mu0 = 4*pi/10; 
B1=0;
for i=1:numturns
B1=B1+bfield((z - firstposition - (i - 1)*spacing), current, diameter);
end