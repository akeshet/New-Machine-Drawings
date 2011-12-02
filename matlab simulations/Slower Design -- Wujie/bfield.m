function [B]=bfield(z,current,diameter)
mu0 = 4*pi/10; 
B=(mu0*current/diameter)./(1 + (2*z/diameter).^2).^(3/2);