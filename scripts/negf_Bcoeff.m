% Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
% Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 
%
% This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
% The software is distributed under the Lesser GNU General Public License (LGPL).
% ANGEL is free software: you can redistribute it and/or modify it under the terms 
% of the Lesser GNU General Public License v3 or later. ANGEL is distributed
% without any warranty; without even the implied warranty of merchantability or 
% fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.

function [x, B, p] = negf_Bcoeff(filebase)

[x, phi, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev([filebase,'_phi_n_p_J_Ec_Ev']);
[xphot, En, RphotG] = negf_read_xE([filebase,'_RphotGAL']);

nm = 1; % units
dens_unit = (1e-7*nm)^3;
time_unit = 1e-12;

nE = size(En,1);
nx = size(xphot,2);
Rtot = zeros(1,nx);
for xx=1:nx
    for ee=1:nE
        hw = En(ee);
        
        if ee==1
            dhw = 0.5*(En(2)-En(1));
        else if ee==nE
                dhw = 0.5*(En(nE-1)-En(nE));
             else
                dhw = 0.5*(En(ee+1)-En(ee-1));
             end
        end
        
        Rtot(xx) = Rtot(xx) + dhw * RphotG(ee,xx)/hw;
    end
end

nx2 = size(x,2);
B = zeros(1,nx2);
for xx=1:nx2
    
    dx = 1;
    if true % DO NOT DIVIDE BY dx! Rphot has units m-3 s-1!
        if xx==1
            dx = 0.5*(x(2)-x(1));
        else if xx==nx2
            dx = 0.5*(x(nx2)-x(nx2-1));
            else
                dx = 0.5*(x(xx+1)-x(xx-1));
            end
        end
    end
    
    
    % find xphot-indices corresponding to x(xx)
    % we assume that the Rphot-discretization is contained in 
    %fprintf(1,'xx=1: x=%e\n',x(xx));
    idx = 0;
    if xphot(1) > x(xx) || xphot(nx-1) < x(xx)
        continue;
    end
    for ii=1:nx
        if abs(xphot(ii)-x(xx)) < 1e-10
            idx = ii;
            break
        end
    end
    assert(idx~=0,'wew');
    %fprintf(1,'x(%d)=%e (n=%e,p=%e) has idx=%d, xphot=%e, Rtot=%e\n',xx,x(xx),n(xx),p(xx),idx,xphot(idx),Rtot(idx));
    
    %lambda = (xphot(idx)-x(xx-1)) / (x(xx)-x(xx-1)) % lambda=0 --> take x(idx-1)
    %Rmix   = lambda * Rtot(idx) + (1-lambda) * Rtot(idx-1);
    % n, p > 1e16[cm-3] = 1e-5[nm-3] --> lower limit 1e-10[nm-6] for n*p
    %B(xx) = Rtot(idx) / (n(xx) * p(xx) + (1e-13*nm^3)^2) / dx;
    B(xx) = Rtot(idx) / (n(xx) * p(xx) + (1e-17*nm^3)^2);
end

% units: everything is saved in simulator units
% --> n, p: nm-3; Rtot: nm3 ps-1 
% to get B[cm3 s-1] we need: 
% n, p are [nm-3]    --> multiply each by 1e-6*1e27=1e21 to have [cm-3]
% Rphot was [nm-3 ps-1], integration dhw Rphot/hw for Rtot is unit-neutral
% Rtot: [nm-3 ps-1] --> multiply by 1e-6*1e27 * 1e12 = 1e33 to get [cm-3 s-1]
% so altogether by 1e-21*1e-21*1e33=1e-9
B_unit = dens_unit / time_unit; 
p = semilogy(x,abs(B*B_unit));
