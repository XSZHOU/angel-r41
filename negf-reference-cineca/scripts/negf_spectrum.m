% Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
% Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 
%
% This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
% The software is distributed under the Lesser GNU General Public License (LGPL).
% ANGEL is free software: you can redistribute it and/or modify it under the terms 
% of the Lesser GNU General Public License v3 or later. ANGEL is distributed
% without any warranty; without even the implied warranty of merchantability or 
% fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.

function [hw,P,p] = negf_spectrum(file, mult)

% multiplicator is optional :-)
if nargin < 2
    mult=1;
end

fprintf(1,'Reading %s...',file);

data = load(file);
format short e;
fprintf(1,'   file has %d columns and %d rows.\n',size(data,2),size(data,1));

% file is in NEGF-internal units eV,ps 
% --> to go from ps-1 to s-1, multiply by 1e12
% --> to go from eV to J, multiply by ec
% --> ec is needed because P = int_d(hw) P(hw) and hw is in eV NO!!!
ec = 1.601e-19;
% conv = 1e12*ec;
% conv = 1;
conv = 1e12;

if (size(data,1)==2 || size(data,1)==3)
    hw = data(1,:);
    P  = data(2,:);
else
    P  = data(:,1);
    hw = linspace(0,size(P,1),size(P,1));
end
hwmult = 1; Pmult = -1;

p = plot(hwmult*hw , Pmult*mult*P*conv,'k','LineWidth',2);

% formatting
xlabel('hw [eV]');
ylabel('P [s^{-1}]');
set(gcf,'Color','w');
fontname = 'Helvetica'; % Arial,Courier,Helvetica,Times,Bookman,lucidabright,...
set(gca,'FontSize',14,'FontName',fontname);
set(get(gca,'XLabel'),'FontName',fontname,'FontSize',20,'FontWeight','bold');
set(get(gca,'YLabel'),'FontName',fontname,'FontSize',20,'FontWeight','bold');
set(get(gca,'Title'),'FontName',fontname,'FontSize',24,'FontWeight','bold');

