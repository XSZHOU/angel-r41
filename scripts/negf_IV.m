% Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
% Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 
%
% This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
% The software is distributed under the Lesser GNU General Public License (LGPL).
% ANGEL is free software: you can redistribute it and/or modify it under the terms 
% of the Lesser GNU General Public License v3 or later. ANGEL is distributed
% without any warranty; without even the implied warranty of merchantability or 
% fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.

function [V,I, p] = negf_IV(file, Imult)

if nargin < 2
    Imult = 1;
end

% read in file
fprintf(1,'Reading %s...',file);

data = load(file);
format short e;
fprintf(1,'  %d cols and %d rows. Taking first col as I, 2nd as V.\n',size(data,2),size(data,1));

V = data(:,1);
I = data(:,2);

Vmult = 1;
%Vmult = -1; Imult = -1;
%Vmult = 1; Imult = 1e14*1e12*1.602e-19;
%Vmult = 1; Imult = 1;

%p = semilogy(Vmult*V,Imult*I,'k','LineWidth',2);
p = plot(Vmult*V,Imult*I,'k','LineWidth',2);
xlabel('V [V]');
ylabel('J [A/cm^2]');
set(gcf,'Color','w');
fontname = 'Helvetica'; % Arial,Courier,Helvetica,Times,Bookman,lucidabright,...
set(gca,'FontSize',14,'FontName',fontname);
set(get(gca,'XLabel'),'FontName',fontname,'FontSize',20,'FontWeight','bold');
set(get(gca,'YLabel'),'FontName',fontname,'FontSize',20,'FontWeight','bold');
set(get(gca,'Title'),'FontName',fontname,'FontSize',24,'FontWeight','bold');
img = 'png';
res = '-r300';
%filename = [file,'.',img];   print(res,['-d',img],filename); fprintf(1,'figure 5 saved to %s\n',filename);
