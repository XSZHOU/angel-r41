% Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
% Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 
%
% This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
% The software is distributed under the Lesser GNU General Public License (LGPL).
% ANGEL is free software: you can redistribute it and/or modify it under the terms 
% of the Lesser GNU General Public License v3 or later. ANGEL is distributed
% without any warranty; without even the implied warranty of merchantability or 
% fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.

function negf_visualize(filebase, window)

if nargin < 2
    window = 1;
end

%close all
hold off
%close 1;
fig1 = figure(window);
nx = 4;
ny = 2;
set(fig1,'Color','w','Position',[50 50 nx*450 ny*400]);
  

% densities
[x, En, nE] = negf_read_xE([filebase,'_nE']);
[x, En, pE] = negf_read_xE([filebase,'_pE']);
subplot(ny,nx,1);contourf(x,En,(abs(nE)+abs(pE)+1e-15),100,'EdgeColor','none');title('n(x,E)+p(x,E)');
%surf(linspace(-10,50,size(nE,2)),linspace(Emin,Emax,size(nE,1)),real(pE),'EdgeColor','none');

% local DOS
[x, En, LDOS] = negf_read_xE([filebase,'_LDOS']);
%[x, En, LDOS_VB] = negf_read_xE([filebase,'_LDOS_VB']);
%[x, En, LDOS_k0] = negf_read_xE([filebase,'_LDOS_k0']);
subplot(ny,nx,2);contourf(x,En,(abs(LDOS)+1e-6),100,'EdgeColor','none');title('LDOS(x,E)');
%subplot(ny,nx,2);contourf(x,En,log(abs(LDOS)-abs(LDOS_VB)+1e-10),100,'EdgeColor','none');title('LDOS_{VB}(x,E)');
%subplot(ny,nx,ny*nx);contourf(x,En,(abs(LDOS_k0)+1e-10),100,'EdgeColor','none');title('LDOS_{k=0}(x,E)');
xidx = 1; %44
subplot(ny,nx,ny*nx);plot(En,LDOS(:,xidx));title(['LDOS(xx=',num2str(xidx),'=',num2str(x(xidx)),'nm,E)']);xlim([min(En) max(En)]);
subplot(ny,nx,ny*nx-1);plot(En);

% current
[x, En, jEe] = negf_read_xE([filebase,'_JEe']);
%[x, En, jEe2] = negf_read_xE([filebase,'_JEe2']);
[x, En, jEh] = negf_read_xE([filebase,'_JEh']);
% note: current plot is wrong because 1st row & col were assumed as axes
subplot(ny,nx,3); 
%contourf(x,En,log(abs(jEe)+abs(jEh)+1e-10),100,'EdgeColor','none');colorbar;title('J_{i,i+1}(x,E)');
contourf(x,En,jEe+jEh,100,'EdgeColor','none');colorbar;title('J_{i,i+1}(x,E)');
%subplot(ny,nx,ny*nx-1);contourf(x(2:62),En,jEe2(:,2:62)+1e-10,100,'EdgeColor','none');title('J^{e}_{2}(x,E)');
%subplot(ny,nx,ny*nx-1);contourf(x,En,(abs(jEe2)+1e-10),100,'EdgeColor','none');title('J^{e}_{2}(x,E)');
%subplot(ny,nx,ny*nx-1);contourf(x,En,(abs(jEh)+1e-15),100,'EdgeColor','none');colorbar;title('J^{h}_{i,i+1}(x,E)');

% phonon current
if exist([filebase,'_Jephon'])
    [x, En, jE3] = negf_read_xE([filebase,'_Jephon']);
    [x, En, jE4] = negf_read_xE([filebase,'_Jhphon']);
    subplot(ny,nx,4);contourf(x,En,(jE3+1e-10),100,'EdgeColor','none');title('J_{phon,e}(x,E)');colorbar;
subplot(ny,nx,7);contourf(x,En,(jE4+1e-10),100,'EdgeColor','none');title('J_{phon,h}(x,E)');colorbar;figure(window);
%figure(3);contourf(x,En,(jEh+1e-15),100,'EdgeColor','none');title('log J_{phon,h}(x,E)');colorbar;figure(window);
end

% photon recombination current
if exist([filebase,'_Jephot']) && 1
    [xJ, EnJ, jEphote]  = negf_read_xE([filebase,'_Jephot']);
    [xJ, EnJ, jEphoth]  = negf_read_xE([filebase,'_Jhphot']);
    subplot(ny,nx,5);contourf(xJ,EnJ,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');title('J_{\gamma}(x,E)');colorbar;
%    subplot(ny,nx,ny*nx);contourf(x,En,((jEphoth)+1e-15),100,'EdgeColor','none');title('J_{\gamma,h}(x,E)');%colorbar;
end

% photon spectrum
if exist([filebase,'_RphotGAL']) && 1
%    [x, En, RphotL]   = negf_read_xE([filebase,'_RphotLAK']);
%    [x, En, RphotLVB] = negf_read_xE([filebase,'_RphotLAK_VB']);
    [x, En, RphotG]   = negf_read_xE([filebase,'_RphotGAL']);
    [x, En, RphotGVB] = negf_read_xE([filebase,'_RphotGAL_VB']);
    
%    subplot(ny,nx,ny*nx);contourf(x,En,(abs(RphotL)+1e-15),100,'EdgeColor','none');title('R_{\gamma,Lake}(x,E)');colorbar;
%    subplot(ny,nx,ny*nx);contourf(x,En,(abs(RphotLVB)+1e-15),100,'EdgeColor','none');title('R_{\gamma,Lake}^{VB}(x,E)');colorbar;
%    subplot(ny,nx,ny*nx);contourf(x,En,(abs(RphotG)+1e-25),100,'EdgeColor','none');title('R_{\gamma,Gal}(x,E)');colorbar;

    subplot(ny,nx,6);contourf(x,En,(abs(RphotGVB)+1e-25),100,'EdgeColor','none');title('R_{\gamma,Gal}^{VB}(x,E)');colorbar;
    subplot(ny,nx,7); [hw,P,p] = negf_spectrum([filebase,'_spectrumGAL_VB']); 
    if (min(hw)~=max(hw)) 
        xlim([min(hw) max(hw)]);
    end
end

% contact current
%if exist([filebase,'_Jecont'])
%    [xJ, EnJ, jEconte]  = negf_read_xE([filebase,'_Jecont']);
%    [xJ, EnJ, jEconth]  = negf_read_xE([filebase,'_Jhcont']);
%    subplot(ny,nx,8);contourf(xJ,EnJ,(abs(jEconte)+abs(jEconth)+1e-10),100,'EdgeColor','none');title('J_{\gamma}(x,E)');colorbar;
%end

[x, phi, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev([filebase,'_phi_n_p_J_Ec_Ev']);

% plot J(x)
subplot(ny,nx,6); plot(x,J);

% formatting
labelsize = 16;
titlesize = 18;
for ii=1:5
    subplot(ny,nx,ii);
    
    % add band edges to plots
    hold on;plot(x,Ec,'k','LineWidth',2);plot(x,Ev,'k','LineWidth',2);hold off;
    
    xlabel('x[nm]');ylabel('E[eV]');
    set(gcf,'Color','w');
    fontname = 'Helvetica'; % Arial,Courier,Helvetica,Times,Bookman,lucidabright,...
    set(gca,'FontSize',14,'FontName',fontname);
    set(get(gca,'XLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight','bold');
    set(get(gca,'YLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight','bold');
    set(get(gca,'Title'),'FontName',fontname,'FontSize',titlesize,'FontWeight','bold');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 15 11]);
end

format = 'png';
res = '-r200';
%figure(1);filename = [filebase,'_npE.',format];    print(res,['-d',format],filename);fprintf(1,'figure 1 saved to %s\n',filename);
%figure(2);filename = [filebase,'_LDOS_k0.',format];print(res,['-d',format],filename);fprintf(1,'figure 2 saved to %s\n',filename);
%figure(3);filename = [filebase,'_LDOS.',format];   print(res,['-d',format],filename);fprintf(1,'figure 3 saved to %s\n',filename);
%figure(4);filename = [filebase,'_JE.',format];     print(res,['-d',format],filename);fprintf(1,'figure 4 saved to %s\n',filename);
