function iwce_prepare_pics()

function cbar()
  % comment the next line out/in to switch between colorbar inclusion
  %colorbar;
end


close all

filebase = '/home/steiger/src/release/amd64/NEGF/results/iwce_GOOD/iwce_Source-';
% fixed-voltage files for n(x,E), J(x,E), LDOS(x,E)
turnonvolt = 1.5;
voltA = [filebase,sprintf('%5.3f',turnonvolt - 0.1),'V'];
voltB = [filebase,sprintf('%5.3f',turnonvolt),'V'];
voltC = [filebase,sprintf('%5.3f',turnonvolt + 0.1),'V'];
% files for spectra, P-I
%volts = [1.100,1.200,1.300,1.400,1.500,1.525,1.550,1.575];
%volts = [1.100,1.200,1.300,1.400,1.500,1.525,1.550,1.575,1.600];
volts = [1.100,1.200,1.250,1.300,1.350,1.400,1.450,1.500,1.550,1.600];

% NEGF-internal units: ps, ec 
% --> power (per cm2 cross-section) is saved in ec/ps rather than J/s
% --> spectrum (per cm2 cross-section) is saved in 1/ps rather than 1/s
spin = 2;
negf_power_fac    = spin * 1; 
negf_spectrum_fac = spin * 1; 

%aquatree = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-12/qw-iwce';
%aquatree = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e12_crate0/qw-iwce';
%aquatree = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-11/qw-iwce';
%aquatree = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-11_crate0/qw-iwce';
%aquatree = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-12_fermi/qw-iwce';
aquatree = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-13_fermi_mh0.2/qw-iwce';
% files for spectra
aquavolts = [1.500,1.525,1.550];
% file for P-V
aquafile = [aquatree '_Source1.850V.xy'];

% AQUA simulation was for 10nm*10nm cross-section; want 1cm2
aqua_power_fac    = 1e12; 
aqua_spectrum_fac = 1e12 / 1.602e-19; 

make_negfpics = 1;

%cm = 1-gray;
%cm = 1-hot;
%cm = [0.9583 1 1;0.7986 1 1;0.6389 1 1;0.4792 1 1;0.3194 1 1;0.1597 1 1;0 1 1;0 0.96 1;0 0.92 1;0 0.88 1;0 0.84 1;0 0.8 1;0 0.76 1;0 0.72 1;0 0.68 1;0 0.64 1;0 0.6 1;0 0.56 1;0 0.52 1;0 0.48 1;0 0.44 1;0 0.4 1;0 0.36 1;0 0.32 1;0 0.28 1;0 0.24 1;0 0.2 1;0 0.16 1;0 0.12 1;0 0.08 1;0 0.04 1;0 0 1;0 0 0.9688;0 0 0.9375;0 0 0.9062;0 0 0.875;0 0 0.8438;0 0 0.8125;0 0 0.7812;0 0 0.75;0 0 0.7188;0 0 0.6875;0 0 0.6562;0 0 0.625;0 0 0.5938;0 0 0.5625;0 0 0.5312;0 0 0.5;0 0 0.4688;0 0 0.4375;0 0 0.4062;0 0 0.375;0 0 0.3438;0 0 0.3125;0 0 0.2812;0 0 0.25;0 0 0.2188;0 0 0.1875;0 0 0.1562;0 0 0.125;0 0 0.09375;0 0 0.0625;0 0 0.03125;0 0 0];
%cm = [1 1 1;0.9813 0.9795 0.976;0.9625 0.959 0.9521;0.9438 0.9386 0.9281;0.9251 0.9181 0.9041;0.9063 0.8976 0.8802;0.8876 0.8771 0.8562;0.8688 0.8566 0.8322;0.8501 0.8362 0.8083;0.8314 0.8157 0.7843;0.7977 0.7836 0.7554;0.7641 0.7516 0.7264;0.7305 0.7195 0.6975;0.6968 0.6875 0.6686;0.6632 0.6554 0.6397;0.6295 0.6234 0.6107;0.5959 0.5913 0.5818;0.5623 0.5592 0.5529;0.5286 0.5272 0.5239;0.495 0.4951 0.495;0.4804 0.4806 0.4804;0.4659 0.466 0.4659;0.4513 0.4514 0.4513;0.4368 0.4369 0.4368;0.4222 0.4223 0.4222;0.4076 0.4078 0.4076;0.3931 0.3932 0.3931;0.3785 0.3786 0.3785;0.364 0.3641 0.364;0.3494 0.3495 0.3494;0.3348 0.3349 0.3348;0.3203 0.3204 0.3203;0.3057 0.3058 0.3057;0.2912 0.2913 0.2912;0.2766 0.2767 0.2766;0.2621 0.2621 0.2621;0.2475 0.2476 0.2475;0.2329 0.233 0.2329;0.2184 0.2184 0.2184;0.2038 0.2039 0.2038;0.1893 0.1893 0.1893;0.1747 0.1748 0.1747;0.1601 0.1602 0.1601;0.1456 0.1456 0.1456;0.131 0.1311 0.131;0.1165 0.1165 0.1165;0.1019 0.1019 0.1019;0.08735 0.08738 0.08735;0.07279 0.07281 0.07279;0.05823 0.05825 0.05823;0.04368 0.04369 0.04368;0.02912 0.02913 0.02912;0.01456 0.01456 0.01456;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0];
cm = [1 1 1;0.9404 0.9561 1;0.8808 0.9122 1;0.8212 0.8682 1;0.7616 0.8243 1;0.702 0.7804 1;0.6651 0.7449 0.9725;0.6281 0.7094 0.9451;0.5912 0.6739 0.9176;0.5542 0.6384 0.8902;0.5173 0.6029 0.8627;0.4803 0.5674 0.8353;0.4434 0.5319 0.8078;0.4064 0.4964 0.7804;0.3695 0.4609 0.7529;0.3325 0.4254 0.7255;0.2956 0.3899 0.698;0.2586 0.3544 0.6706;0.2217 0.3189 0.6431;0.1847 0.2834 0.6157;0.1478 0.2479 0.5882;0.1108 0.2124 0.5608;0.07389 0.1769 0.5333;0.03695 0.1414 0.5059;0 0.1059 0.4784;0.05 0.1506 0.4545;0.1 0.1953 0.4306;0.15 0.24 0.4067;0.2 0.2847 0.3827;0.25 0.3294 0.3588;0.3 0.3741 0.3349;0.35 0.4188 0.311;0.4 0.4635 0.2871;0.45 0.5082 0.2631;0.5 0.5529 0.2392;0.55 0.5976 0.2153;0.6 0.6424 0.1914;0.65 0.6871 0.1675;0.7 0.7318 0.1435;0.75 0.7765 0.1196;0.8 0.8212 0.09569;0.85 0.8659 0.07176;0.9 0.9106 0.04784;0.95 0.9553 0.02392;1 1 0;0.9649 0.9474 0;0.9298 0.8947 0;0.8947 0.8421 0;0.8596 0.7895 0;0.8246 0.7368 0;0.7895 0.6842 0;0.7544 0.6316 0;0.7193 0.5789 0;0.6842 0.5263 0;0.6491 0.4737 0;0.614 0.4211 0;0.5789 0.3684 0;0.5438 0.3158 0;0.5087 0.2632 0;0.4737 0.2105 0;0.4386 0.1579 0;0.4035 0.1053 0;0.3684 0.05263 0;0.3333 0 0];

f = 0;
if make_negfpics
% ---------------------------------------------------
% densities, LDOS, currents, Rphot at voltage A
% ---------------------------------------------------
[x, En, nE]              = negf_read_xE([voltA,'_nE']);
[x, En, pE]              = negf_read_xE([voltA,'_pE']);
[x, En, LDOS]            = negf_read_xE([voltA,'_LDOS']);
%[x, En, LDOS_VB]         = negf_read_xE([voltA,'_LDOS_VB']);
%[x, En, LDOS_k0]        = negf_read_xE([voltA,'_LDOS_k0']);
[xj, Ej, jEe]            = negf_read_xE([voltA,'_JEe']);
[xj, Ej, jEh]            = negf_read_xE([voltA,'_JEh']);
[xphon, Ephon, jEphone]  = negf_read_xE([voltA,'_Jephon']);
[xphon, Ephon, jEphonh]  = negf_read_xE([voltA,'_Jhphon']);
[xphot, Ephot, jEphote]  = negf_read_xE([voltA,'_Jephot']);
[xphot, Ephot, jEphoth]  = negf_read_xE([voltA,'_Jhphot']);
[xR, ER, RphotG]         = negf_read_xE([voltA,'_RphotGAL']);
%[xR, ER, RphotGVB]       = negf_read_xE([voltA,'_RphotGAL_VB']);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(nE+pE)+1e-15),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('n(x,E)+p(x,E)');
ylim([-0.7 1.3]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(nE+pE)+1e-15),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('n(x,E)+p(x,E)');
ylim([-0.3 1.8]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(LDOS)+1e-6),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('LDOS(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xj,Ej,(abs(jEe)+abs(jEh)+1e-10),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('J(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xj,Ej,log(abs(jEe)+abs(jEh)+1e-6),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);%title('log J(x,E)');
%f=f+1; figure(f); hold off;colormap(cm);
%contourf(xphon,Ephon,(jEphone+jEphonh+1e-10),100,'EdgeColor','none');
%xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]); %title('J_{phon,e}(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xphot,Ephot,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');
xlim([13 32]); xlabel('x[nm]');ylabel('E[eV]');colorbar;%title('J_{\gamma}(x,E)');
ylim([turnonvolt-2.05 turnonvolt-1.05]);set(gca,'YTick',[-1.0 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(xphot,Ephot,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');
xlim([13 32]); xlabel('x[nm]');ylabel('E[eV]');colorbar;%title('J_{\gamma}(x,E)');
ylim([turnonvolt-1.05 turnonvolt-0.05]);set(gca,'YTick',[-1.0 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(xR,   ER,   (abs(RphotG)+1e-25),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('$\hbar\omega$ [eV]','interpreter','latex');
xlim([13 32]); ylim([1.33 1.70]);colorbar;%title('R_{\gamma,Gal}(x,E)');

num_voltA = f;

[x, phi, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev([voltA,'_phi_n_p_J_Ec_Ev']);
for ii=1:num_voltA-1
    figure(ii);
    hold on;plot(x,Ec,'k','LineWidth',2);plot(x,Ev,'k','LineWidth',2);hold off;
end

% ---------------------------------------------------
% densities, LDOS, currents, Rphot at voltage B
% ---------------------------------------------------
[x, En, nE]              = negf_read_xE([voltB,'_nE']);
[x, En, pE]              = negf_read_xE([voltB,'_pE']);
[x, En, LDOS]            = negf_read_xE([voltB,'_LDOS']);
%[x, En, LDOS_VB]         = negf_read_xE([voltB,'_LDOS_VB']);
%[x, En, LDOS_k0]        = negf_read_xE([voltB,'_LDOS_k0']);
[xj, Ej, jEe]            = negf_read_xE([voltB,'_JEe']);
[xj, Ej, jEh]            = negf_read_xE([voltB,'_JEh']);
[xphon, Ephon, jEphone]  = negf_read_xE([voltB,'_Jephon']);
[xphon, Ephon, jEphonh]  = negf_read_xE([voltB,'_Jhphon']);
[xphot, Ephot, jEphote]  = negf_read_xE([voltB,'_Jephot']);
[xphot, Ephot, jEphoth]  = negf_read_xE([voltB,'_Jhphot']);
[xR, ER, RphotG]         = negf_read_xE([voltB,'_RphotGAL']);
%[xR, ER, RphotGVB]       = negf_read_xE([voltB,'_RphotGAL_VB']);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(nE+pE)+1e-15),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('n(x,E)+p(x,E)');
ylim([-0.7 1.3]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(nE+pE)+1e-15),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('n(x,E)+p(x,E)');
ylim([-0.3 1.8]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(LDOS)+1e-6),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('LDOS(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xj,Ej,(abs(jEe)+abs(jEh)+1e-10),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('J_{i,i+1}(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xj,Ej,log(abs(jEe)+abs(jEh)+1e-6),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);%title('log J(x,E)');
%f=f+1; figure(f); hold off;colormap(cm);
%contourf(xphon,Ephon,(jEphone+jEphonh+1e-10),100,'EdgeColor','none');
%xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]); %title('J_{phon,e}(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xphot,Ephot,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');
xlim([13 32]); xlabel('x[nm]');ylabel('E[eV]');cbar();%title('J_{\gamma}(x,E)');
ylim([turnonvolt-2.00 turnonvolt-1.00]);set(gca,'YTick',[-1.0 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(xphot,Ephot,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');
xlim([13 32]); xlabel('x[nm]');ylabel('E[eV]');cbar();%title('J_{\gamma}(x,E)');
ylim([turnonvolt-1.00 turnonvolt-0.00]);set(gca,'YTick',[-1.0 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(xR,   ER,   (abs(RphotG)+1e-25),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('$\hbar\omega$ [eV]','interpreter','latex');
xlim([13 32]); ylim([1.33 1.70]);cbar();%title('R_{\gamma,Gal}(x,E)');

num_voltB = f;

[x, phi, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev([voltB,'_phi_n_p_J_Ec_Ev']);
for ii = num_voltA+1 : num_voltB-1
    figure(ii);
    hold on;plot(x,Ec,'k','LineWidth',2);plot(x,Ev,'k','LineWidth',2);hold off;
end


% ---------------------------------------------------
% densities, LDOS, currents, Rphot at voltage C
% ---------------------------------------------------
[x, En, nE]              = negf_read_xE([voltC,'_nE']);
[x, En, pE]              = negf_read_xE([voltC,'_pE']);
[x, En, LDOS]            = negf_read_xE([voltC,'_LDOS']);
%[x, En, LDOS_VB]         = negf_read_xE([voltC,'_LDOS_VB']);
%[x, En, LDOS_k0]        = negf_read_xE([voltC,'_LDOS_k0']);
[xj, Ej, jEe]            = negf_read_xE([voltC,'_JEe']);
[xj, Ej, jEh]            = negf_read_xE([voltC,'_JEh']);
[xphon, Ephon, jEphone]  = negf_read_xE([voltC,'_Jephon']);
[xphon, Ephon, jEphonh]  = negf_read_xE([voltC,'_Jhphon']);
[xphot, Ephot, jEphote]  = negf_read_xE([voltC,'_Jephot']);
[xphot, Ephot, jEphoth]  = negf_read_xE([voltC,'_Jhphot']);
[xR, ER, RphotG]         = negf_read_xE([voltC,'_RphotGAL']);
%[xR, ER, RphotGVB]       = negf_read_xE([voltC,'_RphotGAL_VB']);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(nE+pE)+1e-15),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('n(x,E)+p(x,E)');
ylim([-0.7 1.3]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(nE+pE)+1e-15),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar();%title('n(x,E)+p(x,E)');
ylim([-0.3 1.8]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(x,    En,   (abs(LDOS)+1e-6),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar(); %title('LDOS(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xj,Ej,(abs(jEe)+abs(jEh)+1e-10),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);cbar(); %title('J(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xj,Ej,log(abs(jEe)+abs(jEh)+1e-6),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]);%title('log J(x,E)');
%f=f+1; figure(f); hold off;colormap(cm);
%contourf(xphon,Ephon,(jEphone+jEphonh+1e-10),100,'EdgeColor','none');
%xlabel('x[nm]');ylabel('E[eV]');set(gca,'YTick',[-0.5 0 0.5 1 1.5]); %title('J_{phon,e}(x,E)');
f=f+1; figure(f); hold off;colormap(cm);
contourf(xphot,Ephot,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');
xlim([13 32]); xlabel('x[nm]');ylabel('E[eV]');cbar();%title('J_{\gamma}(x,E)');
ylim([turnonvolt-1.95 turnonvolt-0.95]); set(gca,'YTick',[-1.0 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(xphot,Ephot,(abs(jEphote)+abs(jEphoth)+1e-15),100,'EdgeColor','none');
xlim([13 32]); xlabel('x[nm]');ylabel('E[eV]');cbar();%title('J_{\gamma}(x,E)');
ylim([turnonvolt-0.95 turnonvolt+0.05]);set(gca,'YTick',[-1.0 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
f=f+1; figure(f); hold off;colormap(cm);
contourf(xR,   ER,   (abs(RphotG)+1e-25),100,'EdgeColor','none');
xlabel('x[nm]');ylabel('$\hbar\omega$ [eV]','interpreter','latex');
xlim([13 32]); ylim([1.33 1.70]);cbar();%title('R_{\gamma,Gal}(x,E)');

num_voltC = f;

[x, phi, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev([voltC,'_phi_n_p_J_Ec_Ev']);
for ii = num_voltB+1 : num_voltC-1
    figure(ii);
    hold on;plot(x,Ec,'k','LineWidth',2);plot(x,Ev,'k','LineWidth',2);hold off;
end

% ------------------------------------
% formatting of all pictures so far
% ------------------------------------
labelsize = 16;
titlesize = 18;
for ii=1:num_voltC
    fig = figure(ii);
    set(fig,'Color','w','Position',[ii*50 ii*50 450 400]);
    
    fontname = 'Helvetica'; % Arial,Courier,Helvetica,Times,Bookman,lucidabright,...
    set(gca,'FontSize',14,'FontName',fontname);
    %set(get(gca,'XLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight','bold');
    %set(get(gca,'YLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight','bold');
    set(get(gca,'XLabel'),'FontName',fontname,'FontSize',labelsize);
    set(get(gca,'YLabel'),'FontName',fontname,'FontSize',labelsize);
    set(get(gca,'Title'),'FontName',fontname,'FontSize',titlesize,'FontWeight','bold');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 15 11]);
end
else % if make_negfpics
    num_voltC = 0;
end

% ---------------------------------------------
% P-V plot, I-V plot and spectra (NEGF part)
% ---------------------------------------------
f=f+1; figure(f); hold off; % will be the spectrum figure
lowercut = 5;
uppercut = 9;
%which_current_value = 2; % 2nd column (p-side, calc. from divergence)
which_current_value = 3; % 3rd column (n-side, calc. from divergence)
s = cell(1, min(length(volts),uppercut-1)-lowercut+length(aquavolts));
Power   = zeros(length(volts),1);
Current = zeros(length(volts),1);
for vv = 1:length(volts)
    fname = [filebase,sprintf('%5.3f',volts(vv)),'V.power'];
    fprintf(1,'Reading %s...',fname);
    data = load(fname);
    fprintf(1,'%d cols and %d rows.\n',size(data,2),size(data,1));
    % take last row - col 1 should have voltage, col 2 (or 3) power
    nvolts = size(data,1);
    assert(abs(data(nvolts,1))==volts(vv), 'something is not right.');
    Power(vv) = negf_power_fac * abs(data(nvolts,2));
    
    fname = [filebase,sprintf('%5.3f',volts(vv)),'V.ccurrent'];
    fprintf(1,'Reading %s...',fname);
    data = load(fname);
    fprintf(1,'%d cols and %d rows.\n',size(data,2),size(data,1));
    % take last row - col 1 should have voltage, col 2 (or 3) power
    nvolts = size(data,1);
    assert(abs(data(nvolts,1))==volts(vv), 'something is not right.');
    Current(vv) = abs(data(nvolts,which_current_value)); % no multiplication with any factor
    
    if vv<=lowercut || vv>=uppercut
        continue
    end
    %s{vv-lowercut} = sprintf('%.3g V', volts(vv));
    %s{vv-lowercut} = sprintf('%.3g W/cm^2 (NEGF)', Power(vv));
    s{vv-lowercut} = sprintf('%.3g', Power(vv));
    
    fname = [filebase,sprintf('%5.3f',volts(vv)),'V_spectrumGAL_VB'];
    [hw, P, p] = negf_spectrum(fname, negf_spectrum_fac);
    switch mod(vv-lowercut,3)
        case 1
            set(p, 'Color','black');
        case 2
            set(p, 'Color','red');
        case 0
            set(p, 'Color','blue');
        case 3
            set(p, 'Color','green');
        case 4
            set(p, 'Color','magenta');
        otherwise
            set(p, 'Color','black');
    end    
    hold on;
end
xlim([1.35 1.75]); % spectrum

% ---------------------
% AQUA simulation
% ---------------------
[aquafields, aquadata] = aqua_read(aquafile);
aqua_Vidx = find_idx_in_list(aquafields, 'voltage_Source');
aqua_Iidx = find_idx_in_list(aquafields, 'contact_current_Source');
aqua_Pidx = find_idx_in_list(aquafields, 'total_light_emission_0');

for vv = 1:length(aquavolts)
    fname = [aquatree,'_spectrum_Source',sprintf('%5.3f',aquavolts(vv)),'V.slc'];
    [energy,spectrum,p] = aqua_spectrum(fname,aqua_spectrum_fac);
    set(p, 'LineStyle','--','LineWidth',2);
    switch mod(vv,3)
        case 1
            set(p, 'Color','black');
        case 2
            set(p, 'Color','red');
        case 0
            set(p, 'Color','blue');
        case 3
            set(p, 'Color','green');
        case 4
            set(p, 'Color','magenta');
        otherwise
            set(p, 'Color','black');
    end 
    
    %s{length(volts)+vv} = sprintf('%.3g V (AQUA)', aquavolts(vv));
    % find aqua-voltage vv in .xy-file
    for ii=1:size(aquadata,1)
        if aquadata(ii,aqua_Vidx)==aquavolts(vv)
            %s{min(length(volts),uppercut-1)-lowercut+vv} = sprintf('%.3g W/cm^2 (AQUA)', aqua_power_fac*aquadata(ii,aqua_Pidx));
            s{min(length(volts),uppercut-1)-lowercut+vv} = sprintf('%.03g', aqua_power_fac*aquadata(ii,aqua_Pidx));
            break;
        end
    end
end
legend(s);
xlabel('$\hbar\omega$ [eV]','interpreter','latex');
ylim([0 1e21]);ylabel('Spectrum [s^{-1}/cm^2]');
title('');

% P-V figure
f=f+1; figure(f); hold off;

%Power ./ Current
%aquadata(:,aqua_Pidx) ./ aquadata(:,aqua_Iidx)

semilogy(volts,Power,'k-s','LineWidth',2);hold on;
semilogy(aquadata(:,aqua_Vidx), aqua_power_fac*aquadata(:,aqua_Pidx),'k--','LineWidth',2);
semilogy(volts,Current,'b-s','LineWidth',2);
semilogy(aquadata(:,aqua_Vidx), aqua_power_fac*aquadata(:,aqua_Iidx),'b--','LineWidth',2);
hold off;
ylabel({'Output power [W/cm^2]';'Current [A/cm^2]'});
xlim([1.09 1.71]);ylim([1e-5 1e3]);
set(gca,'XTick',[1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9]);
set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
leg = legend('NEGF P(V)','AQUA P(V)','NEGF I(V)','AQUA I(V)');set(leg,'Location','NorthWest');

xlabel('Voltage [V]');


% I-V figure
%f=f+1; fig=figure(f); hold off;
%semilogy(volts,Current,'k-x','LineWidth',2);hold on;
%semilogy(aquadata(:,aqua_Vidx), aqua_power_fac*aquadata(:,aqua_Iidx),'b--x','LineWidth',2);hold off;
%xlabel('Voltage [V]');ylabel('Current [A/cm^2]');
%set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
%leg = legend('NEGF','AQUA');set(leg,'Location','NorthWest');
%xlim([1.0 1.6]);%ylim([1e-5 1e3]);

% formatting
labelsize = 16;
titlesize = 18;
fontweight = 'normal'; % 'bold' 
for ii=num_voltC+1:f
    fig = figure(ii);
    set(fig,'Color','w','Position',[ii*50 ii*50 450 400]);    
    fontname = 'Helvetica'; % Arial,Courier,Helvetica,Times,Bookman,lucidabright,...
    set(gca,'FontSize',14,'FontName',fontname);
    set(get(gca,'XLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight',fontweight);
    set(get(gca,'YLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight',fontweight);
    set(get(gca,'Title'),'FontName',fontname,'FontSize',titlesize,'FontWeight','bold');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 15 11]);
end

frm = 'epsc'; % png, epsc, ...
ext = 'eps';  % png, eps, ...
res = '-r200';

f = 0;
if make_negfpics
f=f+1;figure(f); filename = [voltA,'_npE1.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_npE2.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_LDOS.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_JE.',ext];     print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_JElog.',ext];  print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
%f=f+1;figure(f); %filename = [voltA,'_JEphon.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_JEphot1.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_JEphot2.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltA,'_Rphot.',ext];  print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);

f=f+1;figure(f); filename = [voltB,'_npE1.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_npE2.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_LDOS.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_JE.',ext];     print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_JElog.',ext];  print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
%f=f+1;figure(f); %filename = [voltB,'_JEphon.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_JEphot1.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_JEphot2.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltB,'_Rphot.',ext];  print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);

f=f+1;figure(f); filename = [voltC,'_npE1.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_npE2.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_LDOS.',ext];   print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_JE.',ext];     print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_JElog.',ext];  print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
%f=f+1;figure(f); %filename = [voltC,'_JEphon.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_JEphot1.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_JEphot2.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [voltC,'_Rphot.',ext];  print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
end

f=f+1;figure(f); filename = [filebase,'spectra.',ext];print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
f=f+1;figure(f); filename = [filebase,'PV.',ext];     print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);
%f=f+1;figure(f); filename = [filebase,'IV.',ext];     print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);

end % iwce_prepare_pics
