

aqua_file = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-13_fermi_mh0.2/qw-iwce_Source1.600V_1d.xy';
negf_file = '/home/steiger/src/release/amd64/NEGF/results/iwce_GOOD/iwce_Source-1.600V_phi_n_p_J_Ec_Ev';
outdir    = '/home/steiger/src/release/amd64/aqua/results/qw-iwce_t1e-13_fermi_mh0.2/';

linewidth = 2;
labelsize = 14;
fontname = 'Helvetica';

xlim   = [0 45];
ylim1  = [-1.0 2.0];
ytick1 = [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
ylim2  = [0 1.1e25];
ytick2 = [0 5e24 1e25];

f = 0;

% ====================================
% NEGF band edges and densities
% ====================================

[xNEGF, phiNEGF, nNEGF, pNEGF, JNEGF, EcNEGF, EvNEGF] = negf_read_phi_n_p_J_Ec_Ev(negf_file);
% Ec, Ev already have potential in them!

f=f+1; fig=figure(f); clf; set(fig,'Color','w','Units','centimeters','Position',[f f xpaper ypaper]);
ax0 = gca; set(ax0,'Box','on','Color','white','XTick',[],'YTick',[]);

% first axes for left y-axis
ax1 = axes('Position',get(ax0,'Position'));
pp1 = plot(ax1, xNEGF, EcNEGF, 'k', xNEGF, EvNEGF,'k'); %hold on;
set(ax1,'Box','off','Color','none','YAxisLocation','left',...
    'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim1,'YTick',ytick1);
set(pp1,'LineWidth',linewidth,'LineStyle','-');

% second axes for right y-axis assuming common x-axis controlled by ax1
ax2 = axes('Position', get (ax0, 'Position'));
pp2 = plot(ax2, xNEGF, 1e27*nNEGF, 'b', xNEGF, 1e27*pNEGF, 'r'); 
set(ax2,'Box','off','Color','none','XTick',[],'YAxisLocation','right', ...
    'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim2,'YTick',ytick2); 
set(pp2,'LineWidth',linewidth,'LineStyle','-');


% ====================================
% AQUA band edges and densities
% ====================================


[aquafields, aquadata] = aqua_read(aqua_file);
aqua_xidx   = find_idx_in_list(aquafields, 'ycoord');
aqua_phiidx = find_idx_in_list(aquafields, 'potential_3d');
aqua_Ecidx  = find_idx_in_list(aquafields, 'conduction_bandedge');
aqua_Evidx  = find_idx_in_list(aquafields, 'valence_bandedge');
aqua_n3idx  = find_idx_in_list(aquafields, 'edensity_3d');
aqua_n2idx  = find_idx_in_list(aquafields, 'edensity_interpolation_well0');
aqua_p3idx  = find_idx_in_list(aquafields, 'hdensity_3d');
aqua_p2idx  = find_idx_in_list(aquafields, 'hdensity_interpolation_well0');

  x = aquadata(aqua_xidx,:);
phi = aquadata(aqua_phiidx,:);
 Ec = aquadata(aqua_Ecidx,:);
 Ev = aquadata(aqua_Evidx,:);
 n3 = aquadata(aqua_n3idx,:);
 n2 = aquadata(aqua_n2idx,:);
 p3 = aquadata(aqua_p3idx,:);
 p2 = aquadata(aqua_p2idx,:);

shift = EcNEGF(1) - (Ec(1)-phi(1));
fprintf(1,'Shift: %.3f\n', shift);

shift_right =  EcNEGF(length(EcNEGF)) - (Ec(length(Ec))-phi(length(phi)));
fprintf(1,'Difference in potential step: %.3f\n', shift - shift_right);

fprintf(1,'NEGF edensity at right boundary: %.2e[cm-3]\n',1e21*nNEGF(length(nNEGF)));

new_figure = 0;
if new_figure
    fig = figure(2); clf; set(fig,'Color','w');
    ax0 = gca; set(ax0,'Box','on','Color','white','XTick',[],'YTick',[]);
end

% first axes for left y-axis
ax3 = axes('Position',get(ax0,'Position'));
pp3 = plot(ax3, x, shift+Ec-phi,'k', x, shift+Ev-phi,'k'); %hold on;
set(ax3,'Box','off','Color','none','YAxisLocation','left',...
    'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim1,'YTick',ytick1);
set(pp3,'LineWidth',linewidth,'LineStyle','--');

% second axes for right y-axis assuming common x-axis controlled by ax1
ax4 = axes('Position', get (ax0, 'Position'));
pp4 = plot(ax4, x,n3+n2,'b',x,p3+p2,'r'); 
set(ax4,'Box','off','Color','none','XTick',[],'YAxisLocation','right', ...
    'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim2,'YTick',ytick2); 
set(pp4,'LineWidth',linewidth,'LineStyle','--');

% =========================
% More Formatting
% =========================
set(get(ax0,'XLabel'),'String',{'';'x [nm]'},'FontName',fontname,'FontSize',labelsize);
set(get(ax3,'YLabel'),'String','E [eV]'     ,'FontName',fontname,'FontSize',labelsize);
set(get(ax4,'YLabel'),'String','n,p [cm^-3]','FontName',fontname,'FontSize',labelsize);

% =======================
% Save to file
% =======================

frm = 'epsc';
ext = 'eps';
dpi        = 300; % dots per inch
dots_percm = 2.54*dpi;
res        = '-r300';
xpaper     = 9;         % cm
ypaper     = 6.5;

filename = [outdir,'npEcEv_AQUANEGF.',ext]; 
print(res,['-d',frm],filename);fprintf(1,'figure %d saved to %s\n',f,filename);

return







% old code with plotyy

fig = figure(2); clf; set(fig,'Color','w');
[AXn,H1n,H2n] = plotyy(10*x, shift + Ec-phi, 10*x, n3+n2); hold on;
[AXp,H1p,H2p] = plotyy(10*x, shift + Ev-phi, 10*x, p3+p2);

set(get(AXn(1),'Ylabel'),'String','E [eV]');
set(get(AXn(2),'Ylabel'),'String','n, p [m^{-3}]');
set(AXn(1),'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim1,'YTick',ytick1);
set(AXn(2),'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim2,'YTick',ytick2);
set(H1n,'LineStyle','-','LineWidth',linewidth,'Color','k');
set(H2n,'LineStyle','-','LineWidth',linewidth,'Color','b');

set(AXp(1),'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim1,'YTick',ytick1);
set(AXp(2),'YColor','k','FontSize',labelsize,'FontName',fontname,'XLim',xlim,'YLim',ylim2,'YTick',ytick2);
set(H1p,'LineStyle','-','LineWidth',linewidth,'Color','k');
set(H2p,'LineStyle','-','LineWidth',linewidth,'Color','r');






