% Matlab script to save .png-images of ANGEL results
%
% call: matlab -nodesktop -nosplash -r save_pngs
% make sure you have the file 'result.files'
% -nodesktop messes up command line in coates, grrrrr...

% filebase = '../results/nanohub/nanohub_Source0.000V';

fid = fopen('result.files');
assert(~feof(fid), 'check result.files');
filedir  = fgetl(fid);
assert(~feof(fid), 'check result.files');
filebase = fgetl(fid);
fclose(fid);

% densities
[x, En, nE] = negf_read_xE([filedir,filebase,'_nE']);
[x, En, pE] = negf_read_xE([filedir,filebase,'_pE']);
f(1) = figure('Visible','off');
contourf(x,En,(abs(nE)+abs(pE)+1e-15),100,'LineColor','none');
title('n(x,E)+p(x,E)');

% local DOS
[x, En, LDOS] = negf_read_xE([filedir,filebase,'_LDOS']);
f(2) = figure('Visible','off');
contourf(x,En,(abs(LDOS)+1e-6),100,'LineColor','none');
title('LDOS(x,E)');

% current
[x, En, jEe] = negf_read_xE([filedir,filebase,'_JEe']);
[x, En, jEh] = negf_read_xE([filedir,filebase,'_JEh']);
% note: current plot is wrong because 1st row & col were assumed as axes
f(3) = figure('Visible','off');
contourf(x,En,jEe+jEh,100,'LineColor','none');colorbar;
title('J_{i,i+1}(x,E)');

[x, phi, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev([filedir,filebase,'_phi_n_p_J_Ec_Ev']);

% formatting
labelsize = 16;
titlesize = 18;
for ii=1:3
	figure(ii); set(gcf, 'Visible', 'off');
    
    % add band edges to plots
	hold(gca);
    plot(x,Ec,'k','LineWidth',2);
	plot(x,Ev,'k','LineWidth',2);
	hold(gca);
    
    xlabel('x[nm]');
	ylabel('E[eV]');
    set(gcf,'Color','w');
    fontname = 'Helvetica'; % Arial,Courier,Helvetica,Times,Bookman,lucidabright,...
    set(gca,'FontSize',14,'FontName',fontname);
    set(get(gca,'XLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight','bold');
    set(get(gca,'YLabel'),'FontName',fontname,'FontSize',labelsize,'FontWeight','bold');
    set(get(gca,'Title'),'FontName',fontname,'FontSize',titlesize,'FontWeight','bold');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 15 11]);
end

% save to file
format = 'png';
res = '-r200';
figure(1); set(gcf, 'Visible', 'off'); filename = [filebase,'_npE.',format];  print(res,['-d',format],filename);fprintf(1,'figure 1 saved to %s\n',filename);
figure(2); set(gcf, 'Visible', 'off'); filename = [filebase,'_LDOS.',format]; print(res,['-d',format],filename);fprintf(1,'figure 2 saved to %s\n',filename);
figure(3); set(gcf, 'Visible', 'off'); filename = [filebase,'_JE.',format];   print(res,['-d',format],filename);fprintf(1,'figure 3 saved to %s\n',filename);

exit

