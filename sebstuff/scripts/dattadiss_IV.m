filetree = '/home/steiger/src/release/amd64/NEGF/results/diss/';

figure(1); clf;
V = zeros(15,5);
I = zeros(15,5);
Lac  = [60 72 84 96 108 120];
Lopt = [60 72 84 96 108 120];
ball  =  6; % max. index for ballistic sims (10...acous) 
acous = 15; % max. index for acoustic sims (10...acous) - important for polyfit
opt   = 25; % index from which acoustic scattering starts - important for polyfit
[V(1,:),I(1,:),p1]=negf_IV([filetree 'ballistic/dattadiss60nm1.2dx/dattadiss60nm1.2dx_Source0.200V.ccurrent']);hold on;
[V(2,:),I(2,:),p2]=negf_IV([filetree 'ballistic/dattadiss72nm1.2dx/dattadiss72nm1.2dx_Source0.200V.ccurrent']);
[V(3,:),I(3,:),p3]=negf_IV([filetree 'ballistic/dattadiss84nm1.2dx/dattadiss84nm1.2dx_Source0.200V.ccurrent']);
[V(4,:),I(4,:),p4]=negf_IV([filetree 'ballistic/dattadiss96nm1.2dx/dattadiss96nm1.2dx_Source0.200V.ccurrent']);
[V(5,:),I(5,:),p5]=negf_IV([filetree 'ballistic/dattadiss108nm1.2dx/dattadiss108nm1.2dx_Source0.200V.ccurrent']);
[V(6,:),I(6,:),p6]=negf_IV([filetree 'ballistic/dattadiss120nm1.2dx/dattadiss120nm1.2dx_Source0.200V.ccurrent']);


[V(10,:),I(10,:),p10]=negf_IV([filetree 'AC/dattadiss60nm1.2dx/dattadiss60nm1.2dx_Source0.200V.ccurrent']);
[V(11,:),I(11,:),p11]=negf_IV([filetree 'AC/dattadiss72nm1.2dx/dattadiss72nm1.2dx_Source0.200V.ccurrent']);
[V(12,:),I(12,:),p12]=negf_IV([filetree 'AC/dattadiss84nm1.2dx/dattadiss84nm1.2dx_Source0.200V.ccurrent']);
[V(13,:),I(13,:),p13]=negf_IV([filetree 'AC/dattadiss96nm1.2dx/dattadiss96nm1.2dx_Source0.200V.ccurrent']);
[V(14,:),I(14,:),p14]=negf_IV([filetree 'AC/dattadiss108nm1.2dx/dattadiss108nm1.2dx_Source0.200V.ccurrent']);
[V(15,:),I(15,:),p15]=negf_IV([filetree 'AC/dattadiss120nm1.2dx/dattadiss120nm1.2dx_Source0.200V.ccurrent']);

[V(20,:),I(20,:),p20]=negf_IV([filetree 'POP/dattadiss60nm1.2dx/dattadiss60nm1.2dx_Source0.200V.ccurrent']);
[V(21,:),I(21,:),p21]=negf_IV([filetree 'POP/dattadiss72nm1.2dx/dattadiss72nm1.2dx_Source0.200V.ccurrent']);
[V(22,:),I(22,:),p22]=negf_IV([filetree 'POP/dattadiss84nm1.2dx/dattadiss84nm1.2dx_Source0.200V.ccurrent']);
[V(23,:),I(23,:),p23]=negf_IV([filetree 'POP/dattadiss96nm1.2dx/dattadiss96nm1.2dx_Source0.200V.ccurrent']);
[V(24,:),I(24,:),p24]=negf_IV([filetree 'POP/dattadiss108nm1.2dx/dattadiss108nm1.2dx_Source0.200V.ccurrent']);
[V(25,:),I(25,:),p25]=negf_IV([filetree 'POP/dattadiss120nm1.2dx/dattadiss120nm1.2dx_Source0.200V.ccurrent']);

set(gcf,'Color','w');
% ballistic
set(p1,'Color','black');
set(p2,'Color','black');
set(p3,'Color','black');
set(p4,'Color','black');
set(p5,'Color','black');
set(p6,'Color','black');
msize = 10;
% AC
set(p10,'Color','blue','Marker','x','MarkerSize',msize);
set(p11,'Color','blue','Marker','o','MarkerSize',msize);
set(p12,'Color','blue','Marker','s','MarkerSize',msize);
set(p13,'Color','blue','Marker','d','MarkerSize',msize);
set(p14,'Color','blue','Marker','p','MarkerSize',msize);
set(p15,'Color','blue','Marker','h','MarkerSize',msize);
% POP
set(p20,'Color','red','Marker','x','MarkerSize',msize);
set(p21,'Color','red','Marker','o','MarkerSize',msize);
set(p22,'Color','red','Marker','s','MarkerSize',msize);
set(p23,'Color','red','Marker','d','MarkerSize',msize);
set(p24,'Color','red','Marker','p','MarkerSize',msize);
set(p25,'Color','red','Marker','h','MarkerSize',msize);
%set(p26,'Color','cyan'   ,'Marker','x','MarkerSize',20);
hold off;


figure(2); clf; set(gcf,'Color','w'); hold on;
xlabel('Voltage [V]');
ylabel('Resistance [Ohm cm2]');
for ii=1:size(V,1)
    r = V(ii,:)./I(ii,:);
    p = plot(V(ii,:),r,'Color',[rand rand rand]);
end
legend('60nm ball.','60nm','72nm','84nm','96nm','108nm','120nm');
hold off;

ec   = 1.602e-19;
dens = 2e18;
figure(3); clf; set(gcf,'Color','w'); hold on;
xlabel('Resistor length [nm]');
ylabel('Resistance V/I [Ohm cm^2]');

fprintf(1,'ACOUSTIC PHONON SCATTERING:\n');
for ii=2:size(V,2) % for all voltages
    r = V(10:acous,ii)./I(10:acous,ii);
    plot(Lac,r,'bx','MarkerSize',10);
    f = polyfit(Lac*1e-7,r',1); % find linear coeff. y=a+bx
    fprintf(1,'V=%4g[V]: r_{ball} = %9.3e[Ohm cm^2], rho_{scatt} = %9.3e[Ohm cm]\n',V(10,ii),f(2),f(1));
    plot(Lac,f(2)+Lac*1e-7*f(1),'b');
    %fprintf(1,'    sigma=1/rho=%9.3e [C/(cmVs)]\n',1/f(1));
    fprintf(1,'    mu from sigma=1/rho=n*e*mu: %9.3e[cm2/Vs]\n', 1/f(1)/ec/dens);
end

fprintf(1,'OPTICAL PHONON SCATTERING:\n');
for ii=2:size(V,2) % for all voltages
    r = V(20:opt,ii)./I(20:opt,ii);
    plot(Lopt,r,'rx','MarkerSize',10);
    f = polyfit(Lopt*1e-7,r',1); % find linear coeff. y=a+bx
    fprintf(1,'V=%4g[V]: r_{ball} = %9.3e[Ohm cm^2], rho_{scatt} = %9.3e[Ohm cm]\n',V(20,ii),f(2),f(1));
    plot(Lopt,f(2)+Lopt*1e-7*f(1),'r');
    %fprintf(1,'    sigma=1/rho=%9.3e [C/(cmVs)]\n',1/f(1));
    fprintf(1,'    mu from sigma=1/rho=n*e*mu: %9.3e[cm2/Vs]\n', 1/f(1)/ec/dens);
end
legend('0.05V','0.1V','0.15V','0.2V');
hold off;



