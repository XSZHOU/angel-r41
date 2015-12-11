function calc_mobility(Volt,Curr,L,offset,num_points,plot_yes_or_no,col)

linewidth  = 1.5;
markersize = 8;

ec   = 1.602e-19;
dens = 2e18;
m0   = 9.1e-31;
me   = 0.067*m0;

markers = ['^','>','<','v'];
s = cell(1, size(Volt,2)-1);

if (plot_yes_or_no)
    for ii=2:size(Volt,2)
        r = Volt(offset+1:offset+num_points,ii) ./ Curr(offset+1:offset+num_points,ii);
        plot(L,r,[col 'x'],'LineWidth',linewidth,'MarkerSize',markersize,'Marker',markers(ii-1));
        s{ii-1} = sprintf('%5.2fV', Volt(1,ii));
    end
    legend(s);
end
for ii=2:size(Volt,2)
    r = Volt(offset+1:offset+num_points,ii) ./ Curr(offset+1:offset+num_points,ii);
    f = polyfit(L*1e-7,r',1); % find linear coeff. y=a+bx
    fprintf(1,'V=%4g[V]: r = %9.3e[Ohm cm^2], rho = %9.3e[Ohm cm]\n',Volt(1,ii),f(2),f(1));
    if (plot_yes_or_no)
        plot(L,f(2)+L*1e-7*f(1),col,'LineWidth',linewidth);
    end
    mu   = 1/f(1)/ec/dens;
    rate = ec / (mu*1e-4 * me);
    fprintf(1,'    mobility      from sigma=1/rho=n*e*mu: %9.3e[cm2/Vs]\n', mu);
    fprintf(1,'    scattering rate from 1/tau = e/(mu*m): %9.3e[1/s]\n', rate);
end

return
