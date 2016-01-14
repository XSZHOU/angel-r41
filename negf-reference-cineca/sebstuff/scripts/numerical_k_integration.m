% compare numerical integration, exact result
% 2D parabolic band structure hbar^2 k^2 / 2m

SIhbar = 1.06E-34;       % Planck constant (h/2pi)
SIm0   = 9.1E-31;        % Bare electron mass
SIec   = 1.60E-19;       % Elementary charge

% units: nm - s - eV - ec
hbar = SIhbar / SIec;
m    = SIm0  *1e-18 / SIec; % 1 kg = 1e-18/SIec eV s^2 nm^-2
EF   = 0.0;
kT   = 0.025;

nE = 100;
E = linspace(-0.3,0.6,nE);

nk   = 50;
k_BZ = 2*pi/5e-10 * 1e-9;
kmax = 0.5*k_BZ;
%kmax = 3.0;

% analytical result
F0 = m*kT/(2*pi*hbar^2) * log(1+exp((EF-E)/kT));

% numerical result - 2D k-integration
k = linspace(0,kmax,nk);
numeric = zeros(1,nE);
for i=1:nk
    if i==1
        dk = (k(i+1)-k(i)) / 2;
    else if i==nk
        dk = (k(i)-k(i-1)) / 2;
        else
        dk = (k(i+1)-k(i-1)) / 2;
        end
    end
    dk = 2*pi * k(i) * dk / (4*pi^2);
    Ek = hbar^2*k(i)^2 / (2*m); 
    tmp = 1.0 ./ (1+exp((E + Ek - EF)/kT));
    numeric(1,:)  = numeric(1,:) + tmp *dk;
end
figure(3);
subplot(2,1,1);
semilogy(E,F0,'k-o',E,numeric,'r',E,abs(F0-numeric),'b',E,abs(F0-numeric)./F0,'b--');
legend('exact', 'numeric (\int dk^2)','difference (abs)', 'difference (rel)');
xlabel('E [eV]'); ylabel('F_0(E) [nm^{-2}]');
titlestring = sprintf('Numerical (nk=%.0f, k_{max}=%.1f nm^{-1}) vs. exact computation',nk,kmax);
title(titlestring);



Efix = 0.3;
% dependence on kmax
Nkmax = 30;
kmax_scaling = zeros(Nkmax,1);
scale = 1.1;
kmax_scaling(1) = k_BZ / 20;
for ii=2:Nkmax
    kmax_scaling(ii) = scale * kmax_scaling(ii-1);
end

% dependence on nk
Nnk = 70;
nk_inc = 10;
nk_vals = linspace(10,Nnk*nk_inc,Nnk);

F0 = m*kT/(2*pi*hbar^2) * log(1+exp((EF-Efix)/kT)) * ones(Nkmax,Nnk);

numeric = zeros(Nkmax,Nnk);
for ii=1:Nkmax
    for jj=1:Nnk
        kmax = kmax_scaling(ii);
        nk   = nk_vals(jj);
        
        k = linspace(0,kmax,nk);
        for i=1:nk
            if i==1
                dk = (k(i+1)-k(i)) / 2;
            else if i==nk
                dk = (k(i)-k(i-1)) / 2;
                else
                dk = (k(i+1)-k(i-1)) / 2;
                end
            end
            dk = 2*pi * k(i) * dk / (4*pi^2);
            Ek = hbar^2*k(i)^2 / (2*m); 
            tmp = 1.0 / (1+exp((Efix + Ek - EF)/kT));
            numeric(ii,jj)  = numeric(ii,jj) + tmp *dk;
        end
    end
end
subplot(2,1,2);

h = surf(kmax_scaling,nk_vals,(abs(F0-numeric)./F0)','EdgeColor','k');
titlestring = sprintf('Numerical vs. exact computation: relative difference, E=%.3g',Efix);
title(titlestring);
xlabel('k_{max} [nm^{-1}]');
ylabel('nk'); 
%legend('difference (rel)');
set(get(h,'Parent'),'ZScale','log'); 
%set(get(h,'Parent'),'XScale','log'); 
set(get(h,'Parent'),'YScale','log'); 
xlim([min(kmax_scaling) max(kmax_scaling)]);
