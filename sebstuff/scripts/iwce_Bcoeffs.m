
filebase = '/home/steiger/src/release/amd64/NEGF/results/iwce_GOOD/iwce_Source-';
volts = [1.200,1.300,1.400];
%filebase = '/home/steiger/src/release/amd64/NEGF/results/pn/photons/pn_Source-';
%volts = [0.800,0.900,1.000,1.100,1.200,1.250,1.300,1.350,1.400,1.450,1.500];

s = cell(1, length(volts));
figure(3); hold off;
for vv = 1:length(volts)
    fname = [filebase,sprintf('%5.3f',volts(vv)),'V'];
    [x, B, p] = negf_Bcoeff(fname);
    switch mod(vv,5)
        case 0
            set(p, 'Color','black');
        case 1
            set(p, 'Color','red');
        case 2
            set(p, 'Color','blue');
        case 3
            set(p, 'Color','green');
        case 4
            set(p, 'Color','magenta');
        otherwise
            set(p, 'Color','black');
    end

    s{vv} = sprintf('%.3g V', volts(vv));
    hold on;
end
legend(s);

