
filebase = '/home/steiger/src/release/amd64/NEGF/results/pn/';
file1 = [filebase,'pn_Source-0.800V'];
%file1 = [filebase,'POP_LuisierSR_broad0.001_Imk0.5nm/pn_Source-0.800V'];
file2 = [filebase,   'POP_EasySR_broad0.001_Imk0.5nm/pn_Source-0.800V'];

% densities

figure(4);
[x1, En1, nE1] = negf_read_xE([file1,'_nE']);
[x1, En1, pE1] = negf_read_xE([file1,'_pE']);
subplot(1,3,1);
contourf(x1,En1,(abs(nE1)+abs(pE1)+1e-15),100,'EdgeColor','none');title('Luisier S^R n(x,E)+p(x,E)');colorbar;
%contourf(x1,En1,(abs(nE1)+1e-15),100,'EdgeColor','none');title('Luisier S^R n(x,E)+p(x,E)');
[x, phi1, n1, p1, J1, Ec1, Ev1] = negf_read_phi_n_p_J_Ec_Ev([file1,'_phi_n_p_J_Ec_Ev']);
hold on;plot(x,Ec1,'k','LineWidth',2);plot(x,Ev1,'k','LineWidth',2);hold off;


[x2, En2, nE2] = negf_read_xE([file2,'_nE']);
[x2, En2, pE2] = negf_read_xE([file2,'_pE']);
subplot(1,3,2);
contourf(x2,En2,(abs(nE2)+abs(pE2)+1e-15),100,'EdgeColor','none');title('Easy S^R n(x,E)+p(x,E)');colorbar;
%contourf(x2,En2,(abs(nE2)+1e-15),100,'EdgeColor','none');title('Easy S^R n(x,E)+p(x,E)');
[x, phi2, n2, p2, J2, Ec2, Ev2] = negf_read_phi_n_p_J_Ec_Ev([file2,'_phi_n_p_J_Ec_Ev']);
hold on;plot(x,Ec2,'k','LineWidth',2);plot(x,Ev2,'k','LineWidth',2);hold off;

subplot(1,3,3);
contourf(x2,En2,(abs(nE1-nE2)+abs(pE1-pE2)+1e-15),100,'EdgeColor','none');title('Difference');colorbar;
hold on;
plot(x,Ec1,'k','LineWidth',2);plot(x,Ev1,'k','LineWidth',2);
plot(x,Ec2,'k','LineWidth',2);plot(x,Ev2,'k','LineWidth',2);
hold off;