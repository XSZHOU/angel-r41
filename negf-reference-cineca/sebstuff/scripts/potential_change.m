function [x, phi, dphi] = potential_change(filebase)

clear phi dphi

stepmin = 10;
stepmax = 50;

for step = stepmin:stepmax
    file = [filebase,'_step',int2str(step),'_phi_n_p_J_Ec_Ev'];
    [x, phi_tmp, n, p, J, Ec, Ev] = negf_read_phi_n_p_J_Ec_Ev(file);
    phi(step-stepmin+1, :) = phi_tmp;
end
size(phi,1)
size(phi,2)
dphi = zeros(size(phi,1)-1,size(phi,2));
for step = 1:(size(phi,1)-1)
    dphi(step,:) = phi(step+1,:) - phi(step,:);
end