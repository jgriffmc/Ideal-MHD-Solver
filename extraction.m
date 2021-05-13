function [rho, mo, E, B, v, rho_e, p] = extraction(ynext)
global gamma
rho = ynext(:,1);
mo = ynext(:,2);
E = ynext(:,3);
B = ynext(:,4);
v = mo./rho;
rho_e = E - 0.5*mo.*v - 0.5*B.^2;
p = rho_e*(gamma-1);
end