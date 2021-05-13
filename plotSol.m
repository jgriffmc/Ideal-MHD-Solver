%Plots the solution of density, velocity, pressure, magnetic field/internal
%energy

function plotSol(y,tc,plotB)
global x

[rho, ~, ~, B, v, rho_e, p] = extraction(y); %extract quantities

subplot(2,2,1)
hold on
plot(x,rho,'.-')
xlabel('x')
title('Density')

subplot(2,2,2)
hold on
plot(x,v,'.-')
xlabel('x')
title('Velocity')

subplot(2,2,3)
hold on
plot(x,p,'.-')
xlabel('x')
title('Pressure')

subplot(2,2,4)
hold on
if(plotB) %plot B or e depending on selection
    plot(x,B*sqrt(4*pi*1e-7),'.-')
    xlabel('x')
    title('Magnetic field')
else
    plot(x,rho_e./rho,'.-')
    xlabel('x')
    title('Internal Energy')
end

sgtitle(sprintf('t = %.2f',tc));

end