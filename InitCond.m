% This function defines the initial conditions for mass density, velocity,
% pressure, and magnetic field (normalized to 1/sqrt(u0))

function [rho, v, p, B] = InitCond(n)
global x
%Initialize to zero
rho = zeros(size(x)); %mass density
v = zeros(size(x)); %velocity
p = zeros(size(x)); %pressure
B = zeros(size(x)); %magnetic field (normalized such that B = sqrt(u0)*B)
switch n     
    case 0 %Sod shock tube from Toro
        a = 0.5;
        rho(x<a) = 1;
        rho(x>=a) = 0.125;
        p(x<a) = 1;
        p(x>=a) = 0.1;
        
    case 1 %Shocktube from Hammett slides
        a = 0.5;
        rho(x<a) = 3;
        rho(x>=a) = 1;
        p(x<a) = 3;
        p(x>=a) = 1;

    case 2 %Acoustic wave from local pressure increase
        a = 0.5;
        rho = ones(size(rho));
        p = 1 + 0.1*exp(-((x-a)/0.1).^2);
        
    case 3 %Acoustic wave with constant B
        a = 0.5;
        rho = ones(size(rho));
        p = 1 + 2*exp(-((x-a)/0.1).^2);
        B = 1*ones(size(B));
        
    case 4 %Sod shock tube with constant B
        a = 0.5;
        rho(x<a) = 1;
        rho(x>=a) = 0.125;
        p(x<a) = 1;
        p(x>=a) = 0.1;
        B = 1*ones(size(B));
        
    case 5 %Acoustic wave with constant B (SI units)
        a = 6e-6;
        rho = 1.2*ones(size(rho)); %kg/m3
        p = 101325*(1 + 5*exp(-((x-a)/0.5e-6).^2)); %Pa
        B = 0.2 * 1/sqrt(4*pi*1e-7)*ones(size(B)); %T/sqrt(u0)
        
    case 6 %Sod Shock Tube (SI units)
        a = 6e-6;
        rho(x<a) = 0.1; %kg/m3
        rho(x>=a) = 2;
        rho(x>=a+4e-6) = 0.1;
        p = 101325*ones(size(p)); %Pa
        p(x>=a & x<=a+4e-6) = 1.6e8;
        B = 0 * 1/sqrt(4*pi*1e-7)*ones(size(B)); %T/sqrt(u0)

    otherwise
        error('Incorrect choice of initial conditions')
end

end
