% This program solves the 1d ideal MHD equations in either cartesian or
% cylindrical coordinate systems
%
% Jesse Griff-McMahon
% APC527: Numerical Algorithms for Scientific Computing
% May 13, 2021

clear
% close all
global dx x gamma dim

%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%
%DOMAIN
xl = 0; %leftmost cell boundary
xr = 1; %rightmost cell boundary
Nx = 200+1; %number of gridpoints (including endpoints)
CFL = 0.1; %CFL number

%TIME
t0 = 0; %start time
tf = 1; %final time in simulation
dtout = 0.05; %intermediate outputs times (integer divisible of tf) (set to 0 to suppress intermediate outputs)

%ADIABATIC INDEX
gamma = 7/5;

%DIMENSION
dim = 1; %choose 0 for cartesian, 1 for cylindrical

%INITIAL CONDITIONS
ICsel = 4; %select initial condition
% 0: Sod shock tube from "Riemann Solvers ..." by Toro
% 1: Shocktube from Hammett slides
% 2: Acoustic wave from local pressure increase
% 3: Acoustic wave with constant B
% 4: Sod shock tube with constant B
% 5: Acoustic wave with constant B (SI units)
% 6: Sod Shock Tube (SI units)

%PLOT INPUTS
plotIC = 1; %choose 1 to plot the initial condition, 0 to omit plot
plotB = 1; %Choose 1 to plot the magnetic field instead of the internal energy

%RUNGE KUTTA SOLVER
RKsel = 1; %select Runge Kutta solver
% 0: RK1 (Forward Euler)
% 1: SSP RK3
% 2: SSP RK3 (4 step)
% 3: SSP RK4

%%%%%%%%%%%%%%%% GRID %%%%%%%%%%%%%%%%
dx = (xr - xl)/Nx; %spatial grid spacing
x = (xl+dx/2:dx:xr-dx/2)'; %spatial grid
tout = t0:dtout:tf; %Intermediate output time grid
if(~isempty(tout)) 
    tout(end) = inf; %prevent last output from displaying as intermediate output
end
t_count = 2; %index for intermediate output

%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%
[rho, v, p, B] = InitCond(ICsel);

%%%%%%%%%%%%%%%% TRANSFORM TO CONSERVED QUANTITIES %%%%%%%%%%%%%%%%
mo = rho.*v; %momentum density
rho_e = p/(gamma-1); %internal energy density
E = rho_e + 0.5*rho.*v.^2 + 0.5*B.^2; %total energy density

y = [rho,mo,E,B]; %initial solution vector

%%%%%%%%%%%%%%%% PLOT INITIAL CONDITION %%%%%%%%%%%%%%%%
if(plotIC)
    figure(1)
    clf
    plotSol(y,t0,plotB);
end

%%%%%%%%%%%%%%%% SELECT RK SOLVER %%%%%%%%%%%%%%%%
switch RKsel
    case 0
        RKsolver = @RK1;
    case 1
        RKsolver = @SSPRK3;
    case 2
        RKsolver = @SSPRK3_4s;
    case 3
        RKsolver = @SSPRK4;
    otherwise
        error('Incorrect RK selection')
end

%%%%%%%%%%%%%%%% MAIN COMPUTATION LOOP %%%%%%%%%%%%%%%%
%guess initial dt
vmax = max(abs(v) + sqrt(gamma*p./rho) + abs(B)./sqrt(rho)); %maximum e'value of the system
dt = min(CFL*dx/vmax,(tf-t0)/20); %if guess is bigger than 1/20 of time interval, choose 1/20 of time interval
tc = t0; %initialize current time

while(tc<tf)
    dt_toobig = 1; %initialize, assume dt is too big
    
    while(dt_toobig==1)
        ynext = RKsolver(@(y) DE_solver(y),y,dt); %calculate next step of solution
        [rho, mo, E, B, v, rho_e, p] = extraction(ynext); %extract quantities
        
        %calculate solution velocity to check if dt is small enough
        vmax = max(abs(v) + sqrt(gamma*p./rho) + abs(B)./sqrt(rho)); %include Alfven velocity
        
        %check if dt satisfies CFL condition
        if(dt>CFL*dx/vmax) % dt is too large
            dt = min((tf-t0)/20,CFL*dx/vmax); %calculate new dt w/o advancing solution
        else
            dt_toobig = 0; %dt is good enough to proceed
            y = ynext; %advance solution
            dt = min((tf-t0)/20,CFL*dx/vmax); %calculate new dt
        end
    end
    
    tc = tc + dt; %advance current time
    
    if(dtout>0 && tc>=tout(t_count)) %plot intermediate solution
        plotSol(y,tc,plotB);
        t_count = t_count + 1;
        pause %press a key to continue
    end
end

%%%%%%%%%%%%%%%% PLOT FINAL CONDITION %%%%%%%%%%%%%%%%
plotSol(y,tc,plotB);
