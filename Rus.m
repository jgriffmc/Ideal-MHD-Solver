%Calculates the Rusanov Flux in cartesian geometry
% Inputs:
%   F(U)  = flux at the left and right cell boundaries
%   U = State variable
%   vmax = maximum solution speed (v + ss + vA)
%   isP = treats pressure term in cylindrical coordinates as cartesian flux

function dUdt = Rus(FL,FR,UL,UR,vmax,isP)
global dx dim x

%calculate global diffusion coefficient
DR = max(vmax)*dx/2; %right boundary
DL = max(vmax)*dx/2; %left boundary

%calculate numerical fluxes
fR = 1/2*(FR+circshift(FL,-1)) - DR/dx.*(circshift(UL,-1)-UR); %right boundary flux
fL = 1/2*(FL+circshift(FR,1)) - DL/dx.*(UL-circshift(UR,1)); %left boundary flux

if(dim==0) %cartesian
    %Flux BCs
    fR(end) = fL(end); %flux in = flux out
    fL(1) = fR(1);
    
    dUdt = -(fR-fL)/dx;
    
elseif(dim==1) %cylindrical
    xL = x - dx/2;
    xR = x + dx/2;
    
    fR = fR(3:end-2);
    fL = fL(3:end-2);
    
    %Flux BCs
    fR(end) = fL(end); %flux in = flux out
    
    if(isP) %pressure term
        dUdt = -(fR-fL)/dx;
    else
        dUdt = -1./x.*(fR.*xR-fL.*xL)/dx;
    end
end