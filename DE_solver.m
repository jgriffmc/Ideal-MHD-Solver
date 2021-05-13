%Solves the Ideal MHD Equations
%Returns the time derivative of the conserved quantities

function dydt = DE_solver(y0)
global gamma dim

[rho, mo, E, B, v, ~, p] = extraction(y0); %extract quantities

a = sqrt(gamma*p./rho); %speed of sound
vA = B./sqrt(rho); %Alfven speed

v_plus_ss = abs(v)+a+abs(vA); %used in calculating the diffusion coefficient (velocity + sound speed + alfven speed)

if(dim==0) %cartesian
    %reconstruct
    [rhoL,rhoR] = recon_MP5(rho); %left and right cell boundaries:  |L   R| 
    [moL,moR] = recon_MP5(mo);
    [EL,ER] = recon_MP5(E);
    [BL,BR] = recon_MP5(B);
    
    %calculate derived quantities
    vL = moL./rhoL;
    vR = moR./rhoR;
    pL = (EL - 0.5*moL.*vL - 0.5*BL.^2)*(gamma-1);
    pR = (ER - 0.5*moR.*vR - 0.5*BR.^2)*(gamma-1);

    %Call rusanov solver
    drho_dt = Rus(moL,moR,rhoL,rhoR,v_plus_ss,0);
    dmo_dt = Rus(vL.*moL+pL+BL.^2/2,vR.*moR+pR+BR.^2/2,moL,moR,v_plus_ss,0);
    dE_dt = Rus(vL.*(EL+pL+BL.^2/2),vR.*(ER+pR+BR.^2/2),EL,ER,v_plus_ss,0);
    dB_dt = Rus(vL.*BL,vR.*BR,BL,BR,v_plus_ss,0);

elseif(dim==1) %cylindrical
    %reconstruct with ghost cells
    rho = [rho(2);rho(1);rho;rho(end);rho(end)];
    mo = [-mo(2);-mo(1);mo;mo(end);mo(end)];
    E = [E(2);E(1);E;E(end);E(end)]; 
    B = [B(2);B(1);B;B(end);B(end)]; 
    v_plus_ss = [v_plus_ss(2);v_plus_ss(1);v_plus_ss;v_plus_ss(end);v_plus_ss(end)];
    
    [rhoL,rhoR] = recon_MP5(rho);
    [moL,moR] = recon_MP5(mo);
    [EL,ER] = recon_MP5(E);
    [BL,BR] = recon_MP5(B);
    
    %calculate derived quantities
    vL = moL./rhoL;
    vR = moR./rhoR;
    pL = (EL - 0.5*moL.*vL - 0.5*BL.^2)*(gamma-1);
    pR = (ER - 0.5*moR.*vR - 0.5*BR.^2)*(gamma-1);
    
    drho_dt = Rus(moL,moR,rhoL,rhoR,v_plus_ss,0);
    dmo_dt = Rus(pL+BL.^2/2,pR+BR.^2/2,moL,moR,v_plus_ss,1) + Rus(vL.*moL,vR.*moR,moL,moR,v_plus_ss,0);
    dE_dt = Rus(vL.*(EL+pL+BL.^2/2),vR.*(ER+pR+BR.^2/2),EL,ER,v_plus_ss,0);
    dB_dt = Rus(vL.*BL,vR.*BR,BL,BR,v_plus_ss,0);
end

dydt = [drho_dt,dmo_dt,dE_dt,dB_dt]; %final solution vector
end