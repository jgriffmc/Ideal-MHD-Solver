% 5th order point-stencil reconstruction of cell boundary with limiters
% Returns reconstructon to the left cell boundary and to the right cell
% boundary:  |L   cell   R|

% Suresh, A. and H. Huynh. “Accurate Monotonicity-Preserving Schemes with 
% Runge-Kutta Time Stepping.” Journal of Computational Physics 136 (1997): 
% 83-99.

function [UL,UR] = recon_MP5(U)
alpha = 4;
eps = 1e-10; %tolerance

%calculate circshifts
Um2 = circshift(U,2);
Um1 = circshift(U,1);
Up1 = circshift(U,-1);
Up2 = circshift(U,-2);

%Reconstruction to right boundary
UR_est = 1/60*(2*Um2 - 13*Um1 + 47*U + 27*Up1 - 3*Up2); %SP Eqn 2.1
U_MP = U + minmod(Up1-U,alpha*(U - Um1)); %SP Eqn 2.12
d = Um1-2*U+Up1; %SP Eqn 2.19
dm1 = circshift(d,1); %SP Eqn 2.19
dp1 = circshift(d,-1); %SP Eqn 2.19
dM4p1 = minmod4(4*d-dp1,4*dp1-d,d,dp1); %SP Eqn 2.27
dM4m1 = minmod4(4*d-dm1,4*dm1-d,d,dm1);
U_UL = U + alpha*(U - Um1); %SP Eqn 2.8
U_AV = 1/2*(U + Up1); %SP Eqn 2.16
U_MD = U_AV - 1/2*dM4p1; %SP Eqn 2.28
U_LC = U + 1/2*(U-Um1) + 4/3*dM4m1; %SP Eqn 2.29
U_min = max(min(min(U,Up1),U_MD),min(min(U,U_UL),U_LC)); %SP Eqn 2.24a
U_max = min(max(max(U,Up1),U_MD),max(max(U,U_UL),U_LC)); %SP Eqn 2.24b
eps_check = ((UR_est-U).*(UR_est-U_MP)<=eps); %SP Eqn 2.30

%If eps_check is true, UR = UR_est
%Else, UR = median(UR_est,Umin,Umax)
UR = UR_est; %set UR = UR_est by default
UR_med = median_vec(UR_est,U_min,U_max);
UR(~eps_check) = UR_med(~eps_check); %set masked elements = UR_median


%Reconstruction to left boundary
UL_est = 1/60*(-3*Um2 + 27*Um1 + 47*U - 13*Up1 + 2*Up2); %SP Eqn 2.1
U_MP = U + minmod(Um1-U,alpha*(U - Up1)); %SP Eqn 2.12
U_UL = U + alpha*(U - Up1); %SP Eqn 2.8
U_AV = 1/2*(U + Um1); %SP Eqn 2.16
U_MD = U_AV - 1/2*dM4m1; %SP Eqn 2.28
U_LC = U + 1/2*(U-Up1) + 4/3*dM4p1; %SP Eqn 2.29
U_min = max(min(min(U,Um1),U_MD),min(min(U,U_UL),U_LC)); %SP Eqn 2.24a
U_max = min(max(max(U,Um1),U_MD),max(max(U,U_UL),U_LC)); %SP Eqn 2.24b
eps_check = ((UL_est-U).*(UL_est-U_MP)<=eps); %SP Eqn 2.30

%If eps_check is true, UL = UL_est
%Else, UL = median(UL_est,Umin,Umax)
UL = UL_est; %set UL = UL_est by default
UL_med = median_vec(UL_est,U_min,U_max);
UL(~eps_check) = UL_med(~eps_check); %set masked elements = UL_median
end