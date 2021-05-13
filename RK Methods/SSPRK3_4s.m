% Four stage SSP-RK3 (allows CFL<=2)
% https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html
function yf = SSPRK3_4s(F,y0,dt)
y1 = 1/2*y0 + 1/2*(y0 + dt*F(y0));
y2 = 1/2*y1 + 1/2*(y1 + dt*F(y1));
y3 = 2/3*y0 + 1/6*y2 + 1/6*(y2 + dt*F(y2));
yf = 1/2*y3 + 1/2*(y3 + dt*F(y3));
end
