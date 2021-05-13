%Strong stability preserving runge kutta - 3rd order
%https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html
function yf = SSPRK3(F,y0,dt)
y1 = y0 + dt*F(y0);
y2 = 3/4*y0 + 1/4*(y1 + dt*F(y1));
yf = 1/3*y0 + 2/3*(y2 + dt*F(y2));
end