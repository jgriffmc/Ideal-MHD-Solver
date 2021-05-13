%Strong stability preserving runge kutta - 4th order
function yf = SSPRK4(F,y0,dt)
y1 = y0;
y2 = y0 + dt/2*F(y1);
y3 = y0 + dt/2*F(y2);
y4 = y0 + dt*F(y3);
yf = y0 + dt/6*(F(y1) + 2*F(y2) + 2*F(y3) + F(y4));
end
