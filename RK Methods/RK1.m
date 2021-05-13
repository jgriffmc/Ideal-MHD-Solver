%Forward euler method

function yf = RK1(F,y0,dt)
yf = y0 + dt*F(y0);
end
