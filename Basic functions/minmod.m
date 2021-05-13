%2 argument minmod
function s = minmod(sL,sR)
s = 0.5*(sign(sL)+sign(sR)).*min(abs(sL),abs(sR));
end