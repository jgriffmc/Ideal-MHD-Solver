%Four argument minmod
function s = minmod4(a,b,c,d)
sa = sign(a);
sb = sign(b);
sc = sign(c);
sd = sign(d);
s = 1/8*(sa + sb).*abs((sa + sc).*(sa + sd)).*min(min(abs(a),abs(b)),min(abs(c),abs(d)));
end