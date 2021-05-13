%median used for vectors x,y,z

function m = median_vec(x,y,z)
m = x + minmod(y-x,z-x);
end
