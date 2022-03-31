function d = minmod(x,y)
d = 0.5*(sign(x)+sign(y)).*min(abs(x),abs(y));
end