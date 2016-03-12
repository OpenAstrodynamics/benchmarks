function el = elements(r, v, mu)
rm = norm(r);
vm = norm(v);
h = cross(r, v);
hm = norm(h);
k = [0.0, 0.0, 1.0];
n = cross(k, h);
nm = norm(n);
xi = vm^2/2 - mu/rm;
ec = ((vm^2 - mu/rm)*r - v*dot(r,v))/mu;
ecc = norm(ec);
if ecc ~= 1.0
    sma = -mu/(2*xi);
else
    sma = hm^2/mu;
end
inc = acos(h(3)/hm);
node = acos(n(1)/nm);
peri = acos(dot(n,ec)/(ecc*nm));
ano = acos(dot(ec,r)/(ecc*rm));
if n(2) < 0
    node = 2*pi - node;
end
if ec(3) < 0
    peri = 2*pi - peri;
end
if dot(r,v) < 0
    ano = 2*pi - ano;
end
el = [sma, ecc, inc, node, peri, ano];
