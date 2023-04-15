function res = E_d(DIG, t, dN, d2N, d4N)
digits(DIG);
xi = vpa(4*d4N+2*dN-6*d2N);
eta = vpa(8*d4N+6*dN-16*d2N);
res = vpa(d2N*(t-1)+d2N*(1+t)+(t-1)*t*xi*(2*t+(xi-eta)/xi)+... 
 (t-1)*xi*(1+t^2+(t*(xi-eta))/xi)+t*xi*(1+t^2+(t*(xi-eta))/xi));
end

