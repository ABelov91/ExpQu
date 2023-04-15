function res = E(DIG, t, dN, d2N, d4N)
digits(DIG);
xi = vpa(4*d4N+2*dN-6*d2N);
eta = vpa(8*d4N+6*dN-16*d2N);
res = vpa(xi*t*(t-1)*(t^(2)+(xi-eta)*t/xi+1)+d2N*(t-1)*(t+1)+d4N);
end

