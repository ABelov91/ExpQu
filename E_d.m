function res = E_d(t, dN, d2N, d4N, DIGITS)
%--------------------------------------------------------------------------
digits(DIGITS);
%--------------------------------------------------------------------------
xi = vpa(vpa(4)*vpa(d4N)+vpa(2)*vpa(dN)-vpa(6)*vpa(d2N));
eta = vpa(vpa(8)*vpa(d4N)+vpa(6)*vpa(dN)-vpa(16)*vpa(d2N));
%--------------------------------------------------------------------------
res = vpa(vpa(d2N).*(vpa(t)-vpa(1))+vpa(d2N).*(vpa(1)+vpa(t))+...
    (vpa(t)-vpa(1)).*vpa(t).*vpa(xi).*(vpa(2).*vpa(t)+(vpa(xi)-...
    vpa(eta))/vpa(xi))+(vpa(t)-vpa(1)).*vpa(xi).*(vpa(1)+vpa(t).*...
    vpa(t)+(vpa(t).*(vpa(xi)-vpa(eta)))./vpa(xi))+vpa(t).*vpa(xi).*...
    (vpa(1)+vpa(t).*vpa(t)+(vpa(t).*(vpa(xi)-vpa(eta)))./vpa(xi)));
%--------------------------------------------------------------------------
end

