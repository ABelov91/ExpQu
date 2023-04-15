function res = extrap2(DIG, dN, d2N, d4N, N, a)
digits(DIG);
t = vpa(find_t(DIG, dN, d2N, d4N));
%--------------------------------------------------------------------------
disp(vpa(N))
disp(vpa(1 / (1 - (a)^N) - t));
%--------------------------------------------------------------------------
A = vpa((2*t-1)^(2)/(t*(t-1))*(-2*dN*t^(2)/(2*t-1)^(2)+d2N*t/(t-1)));
BN = vpa((2*t-1)^(2)/(t*(t-1))*(dN/(2*t-1)-d2N/(t-1)));
%--------------------------------------------------------------------------
res = vpa((-1) * (t - 1)^8 / (t^8 - (t - 1)^8) *...
    (A + BN * 8 * t^8 / (t^8 - (t - 1)^8)));
end
