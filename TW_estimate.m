function res = TW_estimate(M, R, N, DIGITS)
digits(DIGITS);
res = vpa(vpa(2) * vpa(M) * vpa((vpa(R)^vpa(N)-vpa(1))^vpa(-1)));
end

