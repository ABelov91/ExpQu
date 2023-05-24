function res = exact2(L, p1, p2, DIGITS)
digits(DIGITS);
res = vpa(-(((1 - (1 - p1^(L))^(-1))/cosh((p1 - p2)^(-1))^(2))/...
    (p1 - p2)^(2)) - (p1^(-1 + L)*L*tanh((p1 - p2)^(-1)))/(1 - p1^L)^2);
end

