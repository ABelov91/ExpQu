function res = numerator(z, p2, DIGITS)
digits(DIGITS);
res = vpa(tanh(vpa(1) / (vpa(z) - vpa(p2))));
end

