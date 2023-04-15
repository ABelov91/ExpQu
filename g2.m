function res = g2(z, a, q, b)
temp = 1;
for i = 1 : 1 : q
temp = temp .* (z - a);
end
res = tanh(1 ./ (z - b)) ./ temp;
end

