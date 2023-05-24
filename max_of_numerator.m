function res = max_of_numerator(number_of_points, R, p2, DIGITS)

digits(DIGITS);

z_inner = vpa(zeros(number_of_points, 1));
z_outer = vpa(zeros(number_of_points, 1));

temp_inner = vpa(zeros(number_of_points, 1));
temp_outer = vpa(zeros(number_of_points, 1));

h_phi = vpa(vpa(2)*vpa(pi)/vpa(vpa(number_of_points) - vpa(1)));
for i = 1 : 1 : number_of_points
    z_inner(i,1) = vpa(vpa(1)/vpa(R)*vpa(exp(vpa(1j)*...
        vpa(vpa(i) - vpa(1))*vpa(h_phi))));
    temp_inner(i,1) = vpa(abs(vpa(...
        numerator(z_inner(i,1), p2, DIGITS))));
    z_outer(i,1) = vpa(vpa(R)*vpa(exp(vpa(1j)*...
        vpa(vpa(i) - vpa(1))*vpa(h_phi))));
    temp_outer(i,1) = vpa(abs(vpa(...
        numerator(z_outer(i,1), p2, DIGITS))));
end

temp_inner_max = max(temp_inner);
temp_outer_max = max(temp_outer);
res = max(temp_inner_max, temp_outer_max);
disp(res);
end

