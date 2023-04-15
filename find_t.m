function res = find_t(DIG, dN, d2N, d4N)
digits(DIG);
%--------------------------------------------------------------------------
flag = 4;
%--------------------------------------------------------------------------
switch flag
    case 1
        t_prev = vpa(root1(DIG, dN, d2N, d4N));
    case 2
        t_prev = vpa(root2(DIG, dN, d2N, d4N));
    case 3
        t_prev = vpa(root3(DIG, dN, d2N, d4N));
    case 4
        t_prev = 1;%vpa(root4(DIG, dN, d2N, d4N));
end
%--------------------------------------------------------------------------
N = 10;
for k = 1 : 1 : N
    e = vpa(E(DIG, t_prev, dN, d2N, d4N));
    e_d = vpa(E_d(DIG, t_prev, dN, d2N, d4N));
    tmp = abs(e / e_d)
    t_next = vpa(t_prev - e / e_d);
    t_prev = vpa(t_next);
end
%--------------------------------------------------------------------------
res = vpa(t_prev);
%--------------------------------------------------------------------------
disp(vpa(E(DIG, t_prev, dN, d2N, d4N)));
end

