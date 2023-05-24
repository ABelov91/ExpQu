function res = find_t(dN, d2N, d4N, DIGITS)
%--------------------------------------------------------------------------
digits(DIGITS);
%--------------------------------------------------------------------------
t_prev_1 = vpa(root1(DIGITS, dN, d2N, d4N));
t_prev_2 = vpa(root2(DIGITS, dN, d2N, d4N));
t_prev_3 = vpa(root3(DIGITS, dN, d2N, d4N));
t_prev_4 = vpa(root4(DIGITS, dN, d2N, d4N));
%--------------------------------------------------------------------------
t_prev_vec = vpa([vpa(t_prev_1), vpa(t_prev_2), vpa(t_prev_3),...
    vpa(t_prev_4)]);
E_d_vec = vpa(E_d(t_prev_vec, dN, d2N, d4N, DIGITS));
E_d_vec_norm = vpa(vpa(E_d_vec).*conj(vpa(E_d_vec)));
%--------------------------------------------------------------------------
[~, flag] = min(vpa(E_d_vec_norm));
%--------------------------------------------------------------------------
t_prev = vpa(t_prev_vec(1, flag));
%--------------------------------------------------------------------------
res = vpa(t_prev);
disp(flag);
%--------------------------------------------------------------------------
end

