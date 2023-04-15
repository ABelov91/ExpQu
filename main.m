clear; clc;

DIG = 200;
digits(DIG);

% a = vpa(1 / pi);
a = vpa(1 * (pi + 1j * exp(1)) / sqrt(pi^2 + exp(2)) / pi);

q = 2;
b = vpa(1e3 / pi);

% N > 4
N = 9;
r = 2;
K = vpa(r.^((1 : 1 : N).'));

z = vpa(zeros(N, r^(N)));
hz = vpa(zeros(N, r^(N)));

for k = 1 : 1 : N
    L = r^(k);
    h_phi = vpa(2 * pi / L);
    z(k, 1 : L) = vpa(exp(1j * (0 : h_phi : (2 * pi - h_phi))));
    hz(k, 1 : L) = vpa(exp(1j * (0 : h_phi : (2 * pi - h_phi))) *...
        h_phi * 1j);
end

f = vpa(g2(z, a, vpa(q), b));

for k = 1 : 1 : N
    f(k, (r^(k) + 1) : (r^(N))) = vpa(zeros(1, (r^(N) - r^(k))));
end

GN = vpa(sum((f .* hz), 2) / (2 * pi * 1j));
G = vpa(GN(N, 1));
% 1 -> v; ... N-1 -> v; N -> x
DN = vpa(G .* ones((N - 1), 1) - GN(1 : (N - 1), 1));
% 1 -> x; 2 -> x; 3 ->  x; 4 -> v; ... N -> v
DN_extrap = vpa(zeros((N - 4), 1));
GN_extrap = vpa(zeros((N - 4), 1));
for k = 1 : 1 : (N - 4)
    DN_extrap(k, 1) = vpa(extrap2(DIG, DN(k, 1), DN(k + 1, 1),...
        DN(k + 2, 1), vpa(r^k), vpa(a)));
    GN_extrap(k, 1) = vpa(GN(k + 3) + DN_extrap(k, 1));
end
DN_after_extrap = vpa(G .* ones((N - 4), 1) - GN_extrap);

dN = vpa(log10(DN.'' .* DN) / 2);
dN_extrap = vpa(log10(DN_extrap.'' .* DN_extrap) / 2);
dN_after_extrap = vpa(log10(DN_after_extrap.'' .* DN_after_extrap) / 2);

plot(K(1 : (N - 1), 1), dN, 'o-black');
hold on;
plot(K(4 : (N - 1), 1), dN_extrap, '-r');
hold on;
plot(K(4 : (N - 1), 1), dN_after_extrap, '.--b');
xlabel('N');
ylabel('log_{10}R');
title('Graph of convergence');

