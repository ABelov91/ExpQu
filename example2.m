clear; clc;

DIGITS = 150;
digits(DIGITS);

p1 = vpa(vpa(1) / vpa(pi));
q = 2;
p2 = vpa(vpa(1000) / vpa(pi));

number_of_points = 1000;
R = vpa(2.8);
M = vpa(max_of_numerator(number_of_points, R, p2, DIGITS));

N = 7;
r = 2;
K = vpa(zeros(N, 1));
for i = 1 : 1 : N
    K(i, 1) = vpa(vpa(r)^vpa(i));
end

z = vpa(zeros(N, r^(N)));
hz = vpa(zeros(N, r^(N)));

for i = 1 : 1 : N
    L = r^(i);
    h_phi = vpa(2 * pi / L);
    for j = 1 : 1 : L
        z(i, j) = vpa(exp(vpa(1j) * vpa(h_phi) * vpa(j - 1)));
        hz(i, j) = vpa(vpa(z(i, j)) * vpa(h_phi) * vpa(1j));
    end
end

f = vpa(zeros(N, r^(N)));

for i = 1 : 1 : N
    L = r^(i);
    for j = 1 : 1 : L
        f(i, j) = vpa(g2(vpa(z(i, j)), q, vpa(p1), vpa(p2),...
            DIGITS));
    end
end

GN = vpa(zeros(N, 1));

for i = 1 : 1 : N
    L = r^(i);
    temp = vpa(0);
    for j = 1 : 1 : L
        temp = vpa(vpa(temp) + vpa(f(i, j)) * vpa(hz(i, j)));
    end
    GN(i, 1) = vpa(temp) / vpa(vpa(2) * vpa(pi) * vpa(1j));
end

% G_approx = vpa(GN(N, 1));
G_exact = vpa(vpa(-1)/(vpa(cosh(vpa((vpa(p1) -...
    vpa(p2))^vpa(-1)))))^vpa(2)/(p1 - p2)^(2));
G = G_exact;

DN = vpa(zeros(N - 1, 1));

for i = 1 : 1 : (N - 1)
    DN(i, 1) = vpa(vpa(G) - vpa(GN(i, 1)));
end

DN_extrap = vpa(zeros(N - 4, 1));
GN_extrap = vpa(zeros(N - 4, 1));

for i = 1 : 1 : (N - 4)
    DN_extrap(i, 1) = vpa(extrap2(vpa(DN(i, 1)), vpa(DN(i + 1, 1)),...
        vpa(DN(i + 2, 1)), DIGITS));
    GN_extrap(i, 1) = vpa(vpa(GN(i + 3)) + vpa(DN_extrap(i, 1)));
end

DN_after_extrap = vpa(zeros(N - 4, 1));

for i = 1 : 1 : (N - 4)
    DN_after_extrap(i, 1) = vpa(vpa(G) - vpa(GN_extrap(i, 1)));
end

DN_extrap_exact = vpa(zeros(N, 1));
GN_extrap_exact = vpa(zeros(N, 1));

for i = 1 : 1 : N
    DN_extrap_exact(i, 1) = vpa(exact2(vpa(K(i, 1)),...
        vpa(p1), vpa(p2), DIGITS));
    GN_extrap_exact(i, 1) = vpa(vpa(GN(i, 1)) +...
        vpa(DN_extrap_exact(i, 1)));
end

DN_after_extrap_exact = vpa(zeros(N, 1));

for i = 1 : 1 : N
    DN_after_extrap_exact(i, 1) = vpa(vpa(G_exact) -...
        vpa(GN_extrap_exact(i, 1)));
end

DN_TW = vpa(zeros(N, 1));

for i = 1 : 1 : N
    DN_TW(i, 1) = vpa(TW_estimate(M, R, K(i, 1), DIGITS));
end
        
dN = vpa(log10(DN.'' .* DN) / 2);

dN_extrap = vpa(log10(DN_extrap.'' .* DN_extrap) / 2);

dN_after_extrap = vpa(log10(DN_after_extrap.'' .* DN_after_extrap) / 2); 

dN_extrap_exact = vpa(log10(DN_extrap_exact.'' .* DN_extrap_exact) / 2);

dN_after_extrap_exact = vpa(log10(DN_after_extrap_exact.'' .*...
    DN_after_extrap_exact) / 2);

dN_TW = vpa(log10(DN_TW.'' .* DN_TW) / 2);

% graphics

plot(K(1 : (N - 1), 1), dN, 'Color', 'black', 'Marker', 'o',...
    'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white',...
    'MarkerSize', 9, 'LineStyle', '-', 'LineWidth', 1);
hold on;
plot(K(1 : (N - 1), 1), dN_after_extrap_exact(1 : (N - 1), 1),...
    'Color', 'black', 'Marker', 's', 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', 'black', 'MarkerSize', 9,...
    'LineStyle', '-', 'LineWidth', 1);
hold on;
plot(K(4 : (N - 1), 1), dN_after_extrap, 'Color', 'black', 'Marker',...
    's', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white',...
    'MarkerSize', 9, 'LineStyle', '-', 'LineWidth', 1);
hold on;
plot(K(1 : N, 1), dN_TW, 'Color', 'black', 'Marker', 'o',...
    'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black',...
    'MarkerSize', 4, 'LineStyle', '--', 'LineWidth', 1);

xlabel('N');
ylabel('log_{10}|error|');

xlim([0 70]);
ylim([-60 0]);
legend('fact error', 'a priori extrapolation',...
    'a posteriori extrapolation', 'T.-W. estimate');

grid off;
set(gca, 'Box', 'off');
FigureHandle = gcf;
set(findall(FigureHandle,'type','text'),'FontSize',14,...
    'FontName','Times New Roman');
saveas(gcf, 'fig_02.svg');
saveas(gcf, 'fig_02.eps');
saveas(gcf, 'fig_02.png');
