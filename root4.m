function res = root4(DIG, a, b, c)
digits(DIG);
res =...
vpa(-(-vpa(8)*vpa(b) + vpa(4)*vpa(c) + ...
vpa(3)*vpa(a))/(vpa(4)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) + ...
sqrt(-(vpa(15)*vpa(b) - vpa(8)*vpa(c) - ...
vpa(6)*vpa(a))/(vpa(2)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) + ...
(-vpa(8)*vpa(b) + vpa(4)*vpa(c) + ...
vpa(3)*vpa(a))^vpa(2)/(vpa(4)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - ...
vpa(a))^vpa(2)) - (-vpa(15)*vpa(b) + vpa(8)*vpa(c) + ...
vpa(6)*vpa(a))/(vpa(6)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) - ...
(vpa(3)*vpa(b) - ...
vpa(4)*vpa(c))^vpa(2)/(vpa(3)*vpa(2)^(vpa(2)/vpa(3))*(vpa(3)*vpa(b) - ...
vpa(2)*vpa(c) - vpa(a))*(-vpa(54)*vpa(b)^vpa(3) + ...
vpa(288)*vpa(b)*vpa(c)^vpa(2) - vpa(128)*vpa(c)^vpa(3) - ...
vpa(216)*vpa(b)*vpa(c)*vpa(a) + vpa(108)*vpa(c)*vpa(a)^vpa(2) + ...
sqrt(vpa(23328)*vpa(b)^vpa(5)*vpa(c) - ...
vpa(108864)*vpa(b)^vpa(4)*vpa(c)^vpa(2) + ...
vpa(152064)*vpa(b)^vpa(3)*vpa(c)^vpa(3) - ...
vpa(55296)*vpa(b)^vpa(2)*vpa(c)^vpa(4) + ...
vpa(23328)*vpa(b)^vpa(4)*vpa(c)*vpa(a) - ...
vpa(124416)*vpa(b)^vpa(2)*vpa(c)^vpa(3)*vpa(a) + ...
vpa(55296)*vpa(b)*vpa(c)^vpa(4)*vpa(a) - ...
vpa(11664)*vpa(b)^vpa(3)*vpa(c)*vpa(a)^vpa(2) + ...
vpa(46656)*vpa(b)^vpa(2)*vpa(c)^vpa(2)*vpa(a)^vpa(2) + ...
vpa(62208)*vpa(b)*vpa(c)^vpa(3)*vpa(a)^vpa(2) - ...
vpa(27648)*vpa(c)^vpa(4)*vpa(a)^vpa(2) - ...
vpa(46656)*vpa(b)*vpa(c)^vpa(2)*vpa(a)^vpa(3) + ...
vpa(11664)*vpa(c)^vpa(2)*vpa(a)^vpa(4)))^(vpa(1)/vpa(3))) - ...
(-vpa(54)*vpa(b)^vpa(3) + vpa(288)*vpa(b)*vpa(c)^vpa(2) - ...
vpa(128)*vpa(c)^vpa(3) - vpa(216)*vpa(b)*vpa(c)*vpa(a) + ...
vpa(108)*vpa(c)*vpa(a)^vpa(2) + sqrt(vpa(23328)*vpa(b)^vpa(5)*vpa(c) ...
- vpa(108864)*vpa(b)^vpa(4)*vpa(c)^vpa(2) + ...
vpa(152064)*vpa(b)^vpa(3)*vpa(c)^vpa(3) - ...
vpa(55296)*vpa(b)^vpa(2)*vpa(c)^vpa(4) + ...
vpa(23328)*vpa(b)^vpa(4)*vpa(c)*vpa(a) - ...
vpa(124416)*vpa(b)^vpa(2)*vpa(c)^vpa(3)*vpa(a) + ...
vpa(55296)*vpa(b)*vpa(c)^vpa(4)*vpa(a) - ...
vpa(11664)*vpa(b)^vpa(3)*vpa(c)*vpa(a)^vpa(2) + ...
vpa(46656)*vpa(b)^vpa(2)*vpa(c)^vpa(2)*vpa(a)^vpa(2) + ...
vpa(62208)*vpa(b)*vpa(c)^vpa(3)*vpa(a)^vpa(2) - ...
vpa(27648)*vpa(c)^vpa(4)*vpa(a)^vpa(2) - ...
vpa(46656)*vpa(b)*vpa(c)^vpa(2)*vpa(a)^vpa(3) + ...
vpa(11664)*vpa(c)^vpa(2)*vpa(a)^vpa(4)))^(vpa(1)/vpa(3))/(vpa(6)*vpa(...
2)^(vpa(1)/vpa(3))*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))))/vpa(2) ...
+ sqrt(-(vpa(15)*vpa(b) - vpa(8)*vpa(c) - ...
vpa(6)*vpa(a))/(vpa(2)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) + ...
(-vpa(8)*vpa(b) + vpa(4)*vpa(c) + ...
vpa(3)*vpa(a))^vpa(2)/(vpa(2)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - ...
vpa(a))^vpa(2)) + (-vpa(15)*vpa(b) + vpa(8)*vpa(c) + ...
vpa(6)*vpa(a))/(vpa(6)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) + ...
(vpa(3)*vpa(b) - ...
vpa(4)*vpa(c))^vpa(2)/(vpa(3)*vpa(2)^(vpa(2)/vpa(3))*(vpa(3)*vpa(b) - ...
vpa(2)*vpa(c) - vpa(a))*(-vpa(54)*vpa(b)^vpa(3) + ...
vpa(288)*vpa(b)*vpa(c)^vpa(2) - vpa(128)*vpa(c)^vpa(3) - ...
vpa(216)*vpa(b)*vpa(c)*vpa(a) + vpa(108)*vpa(c)*vpa(a)^vpa(2) + ...
sqrt(vpa(23328)*vpa(b)^vpa(5)*vpa(c) - ...
vpa(108864)*vpa(b)^vpa(4)*vpa(c)^vpa(2) + ...
vpa(152064)*vpa(b)^vpa(3)*vpa(c)^vpa(3) - ...
vpa(55296)*vpa(b)^vpa(2)*vpa(c)^vpa(4) + ...
vpa(23328)*vpa(b)^vpa(4)*vpa(c)*vpa(a) - ...
vpa(124416)*vpa(b)^vpa(2)*vpa(c)^vpa(3)*vpa(a) + ...
vpa(55296)*vpa(b)*vpa(c)^vpa(4)*vpa(a) - ...
vpa(11664)*vpa(b)^vpa(3)*vpa(c)*vpa(a)^vpa(2) + ...
vpa(46656)*vpa(b)^vpa(2)*vpa(c)^vpa(2)*vpa(a)^vpa(2) + ...
vpa(62208)*vpa(b)*vpa(c)^vpa(3)*vpa(a)^vpa(2) - ...
vpa(27648)*vpa(c)^vpa(4)*vpa(a)^vpa(2) - ...
vpa(46656)*vpa(b)*vpa(c)^vpa(2)*vpa(a)^vpa(3) + ...
vpa(11664)*vpa(c)^vpa(2)*vpa(a)^vpa(4)))^(vpa(1)/vpa(3))) + ...
(-vpa(54)*vpa(b)^vpa(3) + vpa(288)*vpa(b)*vpa(c)^vpa(2) - ...
vpa(128)*vpa(c)^vpa(3) - vpa(216)*vpa(b)*vpa(c)*vpa(a) + ...
vpa(108)*vpa(c)*vpa(a)^vpa(2) + sqrt(vpa(23328)*vpa(b)^vpa(5)*vpa(c) ...
- vpa(108864)*vpa(b)^vpa(4)*vpa(c)^vpa(2) + ...
vpa(152064)*vpa(b)^vpa(3)*vpa(c)^vpa(3) - ...
vpa(55296)*vpa(b)^vpa(2)*vpa(c)^vpa(4) + ...
vpa(23328)*vpa(b)^vpa(4)*vpa(c)*vpa(a) - ...
vpa(124416)*vpa(b)^vpa(2)*vpa(c)^vpa(3)*vpa(a) + ...
vpa(55296)*vpa(b)*vpa(c)^vpa(4)*vpa(a) - ...
vpa(11664)*vpa(b)^vpa(3)*vpa(c)*vpa(a)^vpa(2) + ...
vpa(46656)*vpa(b)^vpa(2)*vpa(c)^vpa(2)*vpa(a)^vpa(2) + ...
vpa(62208)*vpa(b)*vpa(c)^vpa(3)*vpa(a)^vpa(2) - ...
vpa(27648)*vpa(c)^vpa(4)*vpa(a)^vpa(2) - ...
vpa(46656)*vpa(b)*vpa(c)^vpa(2)*vpa(a)^vpa(3) + ...
vpa(11664)*vpa(c)^vpa(2)*vpa(a)^vpa(4)))^(vpa(1)/vpa(3))/(vpa(6)*vpa(...
2)^(vpa(1)/vpa(3))*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) + ...
(vpa(8) + (vpa(2)*(vpa(15)*vpa(b) - vpa(8)*vpa(c) - ...
vpa(6)*vpa(a))*(-vpa(8)*vpa(b) + vpa(4)*vpa(c) + ...
vpa(3)*vpa(a)))/(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))^vpa(2) - ...
(-vpa(8)*vpa(b) + vpa(4)*vpa(c) + ...
vpa(3)*vpa(a))^vpa(3)/(vpa(3)*vpa(b) - vpa(2)*vpa(c) - ...
vpa(a))^vpa(3))/(vpa(4)*sqrt(-(vpa(15)*vpa(b) - vpa(8)*vpa(c) - ...
vpa(6)*vpa(a))/(vpa(2)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) + ...
(-vpa(8)*vpa(b) + vpa(4)*vpa(c) + ...
vpa(3)*vpa(a))^vpa(2)/(vpa(4)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - ...
vpa(a))^vpa(2)) - (-vpa(15)*vpa(b) + vpa(8)*vpa(c) + ...
vpa(6)*vpa(a))/(vpa(6)*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))) - ...
(vpa(3)*vpa(b) - ...
vpa(4)*vpa(c))^vpa(2)/(vpa(3)*vpa(2)^(vpa(2)/vpa(3))*(vpa(3)*vpa(b) - ...
vpa(2)*vpa(c) - vpa(a))*(-vpa(54)*vpa(b)^vpa(3) + ...
vpa(288)*vpa(b)*vpa(c)^vpa(2) - vpa(128)*vpa(c)^vpa(3) - ...
vpa(216)*vpa(b)*vpa(c)*vpa(a) + vpa(108)*vpa(c)*vpa(a)^vpa(2) + ...
sqrt(vpa(23328)*vpa(b)^vpa(5)*vpa(c) - ...
vpa(108864)*vpa(b)^vpa(4)*vpa(c)^vpa(2) + ...
vpa(152064)*vpa(b)^vpa(3)*vpa(c)^vpa(3) - ...
vpa(55296)*vpa(b)^vpa(2)*vpa(c)^vpa(4) + ...
vpa(23328)*vpa(b)^vpa(4)*vpa(c)*vpa(a) - ...
vpa(124416)*vpa(b)^vpa(2)*vpa(c)^vpa(3)*vpa(a) + ...
vpa(55296)*vpa(b)*vpa(c)^vpa(4)*vpa(a) - ...
vpa(11664)*vpa(b)^vpa(3)*vpa(c)*vpa(a)^vpa(2) + ...
vpa(46656)*vpa(b)^vpa(2)*vpa(c)^vpa(2)*vpa(a)^vpa(2) + ...
vpa(62208)*vpa(b)*vpa(c)^vpa(3)*vpa(a)^vpa(2) - ...
vpa(27648)*vpa(c)^vpa(4)*vpa(a)^vpa(2) - ...
vpa(46656)*vpa(b)*vpa(c)^vpa(2)*vpa(a)^vpa(3) + ...
vpa(11664)*vpa(c)^vpa(2)*vpa(a)^vpa(4)))^(vpa(1)/vpa(3))) - ...
(-vpa(54)*vpa(b)^vpa(3) + vpa(288)*vpa(b)*vpa(c)^vpa(2) - ...
vpa(128)*vpa(c)^vpa(3) - vpa(216)*vpa(b)*vpa(c)*vpa(a) + ...
vpa(108)*vpa(c)*vpa(a)^vpa(2) + sqrt(vpa(23328)*vpa(b)^vpa(5)*vpa(c) ...
- vpa(108864)*vpa(b)^vpa(4)*vpa(c)^vpa(2) + ...
vpa(152064)*vpa(b)^vpa(3)*vpa(c)^vpa(3) - ...
vpa(55296)*vpa(b)^vpa(2)*vpa(c)^vpa(4) + ...
vpa(23328)*vpa(b)^vpa(4)*vpa(c)*vpa(a) - ...
vpa(124416)*vpa(b)^vpa(2)*vpa(c)^vpa(3)*vpa(a) + ...
vpa(55296)*vpa(b)*vpa(c)^vpa(4)*vpa(a) - ...
vpa(11664)*vpa(b)^vpa(3)*vpa(c)*vpa(a)^vpa(2) + ...
vpa(46656)*vpa(b)^vpa(2)*vpa(c)^vpa(2)*vpa(a)^vpa(2) + ...
vpa(62208)*vpa(b)*vpa(c)^vpa(3)*vpa(a)^vpa(2) - ...
vpa(27648)*vpa(c)^vpa(4)*vpa(a)^vpa(2) - ...
vpa(46656)*vpa(b)*vpa(c)^vpa(2)*vpa(a)^vpa(3) + ...
vpa(11664)*vpa(c)^vpa(2)*vpa(a)^vpa(4)))^(vpa(1)/vpa(3))/(vpa(6)*vpa(...
2)^(vpa(1)/vpa(3))*(vpa(3)*vpa(b) - vpa(2)*vpa(c) - vpa(a))))))/vpa(2));
end

