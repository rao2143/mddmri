dat = [0.	0.	0.	0.016666666666666666
3.141592653589793	0.	0.	0.016666666666666666
3.141592653589793	3.141592653589793	0.	0.016666666666666666
0.	3.141592653589793	0.	0.016666666666666666
1.5707963267948966	1.5707963267948966	3.141592653589793	0.016666666666666666
0.	1.5707963267948966	4.71238898038469	0.016666666666666666
3.141592653589793	1.5707963267948966	1.5707963267948966	0.016666666666666666
4.71238898038469	1.5707963267948966	0.	0.016666666666666666
3.141592653589793	1.5707963267948966	4.71238898038469	0.016666666666666666
1.5707963267948966	1.5707963267948966	0.	0.016666666666666666
4.71238898038469	1.5707963267948966	3.141592653589793	0.016666666666666666
0.	1.5707963267948966	1.5707963267948966	0.016666666666666666
1.0172219678978518	0.6283185307179591	4.158814621487645	0.016666666666666666
2.77672882547631	1.0471975511965974	2.77672882547631	0.016666666666666666
2.1243706856919413	1.2566370614359172	4.158814621487645	0.016666666666666666
2.124370685691942	0.6283185307179586	2.124370685691942	0.016666666666666666
1.0172219678978518	1.8849555921538759	4.158814621487645	0.016666666666666666
2.77672882547631	2.094395102393195	3.506456481703276	0.016666666666666666
0.3648638281134833	1.0471975511965974	3.5064564817032764	0.016666666666666666
2.1243706856919413	1.8849555921538759	2.1243706856919413	0.016666666666666666
2.1243706856919413	2.5132741228718345	4.158814621487645	0.016666666666666666
1.0172219678978518	1.2566370614359172	2.124370685691942	0.016666666666666666
0.3648638281134833	2.094395102393195	2.77672882547631	0.016666666666666666
1.0172219678978518	2.5132741228718345	2.124370685691942	0.016666666666666666
4.158814621487645	0.6283185307179582	1.0172219678978518	0.016666666666666666
3.506456481703276	1.0471975511965974	3.506456481703276	0.016666666666666666
1.017221967897851	1.2566370614359172	5.265963339281734	0.016666666666666666
5.265963339281735	0.6283185307179586	5.265963339281735	0.016666666666666666
2.1243706856919413	1.8849555921538759	5.265963339281734	0.016666666666666666
3.506456481703276	2.094395102393195	2.77672882547631	0.016666666666666666
5.918321479066103	1.0471975511965974	2.77672882547631	0.016666666666666666
1.0172219678978518	1.8849555921538759	1.0172219678978518	0.016666666666666666
5.265963339281734	2.5132741228718345	1.0172219678978518	0.016666666666666666
2.1243706856919413	1.2566370614359172	1.0172219678978518	0.016666666666666666
5.918321479066103	2.094395102393195	3.506456481703276	0.016666666666666666
4.158814621487644	2.5132741228718345	5.265963339281734	0.016666666666666666
2.1243706856919413	0.6283185307179586	5.265963339281734	0.016666666666666666
5.918321479066103	1.0471975511965974	5.918321479066103	0.016666666666666666
4.158814621487644	1.2566370614359172	2.124370685691942	0.016666666666666666
1.017221967897851	0.6283185307179586	1.017221967897851	0.016666666666666666
5.265963339281734	1.8849555921538759	2.1243706856919413	0.016666666666666666
5.918321479066103	2.094395102393195	0.3648638281134833	0.016666666666666666
3.506456481703276	1.0471975511965974	0.3648638281134833	0.016666666666666666
4.158814621487644	1.8849555921538759	4.158814621487644	0.016666666666666666
1.0172219678978518	2.5132741228718345	5.265963339281734	0.016666666666666666
5.265963339281735	1.2566370614359172	4.158814621487645	0.016666666666666666
3.506456481703276	2.094395102393195	5.918321479066103	0.016666666666666666
2.1243706856919413	2.5132741228718345	1.0172219678978518	0.016666666666666666
5.265963339281734	0.6283185307179586	2.1243706856919413	0.016666666666666666
0.3648638281134833	1.0471975511965974	0.3648638281134833	0.016666666666666666
5.265963339281735	1.2566370614359172	1.017221967897851	0.016666666666666666
4.158814621487644	0.6283185307179582	4.158814621487644	0.016666666666666666
4.158814621487645	1.8849555921538759	1.0172219678978518	0.016666666666666666
0.3648638281134833	2.094395102393195	5.918321479066103	0.016666666666666666
2.77672882547631	1.0471975511965974	5.918321479066103	0.016666666666666666
5.265963339281734	1.8849555921538759	5.265963339281734	0.016666666666666666
4.158814621487645	2.5132741228718345	2.124370685691942	0.016666666666666666
4.158814621487645	1.2566370614359172	5.265963339281734	0.016666666666666666
2.77672882547631	2.094395102393195	0.36486382811348284	0.016666666666666666
5.265963339281734	2.5132741228718345	4.158814621487645	0.016666666666666666];

alpha = dat(:,1);
beta = dat(:,2);
gamma = dat(:,3);
weight = dat(:,4);

theta = beta;
phi = alpha;

x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);

figure(1), clf
subplot(2,2,1)
plot(alpha,beta,'o')
subplot(2,2,2)
plot3(alpha,beta,gamma,'o')
subplot(2,2,3)
plot3(x,y,z,'o')
axis('square')
xlabel('x'), ylabel('y'), zlabel('z')

save POLYTOPE60 alpha beta gamma weight

