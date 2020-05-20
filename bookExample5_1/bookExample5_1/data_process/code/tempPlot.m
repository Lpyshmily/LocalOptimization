% 读取积分过程中的变量，将春分点轨道根数转化为直角坐标，并画图
clear; clc;
data = dlmread('../source/temp_1100.txt');
edata = data(:, 2:7);
sizeData = size(edata);
rv = zeros(sizeData);
mu = 1;

% 转换过程
for i=1:sizeData(1)
    p = edata(i, 1);
    f = edata(i, 2);
    g = edata(i, 3);
    h = edata(i, 4);
    k = edata(i, 5);
    L = edata(i, 6);
    h_ = sqrt(p/mu);
    n = h_ / (1+f*cos(L)+g*sin(L));
    s2 = 1+h*h+k*k;
    hh_kk = h*h - k*k;
    hk2 = 2*h*k;
    r = n*h_*mu;
    rv(i, 1) = r/s2*(cos(L)+hh_kk*cos(L)+hk2*sin(L));
    rv(i, 2) = r/s2*(sin(L)-hh_kk*sin(L)+hk2*cos(L));
	rv(i, 3) = 2.0*r/s2*(h*sin(L)-k*cos(L));
	rv(i, 4) = -1.0/h_/s2*(sin(L)+hh_kk*sin(L)-hk2*cos(L)+g-hk2*f+hh_kk*g);
	rv(i, 5) = -1.0/h_/s2*(-cos(L)+hh_kk*cos(L)+hk2*sin(L)-f+hk2*g+hh_kk*f);
	rv(i, 6) = 2.0/h_/s2*(h*cos(L)+k*sin(L)+f*h+g*k);
end

figure(1)
plot3(rv(:,1), rv(:,2), rv(:,3));
figure(2)
plot(rv(:,1), rv(:,2))