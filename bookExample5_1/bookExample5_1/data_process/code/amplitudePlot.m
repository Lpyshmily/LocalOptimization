% 推力大小随时间的变化
clear; clc;
data = dlmread('../source/amplitude_1500.txt');
plot(data(:,1), data(:,2),'.-');