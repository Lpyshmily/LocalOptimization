% ������С��ʱ��ı仯
clear; clc;
data = dlmread('../source/amplitude_1500.txt');
plot(data(:,1), data(:,2),'.-');