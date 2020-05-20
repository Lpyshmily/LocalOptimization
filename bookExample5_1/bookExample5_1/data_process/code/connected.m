clear; clc;
data1 = dlmread('../source/temp_1500.txt');
amplitude_data = dlmread('../source/amplitude_1500.txt');

timeList = data1(:,1);
time_length = size(data1, 1);
edata = data1(:, 2:7);
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


% 寻找推力大小跳变的位置
jump_index = 1;
current_state = amplitude_data(jump_index,2);
while(1)
    jump_index = jump_index + 1;
    if (amplitude_data(jump_index,2)) ~= current_state
        break;
    end
end
jump_time = amplitude_data(jump_index,1);

start_index = 1;
stop_index = 0;
flag = 0;
figure(1)
for i=1:time_length
    
    if timeList(i)>jump_time
        % stop_index = i-1;
        stop_index = i;
        if flag==0
            plot(rv(start_index:stop_index,1), rv(start_index:stop_index,2), 'b');
            hold on
            flag = 1;
        else
            plot(rv(start_index:stop_index,1), rv(start_index:stop_index,2), 'r');
            hold on
            flag = 0;
        end
        start_index = i;
        % 寻找下一个跳变点
        current_state = amplitude_data(jump_index,2);
        while(1)
            jump_index = jump_index + 1;
            if (jump_index<=size(amplitude_data,1) && amplitude_data(jump_index,2)) ~= current_state
                break;
            end
        end
        jump_time = amplitude_data(jump_index-1,1);
    end 
end
hold off