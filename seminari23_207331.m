% TDD_2020 Seminar 2 & 3 - 207331
%% Part 1
A = 3;
nMc = 100000;
x = A*(rand(1,nMc)<= 0.5);
z = exprnd(1,1,nMc); 
y = x+z;

eq_distance = zeros(1,nMc);
count = 0;
%Bit probability error 
for i = 1 : nMc
    eq_distance = sqrt((x(:)-y(i)).^2);
    if (min(eq_distance) == sqrt((x(i)-y(i))^2))
        count = count +1;
    end
    
end
Pe = 1 - (count/nMc);

%% Part 2 
A = 3;
nMc = 100000;
x = A*(rand(1,nMc)<= 0.5);
z = exprnd(1,1,nMc); 
y = x+z;
yest = zeros(1,nMc);
count = 0;
threshold = linspace(0,A);
Pe = zeros(1,length(threshold));
%Bit probability error 
for j=1:length(threshold)
    for i = 1 : nMc
        yest(i) = A*(threshold(j)<=y(i));
    end
    Pe(j) = sum(yest ~= x)/nMc;
end

plot(threshold,Pe)
xlabel('Threshold') 
ylabel('Bit error probability') 
title('Estimation of error probability')