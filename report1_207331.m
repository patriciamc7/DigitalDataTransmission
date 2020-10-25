A = 4/6;
n = 10000;
decoded = zeros(1,n);
good_decoded = 0;
min_distance= ones(n,4);
xn = rand(1,n); %vector continuous
zn = zeros(1,n);
out1 = 0; %0+0i
out2 = A; %A+0i
out3 = 1j*A; %0+Ai
out4 =  A+1j*A; %A+Ai

constellation = [out1, out2, out3 , out4];
xn(xn<=0.25) = out1;
xn((xn>0.25) & (xn<=0.5)) = out2;
xn((xn>0.5) & (xn<=0.75)) = out3;
xn((xn>0.75) & (xn<=1)) = out4;
%xn vector discrete with values A or 0

for i=1:n
zn(i) = exp(-i)+ (1j*exp(-i));
end
yn = xn+zn;
%Decision regions
for i=1:n
    for j=1:4
        min_distance(i,j)  = abs(yn(i) - constellation(j))^2;
    end
    [M,I] = min((min_distance(i,:)));
    decoded(i) = constellation(I);
end

for i=1:n
    if (decoded(i) == xn(i))
     good_decoded = good_decoded + 1; 
    end
end
prob_error = 1-(good_decoded/n);