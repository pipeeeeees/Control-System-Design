% David Pipes
% P8
clc, clear, close all
%%
a = 1.005;
m = 360; 
h = 0.1;
c = 6:h:10;
b = ones(1,length(c));


top = (a*c)/(1-a);
bottom = (((a*b)/(a-1)) * a^(m-1)) + (((a*b) + (a*c))/(1-a));

exp = top./bottom;

i = log10(exp)/log10(a);
plot(c, i);
title("Retirement era vs. ratio c/b")
xlabel("c/b")
ylabel("Retirement era (months)")