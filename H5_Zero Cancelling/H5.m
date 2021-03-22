% David Pipes
% H6 Code to do matrix math

clc, clear, close all
%% Set up matricies
syms s K1 K2 K3 K4 K5 %regulator
syms s L1 L2 L3 L4 %estimator
A = [0,0,1,0;0,0,0,1;-1,1,-1,1;1,-1,1,-1];
B = [0;0;1;0];
C = [0,1,0,0];
Ak = [0,0,1,0,0;0,0,0,1,0;-1,1,-1,1,0;1,-1,1,-1,0;0,1,0,0,0];
Bk = [0;0;1;0;0];
sI = [s,0,0,0;0,s,0,0;0,0,s,0;0,0,0,s];
sIk = [s,0,0,0,0;0,s,0,0,0;0,0,s,0,0;0,0,0,s,0;0,0,0,0,s];
K = [K1,K2,K3,K4,K5];
L = [L1;L2;L3;L4];


Ao = Ak-Bk*K; 
Aoo = A-L*C;

dtrK = det(sIk - Ao); 
dtrL = det(sI - Aoo); 

h = 0.001; 
t = 0:h:10;

U = NaN(4,length(t));
Y = NaN(4,length(t));

K1 = [-3094, 1802, 28, 3428];
K2 = 7776;

KK1 = [22,842,23,193];
KK2 = 1296;

row1 = [A-B*K1, -B*K2, -B*K1];
row2 = [C 0 0 0 0 0];
ffz = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
row3 = [ffz A-L*C];
Aoverall = [row1; row2; row3];
Boverall = [B;0;0;0;0;0];
Coverall = row2;

%%
y = [0;
    0;
    0;
    0];
s = 0; 
r = 1;

for n = 0:length(t)-1
    
    u = -K1*y-K2*s; % needed to implement x(hat) for y here
    s = s+h*(y-r);
    
    U(:,n+1) = u;
    Y(:,n+1) = y;
    
    y = y+h*(y+u);
end

%%
UU = NaN(4,length(t));
YY = NaN(4,length(t));

y = [0;
    0;
    0;
    0];
s = 0; 
r = 1;

for n = 0:length(t)-1
    
    u = -KK1*y-KK2*s; % needed to implement x(hat) for y here
    s = s+h*(y-r);
    
    UU(:,n+1) = u;
    YY(:,n+1) = y;
    
    y = y+h*(y+u);
end




%%
subplot(411), plot(t,Y)
title("Y versus t, no zero canceling")
subplot(412), plot(t,U)
title("U versus t, no zero canceling")
subplot(413), plot(t,YY)
title("Y versus t, zero canceling")
subplot(414), plot(t,UU)
title("U versus t, zero canceling")
