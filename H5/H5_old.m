% David Pipes
% H5 Code to do matrix math

clc, clear, close all
%% Set up matricies
syms s K1 K2 K3 K4
syms s L1 L2 L3 L4
A = [0,0,1,0;0,0,0,1;-2,2,0,0;2,-2,0,0];
B = [0;0;1;0];
C = [0,1,0,0];
sI = [s,0,0,0;0,s,0,0;0,0,s,0;0,0,0,s];
K = [K1,K2,K3,K4];
L = [L1;L2;L3;L4];

%% Check for controllability and observability
Ao = A-B*K;
Aoo = A-L*C;

dtrC = det(sI - Ao);    % makes characteristic polynomial for 
dtrO = det(sI - Aoo);   % if not equal to 0, it is observable

c1 = A*B;
c2 = (A^2)*B;
c3 = (A^3)*B;
bigC = [B,c1,c2,c3];        % Sets up the controllable matrix
dtrControllable = det(bigC);% if not zero, system is controllable

o1 = C*A;
o2 = C*(A^2);
o3 = C*(A^3);
bigO = [C;o1;o2;o3];        % Sets up the observable matrix
dtrObservable = det(bigO);  % if not zero, system is observable

