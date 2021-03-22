% David Pipes
% H9

clc, clear, close all
%% Set-up Matricies

syms s K1 K2 K3
syms s L1 L2

A = [0, 1;
     0, 0];
B = [0;
     1];
C = [1, 0];

lambda_e = 12;
lambda_gamma = 3;
T = 0.096;

sI = [s, 0;
      0, s];
sI3x3 = [s, 0, 0;
         0, s, 0
         0, 0, s];
  
K = [K1, K2, K3];
Kone = [K1, K2];
Ktwo = K3;
curly_K = [Kone, Ktwo];

L = [L1;
     L2];

%% Problem 1b code - indirect digital design

polynomial_L = det(sI - (A - L*C));
L_value_polynomial = (s + lambda_e)^2;
% from here we know L2 = 144, L1 = -120.
L1 = 24;
L2 = 144;
L = [L1;
     L2];

curly_A = [0, 1, 0;
           0, 0, 0;
           1, 0, 0];
curly_B = [0;
           1;
           0];

polynomial_K = det(sI3x3 - (curly_A - curly_B*curly_K));
K_value_polynomial = (s + lambda_gamma)^3;
% from here we know K3 = 27, K2 = 9, K1 = 0.
K3 = 27;
K2 = 9;
K1 = 27;

Kone = [K1, K2];
Ktwo = K3;

% Answers to B
F = [1+T*A-T*B*Kone-T*L*C, -T*B*Ktwo;
     zeros(1,2) , 1];

G = [T*L, zeros(2,1);
     T,-T];
H = [-Kone, -Ktwo];

%% Problem 1c

h = T/100;
t = 0:h:4;
x = [0;
     0];
v = [0;
     0;
     0];
U = NaN(1,length(t));
Y = NaN(1,length(t));


for n = 0:length(t)-1
    
    % discrete-time update 
    if mod(n*h,T) == 0
        
        % sensor measurements
        y = C*x;
        r = 1;
        
        % controller output
        u = H*v;
        
        % controller dynamics
        v = F*v+G*[y;r];
    end
    
    % continuous-time update
    
    % data logger
    Y(n+1) = C*x;
    U(n+1) = u;
    
    % plant dynamics
    x = x+h*(A*x+B*u);
end

subplot(211), plot(t,Y), xlabel('t'), ylabel('y')
subplot(212), stairs(t,U), xlabel('t'), ylabel('u')







