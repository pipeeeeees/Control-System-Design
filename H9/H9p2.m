% David Pipes
% H9

clc, clear, close all
%% Set-up Matricies

syms z K1 K2 K3 x L1 L2

A = [0, 1;
     0, 0];
B = [0;
     1];
C = [1, 0];

lambda_e = 12;
lambda_gamma = 3;
T = 0.096;

zI = [z, 0;
      0, z];
zI3x3 = [z, 0, 0;
         0, z, 0
         0, 0, z];
  
K = [K1, K2, K3];
Kone = [K1, K2];
Ktwo = K3;
curly_K = [Kone, Ktwo];

L = [L1;
     L2];

%% Problem 1b code - indirect digital design

Ad = expm(A*T);
Bd = [0.5*(T^2);
      T]; 

polynomial_L = det(zI - (Ad - L*C));
gamma_e = expm(-lambda_e*T);
L1 = 1.368;
L2 = 4.874;
L = [L1;
     L2];

curly_Ad = [Ad, zeros(2,1);
           T*C, 1];
curly_Bd = [Bd;
            0];

polynomial_K = det(zI3x3 - (curly_Ad - curly_Bd*curly_K));
gamma_gamma = expm(-lambda_gamma*T);
K1 = 19.585;  
K2 = 6.88;   
K3 = 19.56; 
Kone = [K1, K2];
Ktwo = K3;
curly_K = [Kone, Ktwo];
%polynomial_K = vpa(det(zI3x3 - (curly_Ad - curly_Bd*curly_K))) % to check
%if the K values are correct


% Answers to B
F = [Ad - Bd*Kone - L*C, -Bd*Ktwo;
     zeros(1,2) , 1];

G = [L, zeros(2,1);
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
sgtitle("Homework 9 - Problem 2C")