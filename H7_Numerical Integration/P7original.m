% David Pipes
% P7

clc, clear, close all
%% Set up matricies
syms s K1 K2 K3 K4 K5
syms s L1 L2 L3 L4 
A = [ 0, 0, 1, 0;
      0, 0, 0, 1;
     -1, 1,-1, 1;
      1,-1, 1,-1];
B = [ 0;
      0;
      1;
      0];
C = [ 0, 1, 0, 0];
E = [ 0;
      0;
      0;
     -1];
Kone = [22,842,23,193];
Ktwo = 1296;
sI = [s,0,0,0;
      0,s,0,0;
      0,0,s,0;
      0,0,0,s];
K = [Kone,Ktwo];
L = [-13870;
     46;
     19966;
     770];

%% Question 1

planttransfer = C*((sI-A)^-1)*B;

M = [C*(A^0);
     C*(A^1);
     C*(A^2)];
N = [1, 1, 0, 1];   % Selected N so NB is zero, det(S) =/= 0
S = [M;N];

testingNB = N*B;        % should be 0 to work
testingdetS = det(S);   % should be non-zero

SInverse = S^-1;

Q = SInverse([1 2 3 4],[1 2 3]);
R = SInverse([1 2 3 4],4);

Ainv = N*A*R;                   % for question 1
Binv = [N*A*Q, 0];              % for question 1

gamma = 1/(C*(A^2)*B);

Cinv = [-gamma*C*(A^3)*R; R];   % for question 1
Dinv = [-gamma*C*(A^3)*Q, gamma;% for question 1
        Q          , [0;0;0;0]];

%% Question 2

% time grid
h = 0.001;
t = 0:h:5;

% memory allocation
uinv = NaN*zeros(4,length(t));
xinv = NaN*zeros(1,length(t)); 
yref = NaN*zeros(1,length(t));
ydotref = NaN*zeros(1,length(t));
ydotdotref = NaN*zeros(1,length(t));
ydotdotdotref = NaN*zeros(1,length(t));
uref = NaN*zeros(1,length(t));
xref = NaN*zeros(4,length(t));

% initial condition
xinv(1) = 0;

for i = 1:length(t)
    yref(i) = (1 - exp(-0.8*t(i)))*sin(2*t(i));
    ydotref(i) = 0.8*exp(-0.8*t(i))*sin(2*t(i)) + (2 - 2*exp(-0.8*t(i)))*cos(2*t(i));
    ydotdotref(i) = exp(-0.8*t(i))*((3.36 - 4*exp(0.8*t(i)))*sin(2*t(i)) + 3.2*cos(2*t(i)));
    ydotdotdotref(i) = exp(-0.8*t(i))*((4.16-8*exp(0.8*t(i)))*cos(2*t(i)) - 9.088*sin(2*t(i)));
    uinv(1,i) = yref(i); 
    uinv(2,i) = ydotref(i);
    uinv(3,i) = ydotdotref(i);
    uinv(4,i) = ydotdotdotref(i);
end

% numerical integration
for i = 1:length(t)-1
    xinv(i+1) = xinv(i) + h*(Ainv*xinv(i)+Binv*uinv(:,i));
end

% output evaluation
yinv = Cinv*xinv + Dinv*uinv;

for i = 1:length(t)
    uref(i) = yinv(1,i);    % uref from yinv
    xref(1,i) = yinv(2,i);
    xref(2,i) = yinv(3,i);
    xref(3,i) = yinv(4,i);
    xref(4,i) = yinv(5,i);
end

%% Question 3

x = NaN*zeros(4,length(t)); 
u = NaN*zeros(1,length(t)); 
y = NaN*zeros(1,length(t)); 
xhat = NaN*zeros(4,length(t)); 
sigma = NaN*zeros(1,length(t)); 

w = 5;          % constant disturbance signal
x(:,1) = [0;    % initialize x
          0;
          0;
          0];
y(1) = C*x(:,1);% initialize y
xhat(:,1) = 0;  % initialize xhat
sigma(1) = (uref(1) + Kone*xref(:,1))/Ktwo;

% numerical integration
for i = 1:length(t)-1
    sigma(i+1) = sigma(i) + h*(y(i)-yref(i));                              % creates new sigma
    u(i) = uref(i) - Kone*(xhat(:,i)-xref(:,i)) - Ktwo*sigma(i);           % creates a new value for u
    xhat(:,i+1) = xhat(:,i) + h*( A*xhat(:,i) + B*u(i) - L* (C*xhat(:,i) - y(i))); % creates new xhat
  
    x(:,i+1) = x(:,i) + h*(A*x(:,i)+B*u(i)+E*w);                           % creates new x 
    y(i+1) = C*x(:,i+1);                                                   % creates new y 
end


%% Plotting

figure(1)
subplot(3,1,1)
plot(t,yref)
title('y_r_e_f(t)')
xlabel('Time (t)')
ylabel('y_r_e_f(t)')

subplot(3,1,2)
plot(t,uref)
title('u_r_e_f(t)')
xlabel('Time (t)')
ylabel('u_r_e_f(t)')

subplot(3,1,3)
plot(t,xref)
title('x_r_e_f(t)')
xlabel('Time (t)')
ylabel('x_r_e_f(t)')

figure(2)
subplot(3,1,1)
plot(t,yref-y)
title('y_r_e_f(t) - y')
xlabel('Time (t)')
ylabel('y_r_e_f(t) - y')

subplot(3,1,2)
plot(t,uref-u)
title('u_r_e_f(t) - u')
xlabel('Time (t)')
ylabel('u_r_e_f(t) - u')

subplot(3,1,3)
plot(t,xref-x)
title('x_r_e_f(t)- x')
xlabel('Time (t)')
ylabel('x_r_e_f(t) - x')