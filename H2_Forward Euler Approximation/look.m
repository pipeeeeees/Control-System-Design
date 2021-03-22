% David Pipes
% ECE 4550 H2 Problem 2 Summer 2020

clear

%% time and constants
pi = 3.14159265;

h = 0.01;           % time increment
t = 0:h:10;         % total time index
sintime = 0:h:1.99; % time spent as a sine wave
zeros = 0:h:8;      % time spent at zero

%% creating u(t)

aa = ones(1,length(sintime));
bb = 0*ones(1,length(zeros));

% filling in 'aa' with the sine values
for k = 1:length(sintime)
    aa(1,k) = sin(pi*sintime(k));
end

% appending the 2 sec sine with the 8 sec zeros
u = [aa,bb];

%% setting up x
x = NaN*ones(4,length(t));

% first column of x set to initial value x
x(:,1) = [1; 0; 1; 0];

A = [0,1,0,0;-1,-1,1,1;0,0,0,1;1,1,-1,-1]; B = [0;1;0;0]; C = [0,0,1,0];

for n = 1:length(t)-1
    x(:,n+1) = x(:,n)+h*(A*x(:,n)+B*u(:,n));
end
y = C*x;

%% plots

subplot(2,1,1);
plot(t,y,'linewidth',1)
xlabel('t'), ylabel('y'), title('Forward Euler Approximation')

subplot(2,1,2);
plot(t,u,'linewidth',1)
xlabel('t'), ylabel('u(t)'), title('u(t) over time (t)')

