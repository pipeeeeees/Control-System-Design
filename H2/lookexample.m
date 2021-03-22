h = 0.001; t = 0:h:20; 2

A = [-1,-1;1,0]; B = [1;0]; C = [1,0];
u = -1*ones(1,length(t));
x = NaN*ones(2,length(t));
x(:,1) = [0; 1];

for n = 1:length(t)-1
    x(:,n+1) = x(:,n)+h*(A*x(:,n)+B*u(:,n));
end

y = C*x;

subplot(2,1,1);
plot(t,y,'linewidth',1)
xlabel('t'), ylabel('y'), title('Forward Euler')

subplot(2,1,2);
plot(t,u,'linewidth',1)
xlabel('t'), ylabel('u(t)'), title('u(t) over time (t)')
