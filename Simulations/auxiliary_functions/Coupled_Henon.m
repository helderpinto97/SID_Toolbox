function [x, y] = Coupled_Henon(N, e, rho)

y = zeros(N + 10000,1);
x = zeros(N + 10000,1);

% Random initial condition
y(1:2) = rand(2,1) * 0.1;
x(1:2) = rand(2,1) * 0.1;

flag = 1;
while flag
  for i=3: N + 10000
    y(i) = 1.4 - y(i-1)^2 + rho(1) * y(i-2);
    x(i) = 1.4 - (e * y(i-1) + (1-e) * x(i-1)) * x(i-1) + rho(2) * x(i-2);
  end
  if sum(isinf(x) +  isinf(y)) > 0    % check if solutions are bounded
    disp('inf')
    y(1:2) = rand(2,1);
    x(1:2) = rand(2,1);
  else
    flag = 0;
  end
end
% Remove the 10000 first data points
y = y(10000:end-1);
x = x(10000:end-1);
