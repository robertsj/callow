n = 5; A = 0.5*speye(n) - 0.2*diag(ones(n-1, 1), -1) - 0.2*diag(ones(n-1, 1), 1); 
L = tril(A)-diag(diag(A));
U = triu(A)-diag(diag(A));
D = diag(diag(A));
b = ones(n, 1);
m = 1;
A

%% richardson
%I = speye(n);
%x = zeros(n, 1);
%for i = 1:m
%  x1 = (I - A)*x + b;
%  r1(i) = norm(x1-x);
%  x = x1;
%end

% jacobi
x = zeros(n, 1);
xa = -L*x;
xb = xa - U*x;
xc = xb + b;
x  = D \ xc; [xa xb xc x]
xa = -L*x;
xb = xa - U*x;
xc = xb + b;
x  = D \ xc;  [xa xb xc x]
xa = -L*x;
xb = xa - U*x;
xc = xb + b;
x  = D \ xc;  [xa xb xc x]
%for i = 1:m
%  x1 = D \ ( -(L+U)*x + b);
%  r2(i) = norm(x1-x);
%  x = x1;
%end

%% gauss seidel
%x = zeros(n, 1);
%for i = 1:m
%  x1 = (D + L) \ ( -U*x + b);
%  r3(i) = norm(x1-x);
%  x = x1;
%end
%r2(1:3)
%semilogy(1:m, r1, 'k', 'LineWidth', 2, 1:m, r2, 'b--', 'LineWidth', 2, 1:m, r3, 'r-.', 'LineWidth', 2), grid on
%x
