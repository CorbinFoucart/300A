% Thomas Algorithm Test Function
% Corbin Foucart

clear all; close all; clc;


N = [4, 400, 1000];

for i = 1:length(N)
    
% Heat equation Matrix
main_diag = -2*ones(N(i) - 1,1);
off_diag = ones(N(i)-2, 1);
A = full(gallery('tridiag', off_diag, main_diag, off_diag));

% RHS
b = zeros(N(i)-1, 1);
b(end) = -2;

x = thomas_algorithm(A, b);
x_check = A\b;

% output to make sure our solutions are equal
% if true, returns N - 1
truth = sum(x == x_check)
fprintf('Should return N - 1')
fprintf(' \n \n')
fprintf('Norm results: ')
norm_thomas = norm(x)
norm_matlab = norm(x_check)

end

 