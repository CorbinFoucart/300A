% HW 6 Problem 2
% Corbin Foucart

close all; clc;

A = [1 2 3;
     4 8 1;
     6 2 9;
     2 7 6]

answer = inv(A'*A)

[Q, R] = qr(A);

R

cols = length(R(1,:))
rows = length(R(:,1))
n = min(cols, rows)

X = zeros(cols, cols);
X(cols, cols) = 1/R(cols, cols).^2;

X = AtAContruct(R,X,n)

difference = X - answer



