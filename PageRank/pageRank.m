% Page Rank
% Corbin Foucart

clear all; close all; clc;

% load data
load('movies.mat');

movies = links(:,1);
actors = links(:,2);

% Constants
alpha = 0.85;

% create matrix B
C = sparse(movies, actors, 1);
B = C*C';

% -------  create matrix P ------- %
% remove diagonal entries, which are self-linking
P = logical(B);
P(logical(eye(size(P)))) = 0;

% normalize by number of entries in each column
sums = 1./sum(P);
sums = diag(sums);

P = P*sums;
P_check = sum(P);

% Plotting Sparsity for P
spy(P);
axis square;

% ---- Solve matrix equation directly ----- %
n = length(P(1,:));
v = ones(n, 1)./n;

A = sparse(eye(n)-alpha*P);
x_direct = A\v;

% normalize direct solution
x_direct = x_direct./norm(x_direct,1);
check_direct = sum(x_direct,1);

% Find the top N values
top_N = 5;
[sorted_dx, sorted_ind] = sort(x_direct, 'descend');
topN_x = sorted_dx(1:top_N)
topN_xind = sorted_ind(1:top_N)
topN_names = movieName(topN_xind)

% ---------------- Jacobi Iteration --------------- %
% Here we have already constructed our pageRank
% matrix, and we simply employ a jacobi iterative
% scheme to come up with our x vector.

% precision
eps = 1e-4;

% initial guess
x_check = v;
error(1) = norm(x_check,1); % error guess to start loop
n_iter = 0;
i = 1;

% iteration
% record error and iteration number
while (error > eps)      
  x_check_next = alpha*P*x_check + v;
  error(i) = abs(sum(x_check_next - x_check))./norm(x_check);
  x_check = x_check_next;
  n_iter = n_iter + 1;
  i = i + 1;
end

last_error = error(end)
x_check_next;
n_iter

% iteration finished
x_jacobi = x_check_next./norm(x_check_next,1);

figure()
plot([1:1:length(error)], error)
title('Error vs. Iteration Number n')
xlabel(' Iteration Number n')
ylabel('Error ||x^{k+1}-x^{k}|| / ||x^{k}||')

% -------- Movie rating etc. plots ---------- %
figure()
scatter(movieRating, x_direct)
xlabel('Movie Rating /10')
ylabel('Page Rank Value')
title('Movie Ratings vs. PageRank Values')

figure()
scatter(movieVotes, x_direct)
xlabel('Movie Votes')
ylabel('Page Rank Value')
title('Movie Votes vs. PageRank Values')










