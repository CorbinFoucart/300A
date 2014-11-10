% Page Rank
% Corbin Foucart

clear all; close all; clc;

load('movies.mat');

movies = links(:,1);
actors = links(:,2);

C = sparse(movies, actors, 1);
B = C*C';
spy(C)
axis square

    