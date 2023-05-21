clear; close all; clc
tic;

%% Design Variables

b = [0.1; 0.2; 0.3; 0.7];

%% Grid with n = 500 nodes
n=500;

[F_500,u_500] = primal(b,n,1);


%% Grid with n = 1000 nodes
n=1000;

[F_1000,u_1000] = primal(b,n,1);


%% Grid with n = 2000 nodes
n=2000;

[F_2000,u_2000] = primal(b,n,1);


%% Grid with n = 3000 nodes
n=3000;

[F_3000,u_3000] = primal(b,n,1);


%% Grid with n = 4000 nodes
n=4000;

[F_4000,u_4000] = primal(b,n,1);

%% Grid with n = 5000 nodes
n=5000;

[F_5000,u_5000] = primal(b,n,1);

%% Grid with n = 6000 nodes
n=6000;

[F_6000,u_6000] = primal(b,n,1);

%% Grid with n = 7000 nodes
n=7000;

[F_7000,u_7000] = primal(b,n,1);


toc;

