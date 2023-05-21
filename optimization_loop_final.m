clear; close all; clc

tic;

%% Initialize design variables

    b0 = [0.1 0.2 0.3 0.7];

%% Choose discretization for primal

    n=1000;

%% Choose method for dfdb
    % 1.FD, 
    % 2.COMPLEX VAR 
    % 3.CONTINUOUS ADJOINT 
    % 4.CONTINUOUS DIRECT DIFFERENTIATION 
    % 5.DISCRETE ADJOINT
    % 6.DISCRETE DIRECT DIFFERENTIATION 

    method = 4;

%% Choose parameters of optimization loop
    % maximum #iterations (maxiter)
    % tolerance for Fnew-Fold (tol_1)
    % tolerance for F (tol_2)
    % h

    maxiter=5000;
    tol_1=1e-9;
    tol_2=1e-7;
    h=1.9;

%% Optimization loop

err=zeros(maxiter,1);

b_new=b0;
iter=1;
while true
    
    b_old=b_new;
    [F_old,~]=primal(b_old,n,0);

    switch method
        case 1
            dfdb=FD(b_old,n);
            dfdb=dfdb';
        case 2
            dfdb=complex_der(b_old,n);
            dfdb=dfdb';
        case 3
            dfdb=continuous_adj(b_old,n);
        case 4
            dfdb=continuous_DD(b_old,n);
        case 5
            dfdb=discrete_adj(b_old,n);
        case 6
            dfdb=continuous_DD(b_old,n);
    end
    
    b_new = b_old-h*dfdb;

    [F_new,~]=primal(b_new,n,0);
    err(iter) = abs(F_new-F_old);
    
    if iter>=maxiter
        break; end
    if err(iter)<tol_1
        break;end
    if F_new<tol_2
    break;end

    iter=iter+1;

end
        
%% Display results

% dfdb comparisson for different methods 

dfdb_show=zeros(5,4);
for i=1:5
    
    dfdb_show(i,:)=der_f(b0,n,i);
    
end
dfdb_show

[F0,~]=primal(b0,n,0);

fprintf('Objective Function for b0: %4.7f \n',F0');
fprintf('Objective Function for b_new: %4.7f \n',F_new');
fprintf('Total number of iterations: %d \n',iter');
toc;

figure(1);
x=1:iter;
y=err(1:iter);
semilogy(x,y,'.r','LineWidth',2)
xlabel('#iteration'); ylabel('|F_{new} - F_{old}|')
grid on; box on; axis tight