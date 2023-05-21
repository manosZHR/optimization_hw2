function dfdb=discrete_adj(b,n)

    L1=0.5;
    L2=2;
    dx=(L2/(n-1));
    k=1/dx^2;
    x=linspace(0,L2,n)';
    
    [F,u] = primal(b,n,0);
    
    j=1;
    while x(j)<=L1
        x_obj(j)=x(j);
        u_obj(j)=u(j);
        j=j+1;
    end
    j=j-1;
    x_obj=x_obj';
    u_obj=u_obj';
    
    %% Discrete Adjoint
    
    f1=zeros(j,1);
    J1=zeros(j,j);
    psi_T=zeros(j,1);
    
    k=1/dx^2;
    for i = 3:j
        a=k-2*x(i)+3*b(4)*u(i)^2;
        J1(i,i)=a;
        J1(i,i-1)=-2*k;
        J1(i,i-2)=k;
    end
    
    
    J1(1,1) = -1;
    J1(2,2) = -1;
    
    f1(1)=dx/2*(u_obj(1)-x_obj(1)^2);
    f1(j)=dx/2*(u_obj(j)-x_obj(j)^2);
    
    for i=2:j-1
    
        f1(i)=dx*(u_obj(i)-x_obj(i)^2);
    
    end
    
    
    psi_T=J1'\(-f1);
    
    for i=j+1:n
    
        psi_T(i)=0;
    
    end
    
    psi=psi_T';
    
    f2=zeros(n,length(b));
    f2(1,1)=1;
    f2(1,2)=0;
    f2(1,3)=0;
    f2(1,4)=0;
    
    f2(2,1)=1;
    f2(2,2)=dx;
    f2(2,3)=0;
    f2(2,4)=0;
    
    for i=3:n
    
        f2(i,1)=0;
        f2(i,2)=0;
        f2(i,3)=(x(i)^2-x(i));
        f2(i,4)=u(i)^3;
    
    end
    
    
    for m=1:4
    
        dfdb(:,m)=psi*f2(:,m);
    
    end
    
end







    
