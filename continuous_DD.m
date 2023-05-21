function dfdb=continuous_DD(b,n)

    format long
    
    [F,u] = primal(b,n,0);
    
    %% Direct Differentiation
    L=2;
    dx=(L/(n-1));
    
    x=linspace(0,L,n)';
    
    
    f1=zeros(n,1);
    f2=zeros(n,1);
    f3=zeros(n,1);
    f4=zeros(n,1);
    
    J=zeros(n,n);
    
    k=1/dx^2;
    for i = 3:n
        a=k-2*x(i)+3*b(4)*u(i)^2;
        J(i,i)=a;
        J(i,i-1)=-2*k;
        J(i,i-2)=k;
    end
    
    
    %% for b1
    J(1,1) = 1;
    J(2,2) = 1;
    
    f1(1) = 1;
    f1(2) = f1(1);
    
    J = sparse(J);
    
    dudb1 = J\f1;
    
    j=1;
        while x(j)<=0.5
            x_obj(j)=x(j);
            u_obj(j)=u(j);
            dudb1_obj(j) = dudb1(j);
            j=j+1;
        end
    
    y_obj=(u_obj-x_obj.^2).*dudb1_obj;
    
    dFdb1 = trapz(x_obj,y_obj);
    
    %% for b2
    
    J(1,1) = 1;
    J(2,2) = 1;
    
    f2(1) = 0;
    f2(2) = f2(1) + dx;
    
    J = sparse(J);
    
    dudb2 = J\f2;
    
    j=1;
        while x(j)<=0.5
            x_obj(j)=x(j);
            u_obj(j)=u(j);
            dudb2_obj(j) = dudb2(j);
            j=j+1;
        end
    
    y_obj=(u_obj-x_obj.^2).*dudb2_obj;
    
    dFdb2 = trapz(x_obj,y_obj);
     
    %% for b3
    
    J(1,1) = 1;
    J(2,2) = 1;
    
    f3(1) = 0;
    f3(2) = 0;
    
    J = sparse(J);
    
    for i = 3:n
        f3(i) = -(x(i)^2 - x(i));
    end
    
    dudb3 = J\f3;
    
    j=1;
        while x(j)<=0.5
            x_obj(j)=x(j);
            u_obj(j)=u(j);
            dudb3_obj(j) = dudb3(j);
            j=j+1;
        end
    
    y_obj=(u_obj-x_obj.^2).*dudb3_obj;
    
    dFdb3 = trapz(x_obj,y_obj); 
    
    %% for b3
    
    J(1,1) = 1;
    J(2,2) = 1;
    
    f4(1) = 0;
    f4(2) = 0;
    
    J = sparse(J);
    
    for i = 3:n
        f4(i) = -u(i)^3;
    end
    
    dudb4 = J\f4;
    
    j=1;
        while x(j)<=0.5
            x_obj(j)=x(j);
            u_obj(j)=u(j);
            dudb4_obj(j) = dudb4(j);
            j=j+1;
        end
    
    y_obj=(u_obj-x_obj.^2).*dudb4_obj;
    
    dFdb4 = trapz(x_obj,y_obj); 
    
    dfdb=[dFdb1 dFdb2 dFdb3 dFdb4];

end
  
 
 
 
 
 

