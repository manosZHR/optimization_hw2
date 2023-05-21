function dfdb=continious_adj(b,nnx)

    L1=0.5;
    L2=2;
    dx=(L2/(nnx-1));
    k=1/dx^2;
    x=linspace(0,L2,nnx)';
    
    [F_prim,u_prim] = primal(b,nnx,0);
    
    j=1;
    while x(j)<=L1
        x_1(j)=x(j);
        u_1(j)=u_prim(j);
        j=j+1;
    end
    j=j-1;
    
    x_1=x_1';
    u_1=u_1';
    
    x_2=x(j+1:nnx);
    u_2=u_prim(j+1:nnx);
    
    f1=zeros(j,1);
    A1=zeros(j,j);
    f2=zeros(nnx-j,1);
    A2=zeros(nnx-j,nnx-j);
    
    m=(nnx-2)-j;
    for i=1:m
        a=k+3*b(4)*u_2(i)^2-2*x_2(i);
        A2(i,i)=a;
        A2(i,i+1)=-2*k;
        A2(i,i+2)=k;
        f2(i)=0;
    end
    
    A2(m+2,m+2)=1;
    f2(m+2)=0;
    A2(m+1,m+1)=1;
    f2(m+1)=f2(m+2);
    
    A2=sparse(A2);
    psi_2=A2\f2;
    
    
    for i=1:j-2
        a=k+3*b(4)*u_1(i)^2-2*x_1(i);
        A1(i,i)=a;
        A1(i,i+1)=-2*k;
        A1(i,i+2)=k;
        f1(i)=x_1(i)^2-u_1(i);
    end
    
    
    A1(j,j)=1;
    f1(j)=( x_1(j)^2-u_1(j)+2*k*psi_2(1)-k*psi_2(2) ) / ( k+3*b(4)*u_1(j)^2-2*x_1(j) );
    
    A1(j-1,j-1)=1;
    f1(j-1)= ( x_1(j-1)^2-u_1(j-1)+2*k*f1(j)-k*psi_2(1) ) / ( k+3*b(4)*u_1(j-1)^2-2*x_1(j-1) );
    
    A1=sparse(A1);
    psi_1=A1\f1;
    
    psi=cat(1,psi_1,psi_2);
    
    dfdb1 = (psi(2)-psi(1))/dx;
    dfdb2 = -psi(1);
    y=(x.^2-x).*psi;
    dfdb3 = trapz(x,y);
    y=psi.*(u_prim.^3);
    dfdb4 = trapz(x,y);

    dfdb=[dfdb1 dfdb2 dfdb3 dfdb4];

end





    
