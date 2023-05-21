function [F,u_new]=primal(b,nnx,flag)
    
    L=2;
    dx=(L/(nnx-1));
    
    x=linspace(0,L,nnx)';
    
    u=zeros(nnx,1);
    u(1)=b(1);
    u(2)=b(1)+dx*b(2);
    
    u_new=u;
    du=zeros(nnx,1);
    r=zeros(nnx,1);
    J=zeros(nnx,nnx);
    
    maxiter=1000;
    err=zeros(maxiter,1);
    
    iter=1;
    while true
        
        u_old=u_new;
        k=1/dx^2;
        for i=3:nnx
            a=k-2*x(i)+3*b(4)*u_old(i)^2;
            J(i,i)=a;
            J(i,i-1)=-2*k;
            J(i,i-2)=k;
            r(i)=-( b(4)*u_old(i)^3+k*u_old(i)-2*x(i)*u_old(i)+k*u_old(i-2)-2*k*u_old(i-1)+b(3)*(x(i)^2-x(i)) );
        end
        
    
        J(1,1) = 1;
        r(1)= - ( u(1)-b(1) );
    
        J(2,2) = 1;
        r(2) = -( u(2)-u(1)-b(2)*dx );
        
        J=sparse(J);
        du=J\r;
        u_new=u_old+du;
        
        err(iter)=norm(du);
        if err(iter)<1e-10
        break;end
    
        iter=iter+1;
        if iter>maxiter
        break;end
    
    
    end
    
    j=1;
    while x(j)<=0.5
        x_obj(j)=x(j);
        u_obj(j)=u_new(j);
        j=j+1;
    end
    

    y_obj=(u_obj-x_obj.^2).^2;
    F=1/2*trapz(x_obj,y_obj);

    if flag==1
        fprintf('Objective Function at %d',nnx');
        fprintf(' nodes: %4.7f \n',F);

        %     figure(1)
        %     plot(x,u_new)
        %     xlabel('x'); ylabel('u');
        %     grid on; box on; axis tight
        %     
        %     figure(2)
        %     plot(err,'-r')
        %     xlabel('||du||'); ylabel('#iteration');
        %     set(gca, 'YScale', 'log')
        %     grid on; box on; axis tight

    end
    

end
