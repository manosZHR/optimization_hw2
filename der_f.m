function dfdb=der_f(b,n,method)
    
    switch method
        %% FD
        case 1
            
            dfdb=FD(b,n);
    
        %% CV
        case 2
    
            dfdb=complex_der(b,n);
    
        %% Continuous Adjoint
        case 3
    
            dfdb=continuous_adj(b,n);
            dfdb=dfdb';
    
        %% Continuous DD
        case 4
    
            dfdb=continuous_DD(b,n);
            dfdb=dfdb';

        %% Discrete Adjoint
        case 5
    
            dfdb=discrete_adj(b,n);
            dfdb=dfdb';


    
    end %switch

end
