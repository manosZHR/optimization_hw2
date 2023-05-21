function dfdb=complex_der(b,n)

    ep=1e-15;

    imag_1 = [complex(b(1),ep); b(2); b(3); b(4)];
    imag_2 = [b(1); complex(b(2),ep); b(3); b(4)];
    imag_3 = [b(1); b(2); complex(b(3),ep); b(4)];
    imag_4 = [b(1); b(2); b(3); complex(b(4),ep)];
    
    dfdb1 = imag(primal(imag_1,n,0))/ep;
    dfdb2 = imag(primal(imag_2,n,0))/ep;
    dfdb3 = imag(primal(imag_3,n,0))/ep;
    dfdb4 = imag(primal(imag_4,n,0))/ep;
    dfdb=[dfdb1;dfdb2;dfdb3;dfdb4];

end