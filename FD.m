function dfdb=FD(b,n)
    
    ep=1e-5;
    e1=[ep 0 0 0];
    e2=[0 ep 0 0];
    e3=[0 0 ep 0];
    e4=[0 0 0 ep];
    
    
    dfdb1 = (primal(b+e1,n,0)-primal(b-e1,n,0))/(2*ep);
    dfdb2 = (primal(b+e2,n,0)-primal(b-e2,n,0))/(2*ep);
    dfdb3 = (primal(b+e3,n,0)-primal(b-e3,n,0))/(2*ep);
    dfdb4 = (primal(b+e4,n,0)-primal(b-e4,n,0))/(2*ep);
    dfdb=[dfdb1;dfdb2;dfdb3;dfdb4];

end