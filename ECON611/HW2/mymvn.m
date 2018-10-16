function pdf = mymvn(x,sigma)
    p   = size(x,1);
    
    r   = chol(sigma);
    
    pdf = (2*pi)^(-p/2)*det(r)*exp(-0.5*x*inv(sigma)*x');
    
end