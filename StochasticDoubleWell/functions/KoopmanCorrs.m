function Kcorr = KoopmanCorrs(v1,v2,lambda,matrix)
    
Kcorr = @(n) ( lambda(2:end).^n .* v1(2:end) ).' * matrix(2:end,2:end) * conj( v2(2:end));

end
