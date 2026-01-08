function Kcorr = KoopmanCorrs_singleterm(v1,v2,lambda,matrix,index1,n)
    
Kcorr = ( lambda(index1).^n .* v1(index1) ).' * matrix(index1,(2:end)) * conj(v2(2:end));

end