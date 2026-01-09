function [h,v] = decompose_genfn_analytical(obs,Kmax,W)
% Compute expansion coefficients of generic function wrt:
% dictionary basis (h) and basis of Koopman eigenfunctions(v)

    kn = nchoosek(Kmax+2,Kmax);

    S = zeros(kn,1);
    G_leb = zeros(kn);
    % Gram matrix of dictionary functions wrt Lebesgue (integrated over
    % regular grid). This is Identity for orthonormal dictionary

    index = 1;
    for i = 0:Kmax
        for j = 0:i
            dic = @(x,y) (x.^j).*(y.^(i-j));
            integrand = @(x,y) obs(x,y) .* conj( dic(x,y) );
            S(index) = integral2(integrand,-2,2,-2,2);
            index2 = 1; 
            for k = 0:Kmax
                for l = 0:k
                    dic2 = @(x,y) (x.^l).*(y.^(k-l));
                    
                    integrand = @(x,y) dic2(x,y) .* conj( dic(x,y) );
                    G_leb(index2,index) = integral2(integrand,-2,2,-2,2);
                    
                    index2 = index2 + 1;
                end
            end
            index = index + 1;
        end
    end

    h = (G_leb.')\S;
    v = W'*h;

end