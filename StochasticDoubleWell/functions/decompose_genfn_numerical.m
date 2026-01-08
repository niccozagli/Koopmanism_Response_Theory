function [hnum,vnum] = decompose_genfn_numerical(x,y,obsnum,Kmax,W)
% Compute expansion coefficients of generic function wrt:
% dictionary basis (h) and basis of Koopman eigenfunctions(v)

    [xgrid,ygrid] = meshgrid(x,y);

    kn = nchoosek(Kmax+2,Kmax);

    S = zeros(kn,1);
    G_leb = zeros(kn);
    % Gram matrix of dictionary functions wrt Lebesgue (integrated over
    % regular grid). This is Identity for orthonormal dictionary   

    index = 1;
    for i = 0:Kmax
        for j = 0:i
            
            dicnum = (xgrid.^j).*(ygrid.^(i-j));
            integrandnum = obsnum .* conj( dicnum );
            S(index) = trapz(y, trapz(x,integrandnum,2) );

            index2 = 1; 
            
            for k = 0:Kmax
                for l = 0:k
                    
                    dicnum2 = (xgrid.^l).*(ygrid.^(k-l));
                    
                    integrand = dicnum2 .* conj( dicnum );
                    G_leb(index2,index) = trapz(y, trapz(x,integrand,2));
                    
                    index2 = index2 + 1;
    
                end
            end
            

            index = index + 1;
        end
    end
    
    hnum = (G_leb.')\S;   
    vnum = W'*hnum;

end

