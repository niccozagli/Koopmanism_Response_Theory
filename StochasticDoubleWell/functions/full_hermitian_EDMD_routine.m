function [Xi,W,lambda,G] = full_hermitian_EDMD_routine(y,Kmax)
    
    [G,A] = EDMD_2D_vectorised(y,Kmax);

    A_h = (A+A')/2;
    A = A_h;

    K = pinv(G)*A;
    [Xi,W,lambda] = get_spectral_properties(K);

end






