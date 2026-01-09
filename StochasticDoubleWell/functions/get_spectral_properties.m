function [XiSorted,WSorted,lambdaSorted] = get_spectral_properties(Kort)
     
    % Find the spectrum of the matrix
    [Xi,Lambda,W] = eig(Kort); lambda = diag(Lambda);

    % Normalise left eigenvectors
    L = W';
    prod = L*Xi;
    Lcor = bsxfun(@rdivide,L,diag(prod)); 
    W = Lcor';
    
    % Order Eigenvalues
    [~,IndSorted] = sort(abs(lambda),'descend');
    lambdaSorted = lambda(IndSorted);
    XiSorted = Xi(:,IndSorted);
    WSorted = W(:,IndSorted);

end    