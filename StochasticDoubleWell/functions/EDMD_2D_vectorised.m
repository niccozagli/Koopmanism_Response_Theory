function [G,A] = EDMD_2D_vectorised(data,Kmax,batch_length)
    % EDMD with tensorised monomial dictionary 
    % Division by M now built-in, see end
  
    % If no batch length specified for vectorised routine, default to 10^7
    if nargin < 3
        batch_length = 10^7;
    end
    
    % Number of X-Y pairs
    M = length(data)-1;

    % Size of dict. array
    kn = nchoosek(Kmax+2,Kmax);
    
    G = zeros(kn,kn);
    A = G;

    num_batches = ceil(M/batch_length);

    fprintf('%d snapshot pairs of data loaded\n',M)
    fprintf('Performing vectorised EDMD using batches of size %d \n', batch_length)

    for m = 1:num_batches-1
        fprintf('chunk %d of %d \n',m,num_batches)
        psi_c = monodict2D(data((m-1)*batch_length+1:m*batch_length+1,:),Kmax);
        G = G + psi_c(1:batch_length,:)'*psi_c(1:batch_length,:);
        A = A + psi_c(1:batch_length,:)'*psi_c(2:batch_length+1,:);
    end
    
    psi_c = monodict2D(data((num_batches-1)*batch_length+1:end,:),Kmax);
    fprintf('final batch, no. %d,length %d \n',num_batches,length(psi_c(:,1))-1)
    G = G + psi_c(1:length(psi_c(:,1))-1,:)'*psi_c(1:length(psi_c(:,1))-1,:);
    A = A + psi_c(1:length(psi_c(:,1))-1,:)'*psi_c(2:end,:);

    G = G./M;
    A = A./M;
end    