function dict = monodict2D(x,n)
    % takes point (x,y) (or mx2 array of pts)
    % and returns row vector [1,x,y,xy,x^2,y^2,x^y,...]
    % or m x (n+2 choose n) matrix given array
    % where n is the highest total order of monomial 
    
    index = 1;
    dict = zeros(length(x(:,1)), nchoosek(n+2,n));
    for i = 0:n
        for j = 0:i
            dict(:,index) = (x(:,1).^j).*(x(:,2).^(i-j));
            index = index + 1;
        end
    end

