function y = integrate2D(x0,num_tsteps,dv,eps,dt)
    y = zeros(num_tsteps+1,2);
    
    % we take initial cond. x0 as (x,y)
    y_current = x0;
    
    y(1,:) = y_current;  

    for i = 2:num_tsteps+1
        y(i,:) = y_current - dv(y_current).*dt + eps.*sqrt(dt).*randn(1,2);
        y_current = y(i,:);
    end
end
