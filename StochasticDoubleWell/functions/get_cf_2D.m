function [cf,lags] = get_cf_2D(Y,obs1,obs2,upto)
    signal1 = obs1(Y(:,1), Y(:,2)) - mean(obs1(Y(:,1), Y(:,2)));
    signal2 = obs2(Y(:,1), Y(:,2)) - mean(obs2(Y(:,1), Y(:,2)));

    [cf,lags] = xcorr(signal1,signal2,upto,'biased');

end

