function [mu_d,gamma_d, A_d] = getParams(cns)
    log_cns = log(cns);
    len_n = size(log_cns,2);
    rep = size(log_cns,1);
    
    lb = ceil(0.8*len_n);
    rb = size(log_cns,2);

    log_mud = zeros(rep,1);

    for i=1:rep
        ifit  = polyfit(lb:rb,log_cns(i,lb:rb),1);
        log_mud(i) = ifit(1);
    end
    mu_d = exp(log_mud)
    m_fits = mean(mu_d);

    % form f_n = log(A_d) + (gamma_d - 1)log(n)
    rb2 = max(ceil(0.1*len_n),10);
    f_n = log(cns(:,1:rb2)) - repmat([1:rb2],rep,1).*log_mud;

    gamma_d = zeros(rep,1);
    log_Ad = zeros(rep,1);
    for i=1:rep
        ifit  = polyfit(log(1:rb2),f_n(i,:),1);
        gamma_d(i) = ifit(1) + 1;
        log_Ad(i) = ifit(2);
    end

    gamma_d
    A_d = exp(log_Ad)

end

