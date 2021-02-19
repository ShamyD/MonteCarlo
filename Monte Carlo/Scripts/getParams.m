function [mu_2,gamma_d, A_d] = getParams(cns)
    log_cns = log(cns);
    len_n = size(log_cns,2);

    lb = ceil(0.8*len_n);
    rb = size(log_cns,2);
    rep = size(log_cns,1);

    log_mu2 = zeros(rep,1);

    for i=1:rep
        ifit  = polyfit(lb:rb,log_cns(i,lb:rb),1);
        log_mu2(i) = ifit(1);
    end
    mu_2 = exp(log_mu2)
    m_fits = mean(mu_2);

    % form f_n = log(A_d) + (gamma_d - 1)log(n)
    rb2 = ceil(0.1*len_n);
    f_n = log(cns(:,1:rb2)) - repmat([1:rb2],rep,1).*log_mu2;

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

