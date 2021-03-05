function [candidate, pdf, revPdf] = rwp(t_prev, rho)
%generates a random walk proposal step candidate for an MH-sampler 
%t_prev is a row vector - Not matrix
    len_t = length(t_prev);
    candidate = zeros(1, len_t-2);
    candidate = [t_prev(1) candidate t_prev(len_t)];
    pdfCum = 1;
    
    %construct candidate and calulate the pdf
    for i = 2:len_t-1
        R = rho*(t_prev(i+1)-candidate(i-1));
        epsilon = 2*(rand-0.5)*R;
        candidate(i) = t_prev(i) + epsilon;
        pdfCum = pdfCum*2*R; 
    end
    pdf = 1/pdfCum;
    
    revPdfCum = 1;
    possible = 1;
    for i = 2:len_t-1
        R = rho*(candidate(i+1)-t_prev(i-1));
        
        %Check if within interval
        if abs(candidate(i)-t_prev(i)) > R
            possible = 0;
        end
        revPdfCum = revPdfCum*2*R;
    end
    revPdf = 1/revPdfCum*possible;
        
%     d1 = diag(-ones(len_t, 1));
%     d2 = diag(ones(len_t, 1), -2);
%     matr = d1(:,1:end-2) + d2(1:end-2,1:end-4);
%     R = t_prev*matr*rho;
%     
%     candidate = t_prev(2:end-1) + (rand(1,len_t-2)-0.5)*2.*R;
%     candidate = [t_prev(1),candidate, t_prev(end)];
%     candidate = sort(candidate, 2);
%     
%     pdf = 1/prod(2*R);
end