function [wk,xhi_k] = gauss(k,n_G)

if n_G == 1
    wk = 2;
    xhi_k = 0;
end

if n_G==2
    wk = 1;
    if k == 1
        xhi_k = -1 / sqrt(3);
    else
         xhi_k = 1 / sqrt(3);
    end
end

if n_G==3
    if k == 1
         wk = 8 / 9;
         xhi_k = 0;
    else
        wk = 5 / 9;
        if  k == 2
            xhi_k = sqrt(3/5);
        else
            xhi_k = -sqrt(3/5);
        end
    end
end