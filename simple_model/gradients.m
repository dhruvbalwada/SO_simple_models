% estimate the first and second derivatives of the quantities. 
% center differencing in the interior
% one sided differencing at the edges

function [gradients] = gradients(M, dl)
    n = length(M)
    d1(1) = (M(2)-M(1))/dl;
    d1(n) = (M(end)-M(end-1))/dl;
    
    for i=2:n-1
        d1(i) = (M(i+1) - M(i-1))/2/dl;
    end
    
    d2(1) = (d1(2)-d1(1))/dl;
    d2(n) = (d1(end)-d1(end-1))/dl;
    
    for i=2:n-1
        d2(i) = (d1(i+1) - d1(i-1))/2/dl;
    end
        
    gradients.d1 = d1;
    gradients.d2 = d2;
end