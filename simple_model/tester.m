for j=2:length(b0)-1
    for i=1:length(l) 
        ylev= chary(i,j);
        idz1 = find(chary(:,j-1)<=ylev,1,'last'); 
        idz2 = find(chary(:,j+1)<=ylev,1,'last'); 
        if (~isempty(idz1) && ~isempty(idz2))
            db0dz(i,j) = (b0(j+1)-b0(j-1))/(charz(idz2,j+1)-charz(idz1,j-1)); 
        else
            db0dz(i,j) = NaN;
        end
    end
end