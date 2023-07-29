function DoFs = DoF(array)
for i=1:length(array)
    for j=1:length(array)
        diff(i,j)=array(i)-array(j);
    end
end
    diff=round(diff(:),6);
    [diff_u,~,~]=unique(diff);
    DoFs=length(diff_u);
end

