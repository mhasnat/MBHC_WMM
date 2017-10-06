function [cf,line1, line2]=l_method(x,y)

numCLust = size(x,2);

ss=Inf(1,numCLust);

for c=2:(numCLust-1)
    x1=x(1:c);
    y1=y(1:c);
        
    linearCoef1 = polyfit(x1,y1,1);
    linearFit1{c} = polyval(linearCoef1,x1);
    err1 = sum((linearFit1{c} - y1).^2);
    
    x2=x(c:end);
    y2=y(c:end);
    
    linearCoef2 = polyfit(x2,y2,1);
    linearFit2{c} = polyval(linearCoef2,x2);
    err2 = sum((linearFit2{c} - y2).^2);
    
    wt(1) = (c-1)/(numCLust-1);
    wt(2) = (numCLust-c)/(numCLust-1);
    
    ss(c)=sum((wt(1)*err1) + (wt(2)*err2));
end

[~,cf]=min(ss);
line1 = linearFit1{cf};
line2 = linearFit2{cf};