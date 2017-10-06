% function [cf,a1f,b1f,a2f,b2f,ss]=reg_pcws(x,y)
function [cf,line1, line2]=reg_pcws2(x,y,wt,deg)

% pcws_reg computes the piecewise linear regression ---with two
% pieces--- to (x,y), for any possible change point, chooses the one leading
% to the smallest least-square error, and returns and plots the
% corresponding parameters.

% a1=inf(1,size(x,2)-1);
% a2=inf(1,size(x,2)-1);
% b1=inf(1,size(x,2)-1);
% b2=inf(1,size(x,2)-1);
ss=Inf(1,size(x,2));

for c=2:(size(x,2)-1)
    x1=x(1:c);
    y1=y(1:c);
        
    if(nargin<3)
        deg(1) = 1;
        deg(2) = 1;
    end
    
    linearCoef1 = polyfit(x1,y1,deg(1));
    linearFit1{c} = polyval(linearCoef1,x1);
    err1 = sum((linearFit1{c} - y1).^2);
    
    x2=x(c:end);
    y2=y(c:end);
    
    linearCoef2 = polyfit(x2,y2,deg(2));
    linearFit2{c} = polyval(linearCoef2,x2);
    
    err2 = sum((linearFit2{c} - y2).^2);
    
    if(nargin<3)
        wt(1) = 1;
        wt(2) = 1;
    end
    
    ss(c)=sum((wt(1)*err1) + (wt(2)*err2));
end

[~,cf]=min(ss);
line1 = linearFit1{cf};
line2 = linearFit2{cf};