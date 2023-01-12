function [tinsurer,ainsurer,pinsurer,uinsurer,wdefender,ydefender,uattacker,proportion,total] = first(I,R,gam,uo,C, theta_o, eps_o,lambda, mu)
% proportion is uinsurer / (uattacker+uinsurer);


    ACC = 50; % maybe we can reduce it to 50 to make the code run quicker? the result seems the same.
    tmax = min ( 1/(gam*R) * log( (1-eps_o)*exp(gam * C) + eps_o*exp( gam*(C+R) ) ) , 1); %t:tau
    tlist = linspace(0,tmax,ACC+1);
    j=0;
    for t = tlist
        j = j+1;
        amax = lambda/gam  * (exp(gam*t*R)-1)*theta_o / ( 1 + (exp(gam*t*R)-1)*theta_o );
        alist = linspace(0,1,ACC);
        i=0;
        for a = alist(2:end)
            i = i+1;
            
            if a < amax
                %{
                cvx_begin 
                    variables w2(i)
                    minimize (exp(gam*a*w2(i)) + (exp(gam*t*R)-1)*exp((gam*a-1)*w2(i)) )
                    subject to
                        w2(i) >= 0
                cvx_end 
                %}

                w(i) = 1/ lambda * log( (lambda-gam*a)*(exp(gam*t*R)-1)*theta_o / gam / a  );
                
                theta = theta_o*exp(-lambda*w(i));
                
                u2(i) = 1/gam*log(-uo) - 1/gam*log( 1 - theta + theta*exp(gam*t*R) ) - w(i) - theta*(1-t)*R;
                u(i) = 1/gam*log(-uo) - 1/gam * log( (1 - theta) * exp(gam*a*w(i)) + theta * exp( gam*(a*w(i)+t*R ) ) ) - theta*(1-t)*R - (1-a)*w(i);
                %u(i) = 1/gam*log(-uo) + (1-gam)/gam * log(1 - gam*a) + log(gam * a) -log(exp(gam*t*R)-1) - (  gam*a*(1-t)*R  /   ((1-gam*a)*(exp(gam*t*R)-1))    );

            else %a >= amax
                w(i) = 0;
                theta = theta_o*exp(-lambda*w(i));
                
                u2(i) = 1/gam*log(-uo) - 1/gam*log( 1 - theta + theta*exp(gam*t*R) ) - w(i) - theta*(1-t)*R;
                u(i) = 1/gam*log(-uo) - 1/gam * log( (1 - theta) * exp(gam*a*w(i)) + theta * exp( gam*(a*w(i)+t*R ) ) ) - theta*(1-t)*R - (1-a)*w(i);
                
            end
            %{
            www=w(i);
            
            cvx_begin
                variables w2(i)
                minimize (exp(gam*a*w2(i)) + (exp(gam*t*R)-1)*exp((gam*a-1)*w2(i)) )
                subject to
                    w2(i)>=0
            cvx_end
            WWW = w2(i);
            %}
        end
        
        %find optimal solution, condition on tau
        us(j) = max(u);
        idx_w=find(u==us(j));
        if length(idx_w)>1%not unique
            FLAG=1;
        end
        ws(j)= w(idx_w(1)); %if not unique, choose the one that has minimum a
        as(j)= alist(1+idx_w(1));
        %{
        if mod(j,100) ==0
            plot( alist(2:end-1) , u);
            hold on
        end
        %}
    end

%{
    figure;
    plot(tlist(2:end), ws);
    hold on 
    plot(tlist(2:end), us);
    hold on
    plot(tlist(2:end), ua);
    hold on
    plot(tlist(2:end), as);
%}
    
    uinsurer=max(us);
    idx_t= find(us==uinsurer);
    if length(idx_t)>1
        FLAG = 0;
    end
    tinsurer = tlist(idx_t(1));% if not unique, choose the one that minimize tau
    wdefender = ws(idx_t(1));
    ainsurer = as(idx_t(1));
    ydefender = 0;
    
    theta = theta_o*exp(-lambda*wdefender);
    uattacker = theta*R;
    proportion = uinsurer / (uattacker+uinsurer);
    total = uattacker+uinsurer;
    pinsurer = 1/gam*log(-uo) - 1/gam * log( (1 - theta) * exp(gam*ainsurer*wdefender) + theta * exp( gam*(ainsurer*wdefender+tinsurer*R ) ) );

end