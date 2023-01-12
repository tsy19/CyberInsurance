function [zinsurer,ainsurer,binsurer,pinsurer,uinsurer,wdefender,ydefender,uattacker,proportion,total] = second(I,R,gam,uo,C, theta_o, eps_o,lambda, mu)

    ACC=50;
    COUNT=0;

    alist = linspace(0,1,ACC+1);
    blist = linspace(0,1,ACC+1);
    zlist = linspace(0,1,ACC);
    ys = zeros( length(alist)-1,length(blist)-1 );
    ws = zeros( length(alist)-1,length(blist)-1 );
    us = zeros( length(alist)-1,length(blist)-1 );

    k=0;
    for z=zlist(1:end)
        k=k+1;

        i=0;
        for a=alist(2:end)
            i=i+1;
            j=0;
            for b=blist(2:end)
                j=j+1;

                p = gam*a;
                q = gam*b;
                u = exp(gam*z*C)-1;
                v = eps_o*(exp(gam*z*(C+R)) - exp(gam*z*C));
                s = theta_o * u;
                r = theta_o * v;
                
                
                if p >= lambda

                    if q>= r*mu/(1+s+r) %C1+C2_1
                        ys(i,j) = 0;
                        ws(i,j) = 0;
                        eps = eps_o * exp(-mu*ys(i,j));
                        theta = theta_o*exp(-lambda*ws(i,j));
                        
                        
                        %check z,a,b
                        if (1-eps)*exp(gam*z*C) + eps*exp(gam*z*(C+R)) >= exp(gam*R) %z,a,b not valid
                            us(i,j) = -Inf;
                        else %z,a,b valid
                            us(i,j) = 1/gam*log(-uo) - 1/gam *log( (1-theta)*exp(gam*a*ws(i,j)+gam*b*ys(i,j)) + theta*(1-eps)*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*C)) + theta*eps*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*(C+R))) ) -theta*(1-z)*(C+eps*R) - (1-a)*ws(i,j) - (1-b)*ys(i,j);
                        end 

                    else %C2_2
                        ys(i,j) = 1 / mu * log( r*(mu - q) / q / (1+s) );
                        ws(i,j) = 0;
                        eps = eps_o * exp(-mu*ys(i,j));
                        theta = theta_o*exp(-lambda*ws(i,j));

                        %check z,a,b
                        if (1-eps)*exp(gam*z*C) + eps*exp(gam*z*(C+R)) >= exp(gam*R) %z,a,b not valid
                            us(i,j) = -Inf;
                        else %z,a,b valid
                            us(i,j) = 1/gam*log(-uo) - 1/gam *log( (1-theta)*exp(gam*a*ws(i,j)+gam*b*ys(i,j)) + theta*(1-eps)*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*C)) + theta*eps*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*(C+R))) ) -theta*(1-z)*(C+eps*R) - (1-a)*ws(i,j) - (1-b)*ys(i,j);
                        end 

                    end
                else % p < lambda 
                    if q >= mu %C3
                        if p >= (s+r)*lambda / (1+s+r) %c3_1
                            ys(i,j) = 0;
                            ws(i,j) = 0;

                            eps = eps_o * exp(-mu*ys(i,j));
                            theta = theta_o*exp(-lambda*ws(i,j));

                            %check z,a,b
                            if (1-eps)*exp(gam*z*C) + eps*exp(gam*z*(C+R)) >= exp(gam*R) %z,a,b not valid
                                us(i,j) = -Inf;
                            else %z,a,b valid
                                us(i,j) = 1/gam*log(-uo) - 1/gam *log( (1-theta)*exp(gam*a*ws(i,j)+gam*b*ys(i,j)) + theta*(1-eps)*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*C)) + theta*eps*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*(C+R))) ) -theta*(1-z)*(C+eps*R) - (1-a)*ws(i,j) - (1-b)*ys(i,j);
                            end 
                            
                            
                        else %c3_2
                            ys(i,j) = 0;
                            ws(i,j) = 1/lambda * log( (r+s)*(lambda - p)/p );
                            
                            eps = eps_o * exp(-mu*ys(i,j));
                            theta = theta_o*exp(-lambda*ws(i,j));
                            %check z,a,b
                            if (1-eps)*exp(gam*z*C) + eps*exp(gam*z*(C+R)) >= exp(gam*R) %z,a,b not valid
                                us(i,j) = -Inf;
                            else %z,a,b valid
                                us(i,j) = 1/gam*log(-uo) - 1/gam *log( (1-theta)*exp(gam*a*ws(i,j)+gam*b*ys(i,j)) + theta*(1-eps)*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*C)) + theta*eps*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*(C+R))) ) -theta*(1-z)*(C+eps*R) - (1-a)*ws(i,j) - (1-b)*ys(i,j);
                            end 

                        end

                    else %c4
                        %COUNT = COUNT+1;
                        
                        %{
                        cvx_begin 
                            variables w y
                            minimize (exp(p*w + q*y) + s*exp((p-lambda)*w+q*y) + r*exp((p-lambda)*w + (q-mu)*y) )
                            subject to
                               w >= 0
                               y >= 0
                        cvx_end 
                        %}
                        [ W , Y ] = solve_c4(gam,a,b,R,C,z,theta_o, eps_o, lambda, mu); 
                        
                        ws(i,j) = W;
                        ys(i,j) = Y;
                        
                        
                        eps = eps_o * exp(-mu*ys(i,j));
                        theta = theta_o*exp(-lambda*ws(i,j));

                        %check z,a,b
                        if (1-eps)*exp(gam*z*C) + eps*exp(gam*z*(C+R)) >= exp(gam*R) %z,a,b not valid
                            us(i,j) = -Inf;
                        else %z,a,b valid
                            us(i,j) = 1/gam*log(-uo) - 1/gam *log( (1-theta)*exp(gam*a*ws(i,j)+gam*b*ys(i,j)) + theta*(1-eps)*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*C)) + theta*eps*exp(gam*(a*ws(i,j)+b*ys(i,j)+z*(C+R))) ) -theta*(1-z)*(C+eps*R) - (1-a)*ws(i,j) - (1-b)*ys(i,j);
                        end 
                    end
                end
                
                
                %check if W,Y correct
                %{
                cvx_begin 
                    variables w y
                    minimize (exp(p*w + q*y) + s*exp((p-lambda)*w+q*y) + r*exp((p-lambda)*w + (q-mu)*y) )
                    subject to
                        w >= 0
                        y >= 0
                cvx_end 
                W = ws(i,j); Y = ys(i,j);
                if abs(ws(i,j)- w ) > 1e-3 || abs(ys(i,j)- y ) > 1e-3
                    flag = 9;
                end
                %}
                
            end   
        end
        
        %get the whole ws,ys,us table, compute best choice condition on z
        uz(k) = max(us,[],"all");
        [idx_a, idx_b]=find(us==uz(k));
        if length(idx_a)>1 %not unique
            FLAG = 1;
        end

        wz(k) = ws(idx_a(1),idx_b(1)); %if not unique, find the one minimize a,b
        yz(k) = ys(idx_a(1),idx_b(1));
        az(k) = alist(idx_a(1)+1);
        bz(k) = blist(idx_b(1)+1);



    end


    uinsurer = max(uz);
    idx_z = find(uz==uinsurer);
    if length(idx_z)>1
        FLAG = 2;
    end
    zinsurer = zlist(idx_z(1)); %if not unique, find the one minimize z 
    ainsurer = az(idx_z(1));
    binsurer = bz(idx_z(1));
    wdefender = wz(idx_z(1));
    ydefender = yz(idx_z(1));
    uattacker = theta_o*exp(-lambda*wdefender)* eps_o *exp(-mu*ydefender) * R ;
    proportion = uinsurer / (uattacker+uinsurer);
    total = uattacker+uinsurer;
    eps = eps_o*exp(-mu*ydefender);
    theta = theta_o*exp(-lambda*wdefender);

    %compute p
    t1 = (1-theta)*exp(gam*ainsurer*wdefender+gam*binsurer*ydefender);
    t2 = theta*(1-eps)*exp(gam*(ainsurer*wdefender+binsurer*ydefender+zinsurer*C));
    t3 = theta*eps*exp(gam*(ainsurer*wdefender+binsurer*ydefender+zinsurer*(C+R)));
    pinsurer = 1/gam*log(-uo) - 1/gam *log( t1 + t2 + t3 );

end


%C1 
% w = 0;
% y = 0;

%C2
% b < 1/gam * (sig*exp(gam*z*R)-sig)/(sig*exp(gam*z*R)+1-sig);
% y = log( sig*(1-gam*b)*(exp(gam*z*R)-1)/(gam*b) );
% w = 0;

%C3
% a < 1/gam *(1- exp(-gam*z*C) / (sig*exp(gam*z*R) - sig +1));
% w = log( (1-gam*a)/(gam*a) * (exp(gam*z*C)*(sig*exp(gam*z*R) - sig +1)-1) );
% y = 0;

%C4
%{
cvx_begin 
    variables w y
    minimize (exp(p*w + q*y) + u*exp((p-1)*w+q*y) + v*exp((p-1)*w + (q-1)*y) )
    subject to
       w >= 0
       y >= 0
cvx_end 
%}
