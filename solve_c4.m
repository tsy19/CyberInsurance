function [ W , Y ] = solve_c4(gam,a,b,R,C,z,theta_o, eps_o, lambda, mu)
%this is solve condition 4: a<1/gam, b<1/gam, using KKT
%I checked the results compared with CVX tool, they are always true and way faster.
    if z==0
        W = 0; Y = 0;
        return;
    end
    
    p = gam*a;
    q = gam*b;
    u = exp(gam*z*C) - 1;
    v = eps_o*(exp(gam*z*(C+R)) - exp(gam*z*C));
    s = theta_o * u;
    r = theta_o * v;
    
    if mu*p ~= lambda*q
        tem1 = s*mu*(p-lambda)/(lambda*q - mu*p);
        tem2 = r*(mu*p-lambda*q)/s/lambda/q;
        if (tem1>=1 && tem2>=1)
            W = 1/ lambda * log(tem1);
            Y = 1/ mu * log(tem2);
            return;
        end
    end
    
    tem = (s+r)*(lambda - p)/p  ;
    if tem >=1
        Y = 0; W = 1/lambda * log(tem);
        lbd2 = q*exp( p*W +q*Y ) + s*q*exp( (p-lambda)*W + q*Y ) + r*(q-mu)*exp( (p-lambda)*W + (q-mu)*Y );
        if lbd2 > 0
            return;
        end
    end
    
    tem = r*(mu - q)/(1+s)/q;
    if tem >=1
        W = 0; Y = 1/mu*log(tem);
        lbd1 = p*exp(p*W +q*Y) + s*(p-lambda)*exp( (p-lambda)*W + q*Y ) + r*(p-lambda)*exp( (p-lambda)*W + (q-mu)*Y );
        if lbd1 > 0
            return;
        end
    end
    
    W = 0; Y = 0;
    lbd1 = p*exp(p*W +q*Y) + s*(p-lambda)*exp( (p-lambda)*W + q*Y ) + r*(p-lambda)*exp( (p-lambda)*W + (q-mu)*Y );
    lbd2 = q*exp( p*W +q*Y ) + s*q*exp( (p-lambda)*W + q*Y ) + r*(q-mu)*exp( (p-lambda)*W + (q-mu)*Y );
    if (lbd1>0 && lbd2>0)
        return;
    end
    
    W = -Inf;
    Y = -Inf;
end