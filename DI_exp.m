function  [uinsurer,tinsurer,zinsurer,ainsurer,binsurer,pinsurer,wdefender,ydefender,uo,uattacker,Pay,proportion,total,Wo,Yo,Ro,ua_o,Pay_o]  = DI_exp(gam, C, I, theta_o, eps_o,lambda, mu)
    

    %run outside contract
    [Wo,Yo,Ro,uo,ua_o,Pay_o] = AD_exp(gam, C, I, eps_o, mu, theta_o, lambda);

    %run first problem
    [t1,a1,p1,u1,w1,y1,ua1,proportion1,total1] = first(I,Ro,gam,uo,C, theta_o, eps_o,lambda, mu);
    %ud1 = -(p1+a1*w1+t1*Ro);
    %theta = exp(-w1);
    %ud1 = (1-theta)*exp(gam*(a1*w1+p1)) + theta*( exp(gam*(a1*w1+p1+t1*Ro)) );

    %run second problem
    [z2,a2,b2,p2,u2,w2,y2,ua2,proportion2,total2] = second(I,Ro,gam,uo,C, theta_o, eps_o,lambda, mu);
    %eps = sig*exp(-y2);
    %ud2 = -(p2+a2*w2+b2*y2+z2*( eps*(C+Ro) + (1-eps)*C ) );
    %theta = exp(-w2);
    %ud2 = (1-theta)*exp(gam*(a2*w2+b2*y2+p2)) + theta*( (1-eps)*exp(gam*(a2*w2+b2*y2+p2+z2*C)) + eps*exp(gam*(a2*w2+b2*y2+p2+z2*(C+Ro))) );

    %compare
    if u1>=u2
        uinsurer = u1;
        tinsurer = t1;
        zinsurer = 1;
        ainsurer = a1;
        binsurer = 1;
        pinsurer = p1;
        wdefender = w1;
        ydefender = y1;
        udefender = uo;
        uattacker = ua1;
        proportion = proportion1;
        total = total1;
        Pay = 1;
    else
        uinsurer = u2;
        tinsurer = 1;
        zinsurer = z2;
        ainsurer = a2;
        binsurer = b2;
        pinsurer = p2;
        wdefender = w2;
        ydefender = y2;
        udefender = uo;
        uattacker = ua2;
        proportion = proportion2;
        total = total2;
        Pay = 0;
    end

end