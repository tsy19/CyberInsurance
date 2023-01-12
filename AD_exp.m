function [W_star,Y_star,R_star,udefender,uattacker,Pay] = AD_exp(gam, C, I, eps_o, mu, theta_o, lambda)

eps_ACC = 1e3; %accuracy of epsilon, using in linspace()
theta_ACC = 1e3; %accuracy of theta, using in linspace()

%C>=I case(by theorem 3.1,3.2)
% if C>=I
%     Y_star = 0;
%     if  (1-gam)*exp(gam*I)>1 %gam < gam_1
%         W_star = log( (1-gam)/gam * (exp(gam*I)-1) );
%     else
%         W_star = 0;
%     end
%     R_star = I;
%     Pay = 1;
%     theta = exp(-W_star);
%     udefender = -(1-theta)*exp(gam*W_star) - theta*exp(gam*(I+W_star));
%     uattacker = theta*R_star;
%     return;
% end

%C<I

%initialization
eps_ls = linspace(eps_o*exp(-mu*I),eps_o,eps_ACC); 
%eps_ls = eps_ls(2:end); %eps = sig * exp(-Y); cannot be 0
theta_ls = linspace(theta_o*exp(-lambda*I),theta_o,theta_ACC);
%theta_ls = theta_ls(2:end); %theta = exp(-W); cannot be 0
[ee,tt] = meshgrid(eps_ls,theta_ls); %broadcasting
e_len = length(eps_ls);
t_len = length(theta_ls);

%get pay/recover list(pay:1,recover:0) and ransom list

%following 2 lines solves region2 in theorem 3.1
pay_ls = (eps_ls >= ( (exp(gam*(I-C))-1) / (exp(gam*I)-1) ) ); 
R_ls = pay_ls * I;
%follwing solves region3
r3_idx = sum(eps_ls < ( (exp(gam*(I-C))-1) / (exp(gam*I)-1) ));
R_pay = C + 1/gam * log( (1-eps_ls(1:r3_idx) )./(1-eps_ls(1:r3_idx) *exp(gam*C)) );
pay_ls(1:r3_idx) = (R_pay >= (eps_ls(1:r3_idx) * I));
R_ls(1:r3_idx) =  pay_ls(1:r3_idx).* R_pay + (1-pay_ls(1:r3_idx))*I ;

%generate cost function of the defender
R_mtx = repmat(R_ls,t_len,1);
pay_mtx = repmat(pay_ls,t_len,1);
ww = 1/lambda * (log(theta_o) - log(tt));
yy = 1/mu * (log(eps_o) - log (ee));
cd_pay = (1 - tt).* exp( gam*(ww+yy) ) + tt.*( exp(gam*(ww+yy+R_mtx)) );
cd_recover = (1 - tt).* exp( gam*(ww+yy) ) + tt.*( (1-ee).*exp(gam*(ww+yy+C)) + ee.*exp(gam*(ww+yy+C+I)) );
cd_mtx = cd_pay.*pay_mtx + cd_recover.*(1-pay_mtx);

%find optimal value
udefender = - min(cd_mtx,[],"all");
[idx_e,idx_t] = find(cd_mtx == -udefender);
Pay = pay_mtx(idx_e,idx_t);
W_star = ww(idx_e,idx_t);
Y_star = yy(idx_e,idx_t);
R_star = R_mtx(idx_e,idx_t);
theta = tt(idx_e,idx_t);
eps = ee(idx_e,idx_t);
uattacker = theta*(Pay*R_star+(1-Pay)*eps*R_star);


end





