sf_u0 = exp( 1/alpha*log(10)/10*ln_mu20 + (1/alpha*log(10)/10)^2*ln_sig20^2/2);
sf_a = exp( 2/alpha*log(10)/10*ln_mu21 + (2/alpha*log(10)/10)^2*ln_sig21^2/2 );
sf_c = exp( 2/alpha*log(10)/10*ln_mu1 + (2/alpha*log(10)/10)^2*ln_sig1^2/2 );

sh_lambda_2 = sf_u0*lambda_2;
sh_lambda_a = sf_a*(pi*lambda_l*lambda_2);
sh_lambda_1 = sf_c*lambda_1;

peu = sh_lambda_2*(1/(sh_lambda_1*((p2*b2)*mlg_2/(p1*b1))^(-2/alpha)))^.5*exp(sh_lambda_2^2/(sh_lambda_1*pi*((p2*b2)*mlg_2/(p1*b1))^(-2/alpha)))*erfc((sh_lambda_2^2/(sh_lambda_1*pi*((p2*b2)*mlg_2/(p1*b1))^(-2/alpha)))^.5 );
pec = 1- peu;

cdfru = @(ru) 1 - exp(-2*sh_lambda_2*ru);
pdfru = @(ru) 2*sh_lambda_2*exp(-2*sh_lambda_2*ru);

cdfrc = @(rc) 1 - exp(-sh_lambda_1*pi*rc.^2);
pdfrc = @(rc) 2*pi*sh_lambda_1*rc.*exp(-sh_lambda_1*pi*rc.^2);

pdfr_g_u = @(r) pdfru(r)/peu.*(1-cdfrc( ((p2*b2)*mlg_2/(p1*b1))^(-1/alpha)*r)) ; %+ 1/peu*( cdfru( ((p1*b1)/(p2*b2)/mlg_2)^(-1/alpha)*r) - cdfru(r) ).*pdfrc(r);

pdfr_g_c = @(r) pdfrc(r)/pec.*(1- cdfru( ((p1*b1)/(p2*b2)/mlg_2)^(-1/alpha)*r)) ;


lapI0_g_u = @(s,r) exp (-2*sh_lambda_2*integral(@(x) 1 - (1+ s*(p2)*mlg_2*x.^(-alpha)/m20).^(-m20) , r, Inf,'ArrayValued', true));
lapIu_g_u = @(s,r) exp(-2*pi*sh_lambda_a*integral( @(x) (1 - (1+s*(p2)*slg_2*x.^(-alpha)/m21).^(-m21)).*x, (mlg_2/slg_2)^(-1/alpha)*r*0, Inf,'ArrayValued', true));
% lapIu_g_u = @(s,r) exp(-2*pi*lambda_l*integral( @(y) 1 - exp(-2*lambda_2*integral(@(x) 1/(1+ (s*(p2*b2)*slg_2)^(-1)*(y.^2+x.^2).^(alpha/2) ),0, Inf,'ArrayValued',true)), 0, Inf,'ArrayValued', true) );
lapIc1_g_u = @(s,r) exp(-2*pi*qc*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*mlg_1)*x.^(-alpha)/m1).^(-m1)).*x, (mlg_2*(p2*b2)/(mlg_1*p1*b1))^(-1/alpha)*r, Inf,'ArrayValued', true));
lapIc2_g_u = @(s,r) exp(-2*pi*(1-qc)*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*slg_1)*x.^(-alpha)/m1).^(-m1)).*x, (mlg_2*(p2*b2)/(mlg_1*p1*b1))^(-1/alpha)*r, Inf,'ArrayValued', true));
lapI_g_u = @(s,r) lapI0_g_u(s,r).*lapIu_g_u(s,r).*lapIc1_g_u(s,r).*lapIc2_g_u(s,r);

lapI0_g_c = @(s,r) exp (-2*sh_lambda_2*integral(@(x) 1 - (1+ s*(p2)*mlg_2*x.^(-alpha)/m20).^(-m20) ,((p1*b1)/(p2*b2)/mlg_2)^(-1/alpha)*r, Inf,'ArrayValued', true));
lapIu_g_c = @(s,r) exp(-2*pi*sh_lambda_a*integral( @(x) (1 - (1+s*(p2)*slg_2*x.^(-alpha)/m21).^(-m21)).*x, (mlg_2/slg_2)^(-1/alpha)*r*0, Inf,'ArrayValued', true));
% lapIu_g_c = @(s,r) exp(-2*pi*lambda_l*integral( @(y) 1 - exp(-2*lambda_2*integral(@(x) 1/(1+ (s*(p2*b2)*slg_2)^(-1)*(y.^2+x.^2).^(alpha/2) ),0, Inf,'ArrayValued',true)), 0, Inf,'ArrayValued', true) );
lapIc1_g_c = @(s,r) exp(-2*pi*qc*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*mlg_1)*x.^(-alpha)/m1).^(-m1)).*x, r, Inf,'ArrayValued', true));
lapIc2_g_c = @(s,r) exp(-2*pi*(1-qc)*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*slg_1)*x.^(-alpha)/m1).^(-m1)).*x, r, Inf,'ArrayValued', true));
lapI_g_c = @(s,r) lapI0_g_c(s,r).*lapIu_g_c(s,r).*lapIc1_g_c(s,r).*lapIc2_g_c(s,r);

cov_prob_u = @(b) integral(@(r) lapI_g_u(m20*b*r.^(alpha)/(p2*mlg_2),r).*pdfr_g_u(r), 0, Inf, 'ArrayValued', true);
cov_prob_c = @(b) integral(@(r) lapI_g_c(m1*b*r.^(alpha)/(p1*mlg_1),r).*pdfr_g_c(r), 0, Inf, 'ArrayValued', true);

%% 

if EVAL_RATE_COVERAGE
    cdfr0 = @(r0) 1-exp(-sh_lambda_2*r0);
    pdfr0 = @(r0) sh_lambda_2*exp(-sh_lambda_2*r0);
    cdft0 = @(t0) 1-exp(-sh_lambda_1*pi*t0.^2);
    pdft0 = @(t0) 2*pi*sh_lambda_1*t0.*exp(-sh_lambda_1*pi*t0.^2);
    k = ((p2*b2)/(p1*b1))^(-1/alpha);
    pelu = integral(@(r0) sh_lambda_2*exp(-sh_lambda_1*pi*r0.^2*k^2/4 - sh_lambda_2*r0), 0, Inf,'ArrayValued', true);
    pelc = 1 - pelu;
    
    cdfz0_g_elu = @(z0) 1 - 1/pelu*integral(@(t0) (cdfr0(2*t0/k)-cdfr0(2*z0)).*pdft0(t0), k*z0, Inf,'ArrayValued', true);
    pdfz0_g_elu = @(z0) 2*pdfr0(2*z0)./pelu.*(1 - cdft0(k*z0));
    cdfz0_g_elc = @(z0) 1-1/pelc*integral(@(r0) (cdft0(k/2*r0) - cdft0(k*z0) ).*pdfr0(r0), 2*z0, Inf,'ArrayValued', true);
    pdfz0_g_elc = @(z0) k*pdft0(k*z0)./pelc.*(1-cdfr0(2*z0));
    
    k = ((p2*b2*mlg_2)/(p1*b1*mlg_1))^(-1/alpha);
    dl = (k+1)/(k-1);
    aa0 = @(z0) k*z0;
    aa1 = @(z1) k*z1;
    d = @(z1, z0) z0+z1;
    th = @(z1, z0) acos( ((z1+z0).^2+aa0(z0).^2-aa1(z1).^2)./(2*aa0(z0).*(z1+z0)));
    ath = @(z1, z0) ( ((z1+z0).^2+aa0(z0).^2-aa1(z1).^2)./(2*aa0(z0).*(z1+z0)));
    ph = @(z1, z0) acos( ((z1+z0).^2+aa1(z1).^2-aa0(z0).^2)./(2*aa1(z1).*(z1+z0)));
    aph = @(z1, z0) ( ((z1+z0).^2+aa1(z1).^2-aa0(z0).^2)./(2*aa1(z1).*(z1+z0)));
    
    a1 = @(z1, z0) pi*aa1(z1).^2 - (aa0(z0).^2.*(th(z1, z0)-.5*sin(2*th(z1, z0))) + aa1(z1).^2.*(ph(z1, z0)-.5*sin(2*ph(z1, z0))));
    a2 = @(z1, z0) pi*aa1(z1).^2 - pi*aa0(z0).^2;
    da1 = @(z1, z0) k^2*z1*(sin(2*acos((z0 + z1 - k^2*z0 + k^2*z1)/(2*k*z1))) - 2*acos((z0 + z1 - k^2*z0 + k^2*z1)/(2*k*z1))) + 2*pi*k^2*z1 - ((k^2 - 1)^2*(- k^2*z0^2 + 2*k^2*z0*z1 - k^2*z1^2 + z0^2 + 2*z0*z1 + z1^2))/(2*k*z0*(((k^2 - 1)*(- k^2*z0^2 + 2*k^2*z0*z1 - k^2*z1^2 + z0^2 + 2*z0*z1 + z1^2))/(k^2*z0^2))^(1/2)) + (z0*(k^2 - 1)^2*(- k^2*z0^2 + 2*k^2*z0*z1 - k^2*z1^2 + z0^2 + 2*z0*z1 + z1^2))/(2*k*z1^2*(((k^2 - 1)*(- k^2*z0^2 + 2*k^2*z0*z1 - k^2*z1^2 + z0^2 + 2*z0*z1 + z1^2))/(k^2*z1^2))^(1/2));
    da2 = @(z1, z0) 2*pi*k*aa1(z1);
    cdfz11_g_z0elu = @(z1, z0) 1 - exp( -2*sh_lambda_2*z1  ); % 0 < z1 < z0/dl
    cdfz12_g_z0elu = @(z1, z0) 1 - exp( -2*sh_lambda_2*z1 - sh_lambda_1*a1(z1,z0) ); % z0/dl < z1 < dl*z0
    cdfz13_g_z0elu = @(z1, z0) 1 - exp( -2*sh_lambda_2*z1 - sh_lambda_1*a2(z1,z0) ); % dl*z0 < z1 < Inf
    pdfz11_g_z0elu = @(z1, z0) exp( -2*sh_lambda_2*z1 ).*(2*sh_lambda_2);
    pdfz12_g_z0elu = @(z1, z0) exp( -2*sh_lambda_2*z1 - sh_lambda_1*a1(z1,z0) ).*(2*sh_lambda_2 + sh_lambda_1*da1(z1, z0));
    pdfz13_g_z0elu = @(z1, z0) exp( -2*sh_lambda_2*z1 - sh_lambda_1*a2(z1,z0) ).*(2*sh_lambda_2 + sh_lambda_1*da2(z1, z0));
    
    cdfz_1 = @(z) integral( @(z0) cdfz11_g_z0elu(z-z0, z0).*pdfz0_g_elu(z0), dl/(1+dl)*z,z,'ArrayValued',true);
    cdfz_2 = @(z) integral( @(z0) cdfz12_g_z0elu(z-z0, z0).*pdfz0_g_elu(z0), 1/(1+dl)*z,dl/(1+dl)*z,'ArrayValued',true);
    cdfz_3 = @(z) integral( @(z0) cdfz13_g_z0elu(z-z0, z0).*pdfz0_g_elu(z0), 0,1/(1+dl)*z,'ArrayValued',true);
    pdfz_1 = @(z) integral( @(z0) pdfz11_g_z0elu(z-z0, z0).*pdfz0_g_elu(z0), dl/(1+dl)*z,z,'ArrayValued',true) - dl/(1+dl)*cdfz11_g_z0elu(1/(1+dl)*z, dl/(1+dl)*z).*pdfz0_g_elu(dl/(1+dl)*z);
    pdfz_2 = @(z) integral( @(z0) pdfz12_g_z0elu(z-z0, z0).*pdfz0_g_elu(z0), 1/(1+dl)*z,dl/(1+dl)*z,'ArrayValued',true) + dl/(1+dl)*cdfz12_g_z0elu(1/(1+dl)*z, dl/(1+dl)*z).*pdfz0_g_elu(dl/(1+dl)*z) - 1/(1+dl)*cdfz12_g_z0elu(dl/(1+dl)*z, 1/(1+dl)*z).*pdfz0_g_elu(1/(1+dl)*z);
    pdfz_3 = @(z) integral( @(z0) pdfz13_g_z0elu(z-z0, z0).*pdfz0_g_elu(z0), 0,1/(1+dl)*z,'ArrayValued',true) + 1/(1+dl)*cdfz13_g_z0elu(dl/(1+dl)*z, 1/(1+dl)*z).*pdfz0_g_elu(1/(1+dl)*z);
    
    mz1 = integral(@(mz) mz.*pdfz_1(mz),0,Inf,'ArrayValued',true);
    mz2 = integral(@(mz) mz.*pdfz_2(mz),0,Inf,'ArrayValued',true);
    mz3 = integral(@(mz) mz.*pdfz_3(mz),0,Inf,'ArrayValued',true);
    meanz = real(mz1 + mz2 + mz3);
    
    pmf_typrsu_1  = @(k) integral( @(z) exp(-lambda_r*z).*(lambda_r*z).^(k)/factorial(k).*pdfz_1(z), 0, Inf,'ArrayValued', true);
    pmf_typrsu_2  = @(k) integral( @(z) exp(-lambda_r*z).*(lambda_r*z).^(k)/factorial(k).*pdfz_2(z), 0, Inf,'ArrayValued', true);
    pmf_typrsu_3  = @(k) integral( @(z) exp(-lambda_r*z).*(lambda_r*z).^(k)/factorial(k).*pdfz_3(z), 0, Inf,'ArrayValued', true);
    pmf_typrsu = @(k) pmf_typrsu_1(k) + pmf_typrsu_2(k) + pmf_typrsu_3(k);
    
    %%
    pdfw_1 = @(w) w.*pdfz_1(w)/meanz;
    pdfw_2 = @(w) w.*pdfz_2(w)/meanz;
    pdfw_3 = @(w) w.*pdfz_3(w)/meanz;
    
    pmf_taggedrsu_1  = @(k) integral( @(w) exp(-lambda_r*w).*(lambda_r*w).^(k)/factorial(k).*pdfw_1(w), 0, Inf,'ArrayValued', true);
    pmf_taggedrsu_2  = @(k) integral( @(w) exp(-lambda_r*w).*(lambda_r*w).^(k)/factorial(k).*pdfw_2(w), 0, Inf,'ArrayValued', true);
    pmf_taggedrsu_3  = @(k) integral( @(w) exp(-lambda_r*w).*(lambda_r*w).^(k)/factorial(k).*pdfw_3(w), 0, Inf,'ArrayValued', true);
    pmf_taggedrsu = @(k) pmf_taggedrsu_1(k) + pmf_taggedrsu_2(k) + pmf_taggedrsu_3(k);
    
    %%
    num_loadv_terms = 100; %max(tagged_RSU_Load)+1;
    pmfload_taggedRSU = zeros(1, num_loadv_terms);
    for kk1=1:num_loadv_terms
        pmfload_taggedRSU(kk1) = real(pmf_taggedrsu(kk1-1));
        if sum(pmfload_taggedRSU)>(1-1e-3)  break; end
    end
    num_loadv_terms = kk1;
    
    mean_taggedarea = 1.2857/sh_lambda_1;
    
    mean_cell_load2 = 1+lambda_r*(pi*lambda_l*mean_taggedarea + 3.216/pi/sh_lambda_1^.5)*(1- sh_lambda_2*meanz);
end


