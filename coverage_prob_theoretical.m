clear all;

% % ============ SPATIAL MODEL ========= %
lambda_l = 10/pi; % (LINE_DENSITY/PI); Line density - km^-1
lambda_r = 15; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)
lambda_2 = 4; %NUMBER OF TIER 2 NODES PER KM
lambda_1 = .5; %NUMBER OF TIER 1 NODES PER KM SQ.

% % ============ TRANSMIT PARAMETERS ========= %
p1_dbm = 43; % Transmit power of tier 1 node in dBm
p1 = 10^(.1*(p1_dbm-30));
p2_dbm = 23; % Transmit power of tier 2 node in dBm
p2 = 10^(.1*(p2_dbm-30));

mlg_2 = 1; % Mainlobe gain TIER 2 node
slg_2 = .01; % Sidelobe gain TIER 2 node

mlg_1 = 1; % Mainlobe gain TIER 1 node
slg_1 = .01; % Sidelobe gain TIER 1 node

b1_db = 0; % Selection bias of TIER 1 node in dB
b1 = 10^(.1*b1_db);
b2_db = 0; % Selection bias of TIER 2 node in dB
b2 = 10^(.1*b2_db);

qc = .05; % Probability with which the mainlobe of interfering TIER 1 node is directed towards the typical receiver


% % ============ PROPAGATION PARAMETERS =================================================== %
alpha = 4; % PATH_LOSS EXPONENT

% ~~~ NAKAGAMI-M FADING PARAMETERS ~~~ %
m1 = 1; % Tier 1
m20 = 2; % Tier 2 typical line
m21 = 1; % Tier 2 other lines
% m = 1;

% ~~~ SHADOWING PARAMETERS ~~~~ %
ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
ln_sig1 = 4; % Std deviation of log-normal shadowing gain in dB for TIER 1

ln_mu20 = 0; % Mean of log-normal shadowing gain in dB for TIER 2 TYP LINE
ln_sig20 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 2 TYP LINE

ln_mu21 = 0;% Mean of log-normal shadowing gain in dB for TIER 2 OTHER LINES
ln_sig21 = 4;% Std deviation  of log-normal shadowing gain in dB for TIER 2 OTHER LINES

% % % ======================================================================================

sir_threshold_dB = 0;%-50:5:50;
sir_threshold = 10.^(.1*sir_threshold_dB);

%---------------END of PARAMETERS ------------------------------------------------------
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

%
zeta21 = (mlg_2*(p2*b2)/(mlg_1*p1*b1))^(-1/alpha);
lapI0_g_u = @(s,r) exp (-2*sh_lambda_2*integral(@(x) 1 - (1+ s*(p2)*mlg_2*x.^(-alpha)/m20).^(-m20) , r, Inf,'ArrayValued', true));
lapIu_g_u = @(s,r) exp(-2*pi*sh_lambda_a*integral( @(x) (1 - (1+s*(p2)*slg_2*x.^(-alpha)/m21).^(-m21)).*x, 0, Inf,'ArrayValued', true));
% lapIu_g_u = @(s,r) exp(-2*pi*lambda_l*integral( @(y) 1 - exp(-2*lambda_2*integral(@(x) 1/(1+ (s*(p2*b2)*slg_2)^(-1)*(y.^2+x.^2).^(alpha/2) ),0, Inf,'ArrayValued',true)), 0, Inf,'ArrayValued', true) );
lapIc1_g_u = @(s,r) exp(-2*pi*qc*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*mlg_1)*x.^(-alpha)/m1).^(-m1)).*x, zeta21*r, Inf,'ArrayValued', true));
lapIc2_g_u = @(s,r) exp(-2*pi*(1-qc)*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*slg_1)*x.^(-alpha)/m1).^(-m1)).*x, zeta21*r, Inf,'ArrayValued', true));
lapI_g_u = @(s,r) lapI0_g_u(s,r).*lapIu_g_u(s,r).*lapIc1_g_u(s,r).*lapIc2_g_u(s,r);

lapI0_g_c = @(s,r) exp (-2*sh_lambda_2*integral(@(x) 1 - (1+ s*(p2)*mlg_2*x.^(-alpha)/m20).^(-m20), 1/zeta21*r, Inf,'ArrayValued', true));
lapIu_g_c = @(s,r) exp(-2*pi*sh_lambda_a*integral( @(x) (1 - (1+s*(p2)*slg_2*x.^(-alpha)/m21).^(-m21)).*x, 0, Inf,'ArrayValued', true));
% lapIu_g_c = @(s,r) exp(-2*pi*lambda_l*integral( @(y) 1 - exp(-2*lambda_2*integral(@(x) 1/(1+ (s*(p2*b2)*slg_2)^(-1)*(y.^2+x.^2).^(alpha/2) ),0, Inf,'ArrayValued',true)), 0, Inf,'ArrayValued', true) );
lapIc1_g_c = @(s,r) exp(-2*pi*qc*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*mlg_1)*x.^(-alpha)/m1).^(-m1)).*x, r, Inf,'ArrayValued', true));
lapIc2_g_c = @(s,r) exp(-2*pi*(1-qc)*sh_lambda_1*integral( @(x) (1 - (1+s*(p1*slg_1)*x.^(-alpha)/m1).^(-m1)).*x, r, Inf,'ArrayValued', true));
lapI_g_c = @(s,r) lapI0_g_c(s,r).*lapIu_g_c(s,r).*lapIc1_g_c(s,r).*lapIc2_g_c(s,r);

pc_e1 = @(b) 0;
f1 = @(s,x) -2*sh_lambda_2*(1 - (1+ s*(p2)*mlg_2*x.^(-alpha)/m20).^(-m20));
f2 = @(s,x) -2*pi*sh_lambda_a*(1 - (1+s*(p2)*slg_2*x.^(-alpha)/m21).^(-m21)).*x;
f3 = @(s,x) -2*pi*qc*sh_lambda_1*(1 - (1+s*(p1*mlg_1).*x.^(-alpha)/m1).^(-m1)).*x -2*pi*(1-qc)*sh_lambda_1*(1 - (1+s*(p1*slg_1)*x.^(-alpha)/m1).^(-m1)).*x;
fg_mf = @(r,s) integral(@(x) f1(s,x), r/zeta21, Inf,'ArrayValued', true) + integral(@(x) f2(s,x),0, Inf,'ArrayValued', true) + integral(@(x) f3(s,x),r,Inf,'ArrayValued', true);

if m1==1
    pc_e1 = @(b) (integral(@(r) lapI_g_c(m1*b*r.^(alpha)/(p1*mlg_1),r).*pdfr_g_c(r), 0, Inf, 'ArrayValued', true));
else
    syms s r x;
    assume(s,'real');
    assume(r,'positive');
    g1 = int( f1(s,x),x,r/zeta21,Inf);
    g2 = int( f2(s,x),x,0, Inf);
    g3 = int( f3(s,x),x,r,Inf);
    fg = g1+g2+g3;
    for n=0:m1-1
        derk1_I_e1 = @(r,s) 0;
        for k=0:n
            for j=0:k
                a1_str = char( diff( (fg)^(k-j), 's', n)); % To convert sym expr to function handle
                a1_str = strrep(a1_str, 'int(', 'integral(@(x)');
                a1_str = strrep(a1_str, 'x,','');
                a1_str = strrep(a1_str, '*','.*');
                a1_str = strrep(a1_str, '/','./');
                a1_str = strrep(a1_str, '^','.^');
                a1_mf = eval(['@(r,s) ' a1_str]); % Converting Sym expr to Matlab Function handle
                derk1_I_e1 = @(r,s) derk1_I_e1(r,s) + (-1)^j/factorial(k)*nchoosek(k,j)*a1_mf(r,s).*(fg_mf(r,s)).^j;
            end
        end
        dern_I_e1 = @(s,r) lapI_g_c(s,r).*derk1_I_e1(r,s);
        pc_e11 = @(b) (integral(@(r) (-m1*b*r.^(alpha)/(p1*mlg_1)).^n/factorial(n).*dern_I_e1(m1*b*r.^(alpha)/(p1*mlg_1),r).*pdfr_g_c(r), 0, Inf, 'ArrayValued', true));
        pc_e1 = @(b) pc_e1(b) +  pc_e11(b);
    end
    assume(s,'clear');
    assume(r,'clear');
    assume(x,'clear');
end

clear fg fg_mf a1_mf

pc_e2 = @(b) 0;
f1 = @(s,x) -2*sh_lambda_2*( 1 - (1+ s*(p2)*mlg_2*x.^(-alpha)/m20).^(-m20));
f2 = @(s,x) -2*pi*sh_lambda_a*(1 - (1+s*(p2)*slg_2*x.^(-alpha)/m21).^(-m21)).*x;
f3 = @(s,x) -2*pi*qc*sh_lambda_1*(1 - (1+s*(p1*mlg_1)*x.^(-alpha)/m1).^(-m1)).*x -2*pi*(1-qc)*sh_lambda_1*(1 - (1+s*(p1*slg_1)*x.^(-alpha)/m1).^(-m1)).*x;
fg_mf = @(r,s) integral(@(x) f1(s,x), r, Inf,'ArrayValued', true) + integral(@(x) f2(s,x),0, Inf,'ArrayValued', true) + integral(@(x) f3(s,x),zeta21*r,Inf,'ArrayValued', true);

if m20==1
    pc_e2 = @(b) (integral(@(r) lapI_g_u(m20*b*r.^(alpha)/(p2*mlg_2),r).*pdfr_g_u(r), 0, Inf, 'ArrayValued', true));
else
    syms s r x;
    assume(s,'real');
    assume(r,'positive');
    g1 = int( f1(s,x),x,r,Inf);
    g2 = int( f2(s,x),x,0, Inf);
    g3 = int( f3(s,x),x,zeta21*r,Inf);
    fg = g1+g2+g3;
    for n=0:m20-1
        derk1_I_e2 = @(r,s) 0;
        for k=0:n
            for j=0:k
                a1_str = char( diff( (fg)^(k-j), 's', n));
                a1_str = strrep(a1_str, 'int(', 'integral(@(x)');
                a1_str = strrep(a1_str, 'x,','');
                a1_str = strrep(a1_str, '*','.*');
                a1_str = strrep(a1_str, '/','./');
                a1_str = strrep(a1_str, '^','.^');
                a1_mf = eval(['@(r,s) ' a1_str]);
                derk1_I_e2 = @(r,s) derk1_I_e2(r,s) + (-1)^j/factorial(k)*nchoosek(k,j)*a1_mf(r,s).*(fg_mf(r,s)).^j;
            end
        end
        dern_I_e2 = @(s,r) lapI_g_u(s,r).*derk1_I_e2(r,s);
        pc_e22 = @(b)  (integral(@(r) (-m20*b*r.^(alpha)/(p2*mlg_2)).^n/factorial(n).*dern_I_e2(m20*b*r.^(alpha)/(p2*mlg_2),r).*pdfr_g_u(r), 0, Inf, 'ArrayValued', true));
        pc_e2 = @(b) pc_e2(b) + pc_e22(b);
    end
end


bval = sir_threshold;
cov_prob_u = pc_e2(bval);
cov_prob_c = pc_e1(bval);
thy_cov_prob = cov_prob_c*pec + cov_prob_u*peu;



