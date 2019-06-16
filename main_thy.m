clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLEASE NOTE THAT THIS CODE IS FOR RAYLEIGH FADING%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ============ SPATIAL MODEL ========= %
lambda_l = 10/pi; % (LINE_DENSITY/PI); Line density - km^-1
lambda_r = 15; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)
lambda_2 = 4; %NUMBER OF TIER 2 NODES PER KM
lambda_1 = .5; %NUMBER OF TIER 1 NODES PER KM SQ.

% % ============ TRANSMIT PARAMETERS ========= %
p1_dbm = 40; % Transmit power of tier 1 node in dBm
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
m20 = 1; % Tier 2 typical line
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

bw = 10;
rate_threshold = 1;
EVAL_RATE_COVERAGE = 1;


eval_analy_expr;

% ---- COVERAGE PROB----------------------------------------------
for k1=1:length(sir_threshold)
    prob_1(k1) = pec;
    prob_2(k1) = peu;
    thy_cov2(k1) = cov_prob_u(sir_threshold(k1));
    thy_cov1(k1) = cov_prob_c(sir_threshold(k1));
    thy_cov_prob(k1) =  prob_2(k1)*thy_cov2(k1) + prob_1(k1)*thy_cov1(k1);
end
%---------------------------------------------------------------


if EVAL_RATE_COVERAGE
    % ========  RATE COVERAGE  ==============================================
    thy_rate_cov = zeros(size(rate_threshold));
    for k1=1:length(rate_threshold)
        % eval_analy_expr;
        ratecov_u = 0;
        for kk1=1:num_loadv_terms
            ratecov_u = ratecov_u + cov_prob_u(2^(rate_threshold(k1)/bw*kk1)-1)*real(pmfload_taggedRSU(kk1));
        end
        ratecov_c = cov_prob_c(2^(rate_threshold(k1)/bw*mean_cell_load2)-1);
        thy_rate_cov(k1) = peu*ratecov_u + pec*ratecov_c;
    end
    %=======================================================================
end


