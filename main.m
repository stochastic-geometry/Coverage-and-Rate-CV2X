clear;

d = 40;  % OBSERVATION WINDOW RADIUS

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

% ~~~ SHADOWING PARAMETERS ~~~~ %
ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
ln_sig1 = 4; % Std deviation of log-normal shadowing gain in dB for TIER 1

ln_mu20 = 0; % Mean of log-normal shadowing gain in dB for TIER 2 TYP LINE
ln_sig20 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 2 TYP LINE

ln_mu21 = 0;% Mean of log-normal shadowing gain in dB for TIER 2 OTHER LINES
ln_sig21 = 4;% Std deviation  of log-normal shadowing gain in dB for TIER 2 OTHER LINES

% % % ======================================================================================

num_iterations = 1e4;

sir_threshold_dB = 0;
sir_threshold = 10^(.1*sir_threshold_dB);

EVAL_RATE = 0;
bw = 10; % BW in MHz
rate_threshold = 10; % in Mbps
%-------------------------------------------------------------------
%======================Simulation=============================

sir_measured = zeros(1, num_iterations);
load_measured = zeros(1, num_iterations);


simulate_veh_hetnet;
sir_measured = sir;
load_measured = loadm;

coverage_prob = nnz( sir_measured > sir_threshold)/num_iterations;
if EVAL_RATE
    rate = bw./load_measured.*log2(1 + sir_measured);
    rate_coverage = nnz( rate > rate_threshold)/num_iterations;
end

%-----------------------------------------------------------------------
% % HERE's the code to compute this for various thresholds
%-------------------------------------------------------------
% sir_threshold_dB = -50:5:50;
% sir_threshold = 10.^(.1*sir_threshold_dB);
% for k1=1:length(sir_threshold)
%     coverage_prob(k1) = nnz(sir_measured > sir_threshold(k1))/num_iterations;
% end

% rate_threshold = 0:20;
% for k1=1:length(rate_threshold)
%     rate_coverage(k1) = nnz( rate > rate_threshold(k1))/num_iterations;
% end



