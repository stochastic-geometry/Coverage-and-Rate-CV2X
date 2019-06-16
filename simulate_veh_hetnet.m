


numc = poissrnd(lambda_1*pi*(d^2),1,num_iterations); % Number of TIER 1 nodes in observation window for all iterations
numl = poissrnd(lambda_l*(2*d*pi),1,num_iterations); % Number of LINES in observation window for all iterations


assocv = zeros(1, num_iterations);
assocv_otherlines = zeros(1, num_iterations);
assocc = zeros(1, num_iterations);
Ru = zeros(1, num_iterations);
Rc = zeros(1, num_iterations);
R = zeros(num_iterations,1);

tagged_BS_Load = ones(1, num_iterations);
tagged_RSU_Load = ones(1, num_iterations);
tagged_cell_length = zeros(1, num_iterations);
sir = zeros(1, num_iterations);
intfv0 = zeros(1, num_iterations);
intfv1 = zeros(1, num_iterations);
intfc = zeros(1, num_iterations);
drate = zeros(1, num_iterations);
loadm = zeros(1, num_iterations);
nnd = zeros(1, num_iterations);

for j1=1:num_iterations
    Nl = numl(j1); % Number of lines
    Nc = numc(j1); % Number of TIER 1 nodes
    
    
    % GENERATING THE RANDOM LINES
    rho1 = [0 -d+2*d*rand(1,Nl) ]; %[ (rv+(d-rv)*urv).*sign(urv-.5) 0];
    [~,ind2] = sort(abs(rho1));% - rho1(end);
    rho = rho1(ind2) - 2*d*(rho1(ind2)>d) + 2*d*(rho1(ind2)<-d);
    theta =pi*rand(1,length(rho));
    len_chord = 2*(d^2-(rho.^2)).^.5; % Length of chord of each line
    
    % ============ vehicular user generation ===================
    numv = poissrnd(lambda_r*len_chord);
    numv = numv + (numv==0);
    total_v = sum(numv);
    
    clens = cumsum(numv);
    idx=zeros(1,clens(end));idx2=idx; idx3=idx;
    idx([1 clens(1:end-1)+1]) = diff([0 len_chord]);
    len_c = cumsum(idx);
    idx2([1 clens(1:end-1)+1]) = diff([0 rho]);
    rho_vec = cumsum(idx2);
    idx3([1 clens(1:end-1)+1]) = diff([0 theta]);
    theta_vec = cumsum(idx3);
    
    x1 = (-len_c/2 + len_c.*rand(1,total_v));
    x = zeros(size(x1));
    for k1=1:length(rho)
        x( sum(numv(1:k1-1))+1: sum(numv(1:k1))) = sort(x1(sum(numv(1:k1-1))+1: sum(numv(1:k1))));
    end
    beta = atan(x./rho_vec);
    gamma = (theta_vec+beta);
    u = (rho_vec.^2+x.^2).^.5;
    signm = sign(rho_vec) + 1*(rho_vec==0);
    pts_v = ([u.*cos(gamma).*signm] + 1j*[u.*sin(gamma).*signm]).';
    
    
    % ============ TIER 2 generation ===================
    numu = poissrnd(lambda_2*len_chord);
    numu = numu + (numu==0);
    total_u = sum(numu);
    
    clens = cumsum(numu);
    idx=zeros(1,clens(end));idx2=idx; idx3=idx;
    idx([1 clens(1:end-1)+1]) = diff([0 len_chord]);
    len_c = cumsum(idx);
    idx2([1 clens(1:end-1)+1]) = diff([0 rho]);
    rho_vec = cumsum(idx2);
    idx3([1 clens(1:end-1)+1]) = diff([0 theta]);
    theta_vec = cumsum(idx3);
    
    x1 = [ -len_c/2 + len_c.*rand(1,total_u)];
    x = zeros(size(x1));
    for k1=1:length(rho)
        x( sum(numu(1:k1-1))+1: sum(numu(1:k1))) = sort(x1(sum(numu(1:k1-1))+1: sum(numu(1:k1))));
    end
    beta = atan(x./rho_vec);
    gamma = (theta_vec+beta);
    u = (rho_vec.^2+x.^2).^.5;
    signm = sign(rho_vec) + 1*(rho_vec==0);
    pts_u = ([ u.*cos(gamma).*signm] + 1j*[ u.*sin(gamma).*signm]).';
    dist_u = abs(pts_u);
    chnl_u = [gamrnd(m20, 1/m20, length(dist_u(1:numu(1))), 1); gamrnd(m21, 1/m21, length(dist_u(numu(1)+1:end)),1) ];
    
    % ================== Cellular BSs generation =======================
    rad_c = sqrt(rand(Nc,1))*d;
    phi_c = rand(Nc,1)*(2*pi);
    xc = [ rad_c.*cos(phi_c)];
    yc = [ rad_c.*sin(phi_c)]; % x =[ -d + 2*d*rand(1,Nc)];% y =[ -d + 2*d*rand(1,Nc)];
    pts_c = xc+1j*yc;
    dist_c = abs(pts_c);
    % Rc(j1) = min(dist_c);
    %  if Nc>0
    [vx,vy] = voronoi(xc,yc);
    [verts, cs] = voronoin([xc(:) yc(:)]);
    chnl_c = gamrnd(m1, 1/m1, length(dist_c), 1);
    %  end
    
    % ==================== Shadowing Generation ===========================
    shad_u0 = 10.^(.1*(ln_mu20 + ln_sig20*randn(numu(1),1)));
    shad_u = 10.^(.1*(ln_mu21 + ln_sig21*randn(sum(numu(2:end)),1)));
    shad_c = 10.^(.1*(ln_mu1 + ln_sig1*randn(length(dist_c), 1)));
    % ==================== ============== ============= ======== ========= ========== =
    new_ptsu = shad_u.^(-1/alpha).*pts_u(numu(1)+1:end);
    nnd_newptsu(j1) = min(abs(new_ptsu));
    
    
    % ==================== ============== ============= ======== ========= ========== =
    %     %     [val, ind] = min(dist_c);
    typ_verts_tmp = verts(cs{1},:);
    [ax,ay] =poly2cw(typ_verts_tmp(:,1) , typ_verts_tmp(:,2) );
    typ_verts =  [ax ay];
    typ_perm(j1) = sum(abs(diff( [typ_verts(:,1)+1j*typ_verts(:,2); typ_verts(1,1)+1j*typ_verts(1,2)])));
    typ_area(j1) = polyarea(typ_verts(:,1), typ_verts(:,2));
    % ==================== ============== ============= ======== ========= ========== =
    %     %     [val, ind] = min(dist_c);
    typ_verts_tmp = verts(cs{1},:);
    [ax,ay] =poly2cw(typ_verts_tmp(:,1) , typ_verts_tmp(:,2) );
    typ_verts =  [ax ay];
    typ_perm(j1) = sum(abs(diff( [typ_verts(:,1)+1j*typ_verts(:,2); typ_verts(1,1)+1j*typ_verts(1,2)])));
    typ_area(j1) = polyarea(typ_verts(:,1), typ_verts(:,2));
    
    Ru(j1) = min(dist_u(1:numu(1)));
    % ================== Load & SIR Computation =======================
    [val, ind] = max( [p2*mlg_2*b2*shad_u0.*dist_u(1:numu(1)).^(-alpha); p1*mlg_1*b1*shad_c.*dist_c.^(-alpha)] );
    if ind<= numu(1)
        if EVAL_RATE
            tagged_rsu = pts_u(ind);
            neigh_rsu = pts_u([find(real(pts_u(1:numu(1))-real(pts_u(ind)))<0,1,'first');find(real(pts_u(1:numu(1))-real(pts_u(ind)))>0,1,'last') ]); % pts_u([ind-1; ind+1]);%
            boundary_pts  = [ mean([neigh_rsu(1)]) ; mean([neigh_rsu(2)]) ];
            
            %         cell_end_pts  = [ mean([tagged_rsu neigh_rsu(1)]) ; mean([tagged_rsu neigh_rsu(2)]) ];
            pts_v0 = pts_v(1:numv(1));
            xdiff1 = sign(real(boundary_pts(1)) - real(pts_v0));
            xdiff2 = sign(real(boundary_pts(2)) - real(pts_v0));
            ydiff1 = sign(imag(boundary_pts(1)) - imag(pts_v0));
            ydiff2 = sign(imag(boundary_pts(2)) - imag(pts_v0));
            vind1 = find( ((xdiff1.*xdiff2)<=0) & ((ydiff1.*ydiff2)<=0) ) ;
            pts_v0_tagged1 = pts_v0(vind1);
            for kk1=1:length(pts_v0_tagged1)
                shad_from_rsu = (10.^(.1*(ln_mu20 + ln_sig20*randn(numu(1),1))) );
                shad_from_bs = 10.^(.1*(ln_mu1 + ln_sig1*randn(Nc, 1)));
                [maxpwr, maxpwr_ind] = max([p2*mlg_2*b2*shad_from_rsu.*abs(pts_v0_tagged1(kk1)-pts_u(1:numu(1))).^(-alpha); ...
                    p1*mlg_1*b1*shad_from_bs.*abs(pts_v0_tagged1(kk1)-pts_c).^(-alpha)]); %p2*mlg*shad_from_rsu(1)*abs(pts_v(kk1)-tagged_rsu)^(-alpha); ...
                tagged_RSU_Load(j1) = tagged_RSU_Load(j1) + (maxpwr_ind==ind) ;
            end
        end
        bf_gains = [mlg_2*ones(numu(1),1); slg_2*ones(total_u-numu(1),1)];
        shad_2 = [shad_u0; shad_u];
        sig_pwr = p2*bf_gains(ind)*shad_2(ind)*chnl_u(ind).*dist_u(ind).^(-alpha); %abs(tagged_rsu).^(-alpha);
        tmp2 = rand(length(dist_c),1)<qc;slg_1c = mlg_1*(tmp2==1) + slg_1*(tmp2==0);
        total_pwr = sum( [p2*mlg_2*chnl_u(1:numu(1)).*shad_u0.*dist_u(1:numu(1)).^(-alpha);...
            p2*slg_2*chnl_u(numu(1)+1:end).*shad_u.*dist_u( numu(1)+1:end).^(-alpha);  ...
            p1*slg_1c.*chnl_c.*shad_c.*dist_c.^(-alpha)] );
        sir(j1) = sig_pwr/(total_pwr - sig_pwr);
        ind_intfs = setdiff(1:numu(1),ind);
        intfv0(j1) = sum([ p2*mlg_2*chnl_u(ind_intfs).*shad_u0(ind_intfs).*dist_u(ind_intfs).^(-alpha)  ] ) ;
        intfv1(j1) = sum([  p2*slg_2*chnl_u(numu(1)+1:end).*shad_u.*dist_u( numu(1)+1:end).^(-alpha)] ) ;
        intfc(j1) = sum([p1*slg_1c.*chnl_c.*shad_c.*dist_c.^(-alpha)]);
        loadm(j1) = tagged_RSU_Load(j1);
        %v% drate(j1) = 1/tagged_RSU_Load(j1)*log2(1 + sir(j1));
        assocv(j1) = 1;
        assocv_otherlines(j1) = (ind>numu(1));
        R(j1) = dist_u(ind);
    else
        if EVAL_RATE
            ind_tagged_bs = ind-numu(1);
            tagged_bs = pts_c(ind_tagged_bs);
            croft_verts = verts(cs{ind-numu(1)},:);
            scale_factor = 1.5;
            indv1 = find(inpolygon(real(pts_v),imag(pts_v),croft_verts(:,1),croft_verts(:,2)) );
            %indv1 = find( abs(pts_v-tagged_bs)< scale_factor*max( abs([croft_verts(:,1)+1j*croft_verts(:,2)]-tagged_bs) ) );
            
            shad_from_bs = 10.^(.1*(ln_mu1 + ln_sig1*randn(Nc, 1)));
            for k1=1:length(indv1)
                line_ind = find( floor(cumsum(numv)./indv1(k1)) , 1 );
                shad_from_rsu = (10.^(.1*(ln_mu20 + ln_sig20*randn(numu(line_ind),1))) );
                % dist_nearest_rsu = min( abs(pts_u(sum(numu(1:line_ind-1))+1:sum(numu(1:line_ind)) ) - pts_v(indv1(k1)) ) );
                % num_pts_v_taggedcell(j1) = num_pts_v_taggedcell(j1)+( p2*dist_nearest_rsu^(-alpha) < p1*abs(pts_v(indv1(k1))-tagged_bs)^(-alpha) );
                [maxpwr, maxpwr_ind] = max([p1*mlg_1*b1*shad_from_bs(ind_tagged_bs).*abs(pts_v(indv1(k1))-tagged_bs).^(-alpha);...
                    p1*slg_1*b1*shad_from_bs(setdiff(1:Nc,ind_tagged_bs)).*abs(pts_v(indv1(k1))-pts_c(setdiff(1:Nc,ind_tagged_bs))).^(-alpha);...
                    p2*mlg_2*shad_from_rsu.*abs(pts_v(indv1(k1))-pts_u(sum(numu(1:line_ind-1))+[1:numu(line_ind)]) ).^(-alpha); ]);
                tagged_BS_Load(j1) = tagged_BS_Load(j1) + (maxpwr_ind==1) ;
            end
            
        end
        subind = numu(1);
        sig_pwr = p1*mlg_1*chnl_c(ind-subind)*shad_c(ind-subind)*dist_c(ind-subind).^(-alpha);%   abs(tagged_bs).^(-alpha);
        tmp2 = rand( length(dist_c),1)<qc;slg_1c = mlg_1*(tmp2==1) + slg_1*(tmp2==0); slg_1c(ind-subind) = mlg_1;
        total_pwr = sum( [p2*mlg_2*chnl_u(1:numu(1)).*shad_u0.*dist_u(1:numu(1)).^(-alpha); p2*slg_2*chnl_u(numu(1)+1:end).*shad_u.*dist_u( numu(1)+1:end).^(-alpha);  p1*slg_1c.*chnl_c.*shad_c.*dist_c.^(-alpha)] );
        sir(j1) = sig_pwr/(total_pwr - sig_pwr);
        intfv0(j1) = sum([ p2*mlg_2*chnl_u(1:numu(1)).*shad_u0.*dist_u(1:numu(1)).^(-alpha)  ] );
        intfv1(j1) = sum([  p2*slg_2*chnl_u(numu(1)+1:end).*shad_u.*dist_u( numu(1)+1:end).^(-alpha)] ) ;
        intfc(j1) = sum([p1*slg_1c.*chnl_c.*shad_c.*dist_c.^(-alpha)])-sig_pwr;
        loadm(j1) = tagged_BS_Load(j1);
        %v% drate(j1) = 1/tagged_BS_Load(j1)*log2(1 + sir(j1));
        assocc(j1) = 1;
        R(j1) = dist_c(ind-subind);
    end
end
prob_E1 = sum(assocc)/num_iterations;
prob_E2 = sum(assocv)/num_iterations;


% covv = zeros(size(threshold_dB));
% covc = covv;
% for k1=1:length(thr)
% % coverage_2 = nnz( sir(find(assocv))>sir_threshold )/max(1,nnz(assocv));
% % coverage_1 = nnz( sir(find(assocc))>sir_threshold )/max(1,nnz(assocc));
% end
% % coverage_probability = coverage_2*prob_E2 + coverage_1*prob_E1;
