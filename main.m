%% include code package to path
init_path;
% compile the CPD mex file for the first time
cpd_make;

%% Read Source Point
RGB = imread('scissor1.bmp'); I = rgb2gray(RGB); [row,col] = find(edge(I)~=0); %[row,col] = find(I==0); %
X_full = [row, col];
temp = pcdownsample(pointCloud([X_full,zeros(size(X_full,1),1)]),'gridAverage',7);
X = temp.Location(:,1:2);
X = X + repmat([-200,100],size(X,1),1);
%plot(X(:,1),X(:,2),'*')

%% Read Template Point
RGB = imread('scissor2.bmp'); I = rgb2gray(RGB); [row,col] = find(edge(I)~=0); %[row,col] = find(I==0); %
Y_full = [row, col];
temp = pcdownsample(pointCloud([Y_full,zeros(size(Y_full,1),1)]),'gridAverage',7);
Y = temp.Location(:,1:2);
%plot(Y(:,1),Y(:,2),'*')

%% Parameter Setting and Run Registration
% Change RegisterMethod to switch registration methods
RegisterMethod = 'SPR';  % SPR, CPD, GLTP, TPS-RPM

switch RegisterMethod
    case 'SPR'
        SPR_opt.method = 'nonrigid';   %'nonrigid','nonrigid_lowrank'
        SPR_opt.viz = 1;   % disable visualization to speed up!
        SPR_opt.max_it = 150; SPR_opt.tol = -1;  % disable tolerance check, only max_it --> same iterations
        SPR_opt.outliers = 0;
        SPR_opt.knn = 20;
        SPR_opt.tau = 500;
        SPR_opt.beta = 2;
        SRP_opt.lambda = 3;
        SPR_opt.tau_annealing_handle = @(iter, max_it)  0.97^iter; 
        SPR_opt.lambda_annealing_handle = @(iter, max_it) 0.97^iter;
        [SPR_Transform, ~] = SPR_register(Y, X, SPR_opt); % CPD warp Y to X, fliped!
        X_warp = SPR_Transform.Y;
        
    case 'CPD'
        CPD_opt.method = 'nonrigid';   %'nonrigid','nonrigid_lowrank'
        CPD_opt.viz = 1;
        CPD_opt.max_it = 150; CPD_opt.tol = -1;  % disable tolerance check, only max_it
        CPD_opt.outliers = 0;
        CPD_opt.beta = 2;
        CPD_opt.lambda = 3;
        [CPD_Transform, ~] = cpd_register(Y, X, CPD_opt);
        X_warp = CPD_Transform.Y;
        
    case 'GLTP'
        GLTP_opt.method = 'nonrigid';   %'nonrigid','nonrigid_lowrank'
        GLTP_opt.viz = 1;
        GLTP_opt.max_it = 150; GLTP_opt.tol = -1;  % disable tolerance check, only max_it
        GLTP_opt.outliers = 0;
        GLTP_opt.knn = 20;
        GLTP_opt.tau = 500;
        GLTP_opt.beta = 2;
        GLTP_opt.lambda = 3;
        GLTP_opt.tau_annealing_handle = @(iter, max_it) 0.97^iter;
        GLTP_opt.lambda_annealing_handle = @(iter, max_it) 0.97^iter;
        [GLTP_Transform, ~] = GLTP_register(Y, X, GLTP_opt);
        X_warp = GLTP_Transform.Y;
        
    case 'TPS-RPM'
        TPS_opt.frac = 1;
        TPS_opt.T_init     = 0.5;
        TPS_opt.T_finalfac = 500;
        [~, ~, ~, X_warp] = cMIX(X, Y, TPS_opt.frac, TPS_opt.T_init, TPS_opt.T_finalfac, 0);
end


%% Plot
figure;
plot(X(:,1), X(:,2),'bo', Y(:,1), Y(:,2),'r*'); 
legend('X', 'Y');
xlim([min([X(:,1); Y(:,1)]) - 50, max([X(:,1); Y(:,1)])+50]);
ylim([min([X(:,2); Y(:,2)]) - 50, max([X(:,2); Y(:,2)])+50]);

figure
plot(X_warp(:,1), X_warp(:,2),'bo', Y(:,1), Y(:,2),'r*'); 
xlim([min([X_warp(:,1); Y(:,1)]) - 50, max([X_warp(:,1); Y(:,1)])+50]);
ylim([min([X_warp(:,2); Y(:,2)]) - 50, max([X_warp(:,2); Y(:,2)])+50]);
legend( [RegisterMethod, ' X warp'], 'X true');

