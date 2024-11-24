clear;
close all
clc

format short



%% set parameters
registered_anomalies_1 = [];
registered_anomalies_2 = [];

genotype = 'AD';
design = 4.1; % reinstatement (2 days extinction)
estimation  = 'hist'; % 'kde', 'hist', 'mdn'
file_path = '\\10.242.91.178\suzukiz\suzuki\FC\reinst\AD\12mo\';
% file_name = 'tbl_12mo_IQR_concatenate_2.csv';
file_name = 'mdn.csv';

parameters = set_parameters(design);

timestamp = char(datetime('now','Format','yyyyMMddHHmmss'));
str = ['parameters_',timestamp,'.mat'];
save([file_path,str],'parameters')

max_itr = 10;


%% load data
tbl = readtable([file_path,file_name],'delimiter',',');
tbl = sortrows(tbl,{'SN','Phase_','Trial','Day'});
[SN,ia,~] = unique(tbl.SN);
N = length(SN);
groups = tbl.Group(ia);
[G,~,ig] = unique(groups);




%% mcmc
T = parameters.T;
TF = parameters.TF;
k = find(TF); % index of estimated parameter
mu = parameters.mu;
w = parameters.w;
varnames = parameters.varnames;
M = parameters.M0;
L = parameters.L;
burnin = parameters.burnin;
m = parameters.m;
ub = parameters.ub;
lb = parameters.lb;
pd_p = parameters.pd_p;
Z = cell(1,N);
u = {parameters.u};

anomalies = false(1,N);
for n = 1:N

    fprintf('\n\n####### SN%d, %s #######\n',SN(n),[G{ig(n)}])
    if ismember(SN(n),registered_anomalies_1)
        fprintf('SN%d has been registered as an anomaly.\n\n',SN(n))
        continue
    else
        fprintf('detected anomalies: %s\n',strjoin(strsplit(num2str(SN(anomalies)'))))
    end

    % get observation
    y_obs = tbl.FreezingRate((tbl.SN==SN(n))&(strcmp(tbl.CS,'TRUE')));
    y_obs = (y_obs-min(y_obs))./(max(y_obs)-min(y_obs));

    t = 1:T;

    % create function handle
    pd = @(z_star) likelihood(z_star,y_obs,parameters);

    % Slice sampling until convergence
    Z_sample = [];
    itr = 1; % no. of sampling iteration
    while true
        disp(datetime)        
        fprintf('start %d-th itr...',itr)

        if itr==1
            z0 = mu(TF);
        else
            z0 = Z_sample_i(end,:);
        end
        tic
        [Z_sample_i, neval] = slicesample(z0,L,'burnin',burnin,'thin',M,'width',w(TF)','logpdf',pd);
        toc
        Z_sample = [Z_sample; Z_sample_i];
        
        R_hat = nan(1,length(TF)); % Gelman-Rubin statistic
        R_hat(TF) = GelmanRubin(size(Z_sample,1),m,Z_sample);
        z_hat = mu;
        
        % estimate parameters by histogram with optimal bin width
        for i = 1:sum(TF)
            [optN, ~, ~] = sshist_v2(sort(Z_sample(:,i)));
            edges = linspace(lb(k(i)),ub(k(i)),optN+1);
            [count,~] = histcounts(Z_sample(:,i),edges);
            [~,I] = max(count);
            bins = edges(1:end-1)+(edges(2)-edges(1))/2;
            z_hat(k(i)) = bins(I);
        end

        z_hat(3) = round(z_hat(3)); z_hat([5,11]) = ceil(z_hat([5,11]));

        R_hat = array2table(R_hat); R_hat.Properties.VariableNames = varnames;
        z_hat = array2table(z_hat); z_hat.Properties.VariableNames = varnames;
        ci = nan(2,11); ci(:,TF) = quantile(Z_sample,[0.025 0.975],1); ci = array2table(ci); ci.Properties.VariableNames = varnames;


        % simulate CR
        results = imm_run_tune(design, table2struct(z_hat));
        y_hat = results.V;
        % prediction errors
        err = (y_obs-y_hat).^2;
        p_err = arrayfun(@(x,y) cdf(x,y), pd_p, err); % cumulative density
        

        % show result
        fprintf('mean number of evaluations: %0.2f\n', neval)
        disp('Gelman-Rubin statistic')
        disp(R_hat)
        if sum(p_err>=u{2})==0; str = ''; else; str = strjoin(strsplit(num2str(find(p_err>=u{2})'))); end
        fprintf('prediction errors in %d trials exceed the threshold: %s\n',sum(p_err>=u{2}),str)
        figure;
        plot(1:T,y_obs,'o-',1:T,y_hat,'+-')
        drawnow
        
        if all(R_hat{1,TF}<u{1}) && all(p_err<u{2})
            disp('estimnated parameters')
            disp(z_hat)
            disp('credible interval of parameters')
            disp(ci)
            % save result
            Z{n} = Z_sample;
            break
        elseif size(Z_sample,1)>(L*max_itr)
            anomalies(n) = true;
            Z{n} = Z_sample;
            fprintf('SN%d is registered as an anomaly.\n',SN(n))
            break
        else
            itr = itr + 1;            
            disp('continue sampling...')
        end
    end       

end

str = ['Z_',timestamp,'.mat'];
save([file_path,str],'Z')




%% stats
Z_hat = nan(N,length(TF));
for n = 1:N
    
    if n>size(Z,2)
        break
    elseif isempty(Z{n})
        continue
    end

    Z_sample = Z{n};
    if strcmp(estimation,'mdn')
        Z_hat(n,TF) = median(Z_sample,1);
    elseif strcmp(estimation,'kde')        
        for i = 1:sum(TF)
            tin = linspace(lb(k(i)),ub(k(i)),1000);
            [y,t,optw,W,C,confb95,yb] = sskernel_ref_rate_v4(Z_sample(:,i),tin);
            [~,I] = max(y);
            x = linspace(lb(k(i)),ub(k(i)),length(y));
            Z_hat(n,k(i)) = x(I);
        end
    elseif strcmp(estimation,'hist')
        for i = 1:sum(TF)
            [optN, ~, ~] = sshist_v2(sort(Z_sample(:,i)));
            edges = linspace(lb(k(i)),ub(k(i)),optN+1);
            [count,~] = histcounts(Z_sample(:,i),edges);
            [~,I] = max(count);
            bins = edges(1:end-1)+(edges(2)-edges(1))/2;
            Z_hat(n,k(i)) = bins(I);
        end
    end
end

for i = 1:length(TF)
    if ~TF(i)
        continue
    end
    A = Z_hat(contains(groups,genotype),i);
    B = Z_hat(~contains(groups,genotype),i);
    p = ranksum(A,B);
    fprintf('%s, p = %0.2f, Mdn = %0.2f(%s), %0.2f(%s)\n',varnames{i},p,...
        median(A,'omitnan'),G{contains(G,genotype)},median(B,'omitnan'),G{~contains(G,genotype)})
end



%% save
estimation = [table(groups),table(SN),array2table(Z_hat)];
estimation.Properties.VariableNames = [{'Group','SN'},varnames];
str = ['estimation_',timestamp,'.csv'];
writetable(estimation,[file_path,str])



%% define likelihood
function p = likelihood(z_star,y_obs,parameters)

TF = parameters.TF;
mu = parameters.mu;
lb = parameters.lb;
ub = parameters.ub;
pd_p = parameters.pd_p;
design = parameters.design;
u = {parameters.u};

% simulate CR
z = zeros(1,length(TF));
z(TF) = z_star; z(~TF) = mu(~TF);

opts = struct('alpha',z(1),'g',z(2),'psi',round(z(3)),'eta',z(4),'maxIter',ceil(z(5)),...
    'w0',z(6),'sr',z(7),'sx',z(8),'theta',z(9),'lambda',z(10),'K',ceil(z(11)));

results = imm_run_tune(design, opts);
y_hat = results.V;


% prediction errors
err = (y_obs-y_hat).^2;
p_err = arrayfun(@(x,y) cdf(x,y), pd_p, err); % cumulative density

if any(isnan(err)) % inappropriate output of the model
    p = log(0);
elseif any(z_star<=lb(TF)')||any(z_star>(ub(TF)')) % inequality constraint
    p = log(0);
elseif any(p_err>=u{2})
    p = -1e6;
else
    p = log(sum(arrayfun(@(x,y) pdf(x,y), pd_p,err)));
end

end



