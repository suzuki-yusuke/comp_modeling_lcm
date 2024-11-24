function parameters = set_parameters(design)


% experimental parameters
T = 36; % number of trials
test_index = false(T,1); test_index([4,34,36]) = true;

% MCMC parameters
l = 5e2; % sampling length in each chains [4e3]
m = 2; % number of chains [5]
M0 = 100; % sampling step [10]
L = l*m; % sample size finally obtained
burnin = floor(L*1); % set 20~100% of trials to calculate mean of parameters [floor(L*0.5)]
Sigma = ones(1,T).*0.3; % set sigma of distribution of prediction error
pd_p = []; for i = 1:T; pd_p = [pd_p; makedist('HalfNormal','mu',0,'sigma',Sigma(i))]; end % make distribution of prediction error

u = zeros(T,1)+0.5; u(1:3) = 1; u(test_index) = 0.25; % error threshold for each trial
u = {1.1, u}; % criteria of Gelman-Rubin statistic [1.1]

% LCM parameters
varnames = {'alpha','g','psi','eta','maxIter','w0','sr','sx','theta','lambda','K'}; % names of parameters
TF = logical([...
    1,... % concentration parameter alpha (default: 1)
    1,... % temporal scaling parameter g (default: 0)
    0,... % psi - [N x 1] binary vector specifying when protein synthesis inhibitor is injected (default: all 0)
    1,... % eta, learning rate (default: 1)
    1,... % maxIter, maximum number of iterations between each trial (default: 1)
    1,... % initial weight value w0 (default: 0)
    1,... % US variance sr (default: 0)
    1,... % stimulus variance sx (default: 0)
    1,... % response threshold theta (default: 1)
    1,... % response gain lambda (default: 1)
    1]); %  max number of latent sources K (default: 1)


% mean of multivariate Gaussian as a proposal dis tribution
mu = [...
    1.5,... % concentration parameter alpha (default: 0.1)
    1,... % temporal scaling parameter g (default: 1)
    0,... % psi - [N x 1] binary vector specifying when protein synthesis inhibitor is injected (default: all 0)
    0.1,... % eta, learning rate (default: 0.2)
    3,... % maxIter, maximum number of iterations between each trial (default: 3)
    0,... % initial weight value w0 (default: 1e-3)
    0.4,... % US variance sr (default: 0.4)
    1,... % stimulus variance sx (default: 1)
    0.01,... % response threshold theta (default: 0.008)
    0.01,... % response gain lambda (default: 0.007)
    33]; %  max number of latent sources K (default: 15)


% lower boundarlies
lb = [...
    1;... % concentration parameter alpha (default: 0.01,0.1,1)
    0.01;... % temporal scaling parameter g (default: 0.01)
    0;... % psi - [N x 1] binary vector specifying when protein synthesis inhibitor is injected (default: all 0)
    0.01;... % eta, learning rate (default: 0.01)
    1;... % maxIter, maximum number of iterations between each trial (default: 1)
    -0.01;... % initial weight value w0 (default: -1)
    0.01;... % US variance sr (default: 0.01)
    0.01;... % stimulus variance sx (default: 0.01)
    0.002;... % response threshold theta (default: 0.002)
    0.002;... % response gain lambda (default: 0.002)
    1]; %  max number of latent sources K (default: 1)


% upper boundaries
ub = [...
    3;... % concentration parameter alpha (default: 10)
    2;... % temporal scaling parameter g (default: 2)
    1;... % psi - [N x 1] binary vector specifying when protein synthesis inhibitor is injected (default: all 1)
    1;... % eta, learning rate (default: 0.5)
    5;... % maxIter, maximum number of iterations between each trial (default: 10)
    0.01;... % initial weight value w0 (default: 3)
    3;... % US variance sr (default: 3)
    3;... % stimulus variance sx (default: 3)
    0.02;... % response threshold theta (default: 0.01)
    0.02;... % response gain lambda (default: 0.01)
    T]; %  max number of latent sources K (default: 20)


% windows for slice sampling
w = (ub-lb)./50;


parameters = struct('design',design,'l',l,'m',m,'M0',M0,'L',L,'burnin',burnin,'Sigma',Sigma,'pd_p',pd_p,...
    'varnames',{varnames},'TF',TF,'mu',mu,'lb',lb,'ub',ub,'w',w,'test_index',test_index,'T',T,'u',u);

