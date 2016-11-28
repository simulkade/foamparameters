% Estimating foam parameters using foam scan data
% Written by A. A. Eftekhari
% Revised on September 4th, 2014
% Revised on September 24th, 2014: including a faster objective function,
% few more optimization algorithms, and a worse visualization of the
% results!
% Explanaton: This code first ignores the epcap, fmcap parameters and fits
% fmmob, fmdry, and epdry, which will be used as initial estimates for the
% next optimization step for all parameters
% note:
% F2(1) = fmmob
% F2(2) = epdry
% F2(3) = fmdry
% F5(1) = fmcap
% F5(2) = epcap
% some more clarification: Frist, I calculate the saturation values using
% the foam flood data (fg, and mu_app) and the liquid relperm curve. Then I
% use those saturation values in the FM term (stars model) and mobility
% terms to define a very efficient objective function.
clc; clear;
%% Physical properties
% muw = 9.5430e-04; % Pa.s
muw = 0.0009523072117647059; % Pa.s
% mug = 2e-5; % Pa.s
mug = 1.976701521987679e-5; % Pa.s
visc = [mug muw];
% A = pi()*0.038^2/4; % [m^2] cross section area of the core
% u = 1e-6/60/A; % [m/s]
% u = 10*0.3048/(24*3600); % [m/s] **** IT'S NEVER CONSTANT ****
% sigma_wg = 0.03; % [N/m]
%% data
fg_exp = [0.0492309
 0.0984484
 0.196418
 0.293757
 0.391969
 0.487949
 0.586442
 0.637661
 0.685237
 0.737687
 0.789534
 0.83579
 0.84349
 0.896027
 0.896229
 0.94964
 0.980033];

muf_exp = [0.469885
 0.546869
 0.678961
 0.791027
 0.869952
 1.02302
 1.21448
 1.28277
 1.37034
 1.46518
 1.4465
 1.07299
 1.15943
 0.896267
 0.900058
 0.299594
 0.156026]; % [Pa.s]
u = [1.46839e-5
 1.46705e-5
 1.46302e-5
 1.45658e-5
 1.45016e-5
 1.43499e-5
 1.4214e-5
 1.41953e-5
 1.40065e-5
 1.40059e-5
 1.3965e-5
 1.78987e-5
 1.40845e-5
 1.41342e-5
 1.41617e-5
 1.45908e-5
 1.47197e-5]; %
[fg_exp, fg_ind] = sort(fg_exp);
muf_exp = muf_exp(fg_ind);
u = u(fg_ind);
%% define the functions
% swc = 0.20;
% sgr = 0.05;
% krg0 = 0.94;
% ng = 4.8;
% krw0 = 0.12;
% nw = 4;
% relperm = [krg0 ng krw0 nw swc sgr];
relperm = [ 0.94    1.8    0.46453264685357337,1.9774829446424782,0.2,0.04968307320292199]; % foam ct 0.03%
% relperm = [ 0.9560    1.4293    0.2843    3.9173    0.068    0.01]; % nanofoam
krg0 = relperm(1);
ng = relperm(2);
krw0 = relperm(3);
nw = relperm(4);
swc = relperm(5);
sgr = relperm(6);
% swc = 0.10;
% sgr = 0.05;
% krg0 = 0.94;
% ng = 1.8;
% krw0 = 0.22;
% nw = 4;
sigma_wg = 0.03;

sws = @(sw)((sw>swc).*(sw<1-sgr).*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr).*ones(size(sw)));
kr = @(sw)((sw>=swc).*krg0.*(1-sws(sw)).^ng+(sw<swc).*(1+(krg0-1)/swc*sw));
krw = @(sw)((sw<=1-sgr).*krw0.*sws(sw).^nw+(sw>1-sgr).*(-(1-krw0)/sgr.*(1-sw)+1));

fm = @(sw, F2)(1+F2(1)*(0.5+atan(F2(2).*(sw-F2(3)))/pi()));
krg = @(sw, F2)(kr(sw)./fm(sw, F2));


%% Define the main functions
fg = @(sw, F2)((krg(sw,F2)/mug)./(krw(sw)/muw+krg(sw,F2)/mug));
mu_foam = @(sw, F2)(1./(krw(sw)/muw+krg(sw, F2)/mug));
fm2 = @(sw, F2, F5, mu_f)(1+F2(1)*(0.5+atan(F2(2).*(sw-F2(3)))/pi()).* ...
    (F5(1)./(mu_f.*u/sigma_wg)).^F5(2));
krg2 = @(sw, F2, F5, mu_f)(kr(sw)./fm2(sw, F2, F5, mu_f));
fg2 = @(sw, F2, F5, mu_f)((krg2(sw,F2, F5, mu_f)/mug)./(krw(sw)/muw+...
    krg2(sw,F2, F5, mu_f)/mug));
mu_foam2 = @(sw, F2, F5, mu_f)(1./(krw(sw)/muw+krg2(sw, F2, F5, mu_f)/mug));
%% find sw for a known liquid relperm
% sw = zeros(size(muf_exp));
% n_data = length(muf_exp);
% for i = 1:n_data
% %     fsw = @(sw)(krw(sw)-(muw*(1-fg_exp(i))./muf_exp(i)));
% %     sw(i) = fzero(fsw, 0.5);
%     sw(i) = (1-swc-sgr)*(muw*(1-fg_exp(i))./(krw0*muf_exp(i))).^(1/nw)+swc;
% end
sw = (1-swc-sgr)*(muw*(1-fg_exp)./(krw0*muf_exp)).^(1/nw)+swc;

plot(sw, fg_exp, 'o'); axis([0 1 0 1]);

%% calculate FM
FM = mug/muw*fg_exp./(1-fg_exp).*(krw(sw)./kr(sw));
FM2 = mug./muf_exp.*fg_exp./kr(sw);
f_foam = 1./FM-1;
k_dimless = muf_exp.*(krw(sw)/muw+kr(sw)/mug);
%% The initial three-parameter objective function
w = ones(size(fg_exp));

[mumax, muind] = max(muf_exp);
w(muind-1:muind+1)=5;
%     Fun5 = @(x)(w.*abs(fm(sw(2:end-1,i),x)-f_foam(2:end-1,i)));
Fun5 = @(x)(w.*abs(mu_foam(sw,x)-muf_exp));
x_guess5 = [1000 500 mean(sw)];
[x_new5, fval5]=lsqnonlin(Fun5, x_guess5,...
    [10 10 swc], [500000 100000 swc+0.6]);

sw_plot = linspace(0,1,500)';
figure(1); subplot(2,2,1);plot(fg(sw_plot, x_new5), mu_foam(sw_plot, ...
    x_new5), fg_exp, muf_exp, 'o');
title('Least square, no cap');
%% including capillary term; least square
options = optimset('TolX', 1e-30, 'TolFun', 1e-30, 'TolCon', 1e-30, ...
            'TolGradCon', 1e-30);
[mumax, muind] = max(muf_exp);
w(muind-1:muind+1)=10;
Fun = @(x)(w.*abs(mu_foam2(sw,x(1:3),x(4:5), muf_exp)-muf_exp));
% x_guess = [x_new5 1e-5 -0.7];
    x_guess = [x_guess5 1e-5 1];
[x_new_ls, fval_ls]=lsqnonlin(Fun, x_guess,...
    [10 10 swc 1e-6 -2], [500000 500000 swc+0.6 1e-3 5], options);
%% visualization; least square
mu_app2 = sw;
for j = 1:length(sw)
    options = optimset('TolX', 1e-30, 'TolFun', 1e-30, 'TolCon', 1e-30, ...
        'TolGradCon', 1e-30);
%         mu_app2(j) = fzero(@(x)(x-mu_foam2(sw(j), x_new_ls(1:3), x_new_ls(4:5), x)), 0.5);
    mu_app2(j) = fsolve(@(x)(x-mu_foam2(sw(j), x_new_ls(1:3), x_new_ls(4:5), x)), 0.3, options);
%         mu_app2(j,i) = fmincon(@(x)abs(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), 0.1, ...
%             [],[],[],[], 0, 10, [], options);
end

    figure(1);subplot(2,2,2);plot(fg2(sw, x_new_ls(1:3), x_new_ls(4:5), mu_app2), ...
        mu_foam2(sw, x_new_ls(1:3), x_new_ls(4:5), mu_app2), ...
        fg_exp, muf_exp, 'o');
    title('Least-square, all 5');

%% including capillary term; constrained optimization
[mumax, muind] = max(muf_exp);
w(muind-1:muind+1)=10;

Fun = @(x)sum(w.*abs(mu_foam2(sw,x(1:3),x(4:5), muf_exp)-muf_exp));
x_guess = x_new_ls;
%     x_guess = [x_guess5 1e-4 0.5];
[x_new_con, fval_con]=fmincon(Fun, x_guess, [],[],[],[],...
    [10 10 swc 1e-6 -2], [500000 100000 swc+0.6 1e-3 5]);


%% visualization; constrained
mu_app2 = sw;
for j = 1:length(sw)
    options = optimset('TolX', 1e-30, 'TolFun', 1e-30, 'TolCon', 1e-30, ...
        'TolGradCon', 1e-30);
%         mu_app2(j,i) = fzero(@(x)(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), [eps*100, 1]);
    mu_app2(j) = fsolve(@(x)(x-mu_foam2(sw(j), x_new_con(1:3), x_new_con(4:5), x)), 0.5, options);
%         mu_app2(j,i) = fmincon(@(x)abs(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), 0.1, ...
%             [],[],[],[], 0, 10, [], options);
end

figure(1);subplot(2,2,3);plot(fg2(sw, x_new_con(1:3), x_new_con(4:5), mu_app2), ...
    mu_foam2(sw, x_new_con(1:3), x_new_con(4:5), mu_app2), ...
    fg_exp, muf_exp, 'o');
title('sum of errors, all 5');

%% including capillary term; search based
[mumax, muind] = max(muf_exp);
w(muind-1:muind+1)=10;

Fun = @(x)sum(w.*abs(mu_foam2(sw,x(1:3),x(4:5), muf_exp)-muf_exp));
x_guess = x_new_con;
%     x_guess = [x_guess5 1e-4 0.5];
[x_new_src, fval_src]=fminsearch(Fun, x_guess, options);


%% visualization; search based
mu_app2 = sw;
for j = 1:length(sw)
    options = optimset('TolX', 1e-30, 'TolFun', 1e-30, 'TolCon', 1e-30, ...
        'TolGradCon', 1e-30);
%         mu_app2(j,i) = fzero(@(x)(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), [eps*100, 1]);
    mu_app2(j) = fsolve(@(x)(x-mu_foam2(sw, x_new_src(1:3), x_new_src(4:5), x)), 0.5, options);
%         mu_app2(j,i) = fmincon(@(x)abs(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), 0.1, ...
%             [],[],[],[], 0, 10, [], options);
end

figure(1);subplot(2,2,4);plot(fg2(sw, x_new_src(1:3), x_new_src(4:5), mu_app2), ...
    mu_foam2(sw, x_new_src(1:3), x_new_src(4:5), mu_app2), ...
    fg_exp, muf_exp, 'o');
title('sum of errors, all 5, search');

%% including capillary term; genetic (needs tuning)
% options = gaoptimset('PopulationSize', 100, 'TolFun', 1e-20, ...
%     'Generations', 5000);
% for i =1:n_data
%     [mumax, muind] = max(muf_exp_all(2:end-1,i));
%     w(muind-1:muind+1)=1;
%
%     Fun = @(x)sum(w.*abs(mu_foam2(sw(2:end-1,i),x(1:3),x(4:5), muf_exp_all(2:end-1,i)) ...
%         -muf_exp_all(2:end-1,i)));
%     x_guess = x_new(:,i);
% %     x_guess = [x_guess5 1e-4 0.5];
%     [x_new(:,i), fval(i)] = ga(Fun,5,[],[],[],[], ...
%         [10 10 swc 1e-6 -2], [500000 100000 swc+0.6 1e-3 5], [], options);
% end
%
% %% visualization; genetic
% mu_app2 = muf_exp_all;
% for i =1:n_data
%     for j = 2:length(fg_exp)-1
%         options = optimset('TolX', 1e-30, 'TolFun', 1e-30, 'TolCon', 1e-30, ...
%             'TolGradCon', 1e-30);
% %         mu_app2(j,i) = fzero(@(x)(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), [eps*100, 1]);
%         mu_app2(j,i) = fsolve(@(x)(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), 0.5, options);
% %         mu_app2(j,i) = fmincon(@(x)abs(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), 0.1, ...
% %             [],[],[],[], 0, 10, [], options);
%     end
%
%     figure(5);subplot(2,4,i);plot(fg2(sw(:,i), x_new(1:3,i), x_new(4:5,i), mu_app2(:,i)), ...
%         mu_foam2(sw(:,i), x_new(1:3,i), x_new(4:5,i), mu_app2(:,i)), ...
%         fg_exp_all(:,i), muf_exp_all(:,i), 'o');
%     title(['k = ' num2str(k_all(i))]);
% end

%% including capillary; pattern search
% x_new_src = zeros(5, n_data);
% fval_src=zeros(1,n_data);
% for i =1:n_data
%     [mumax, muind] = max(muf_exp_all(2:end-1,i));
%     w(muind-1:muind+1)=1;
%
%     Fun = @(x)sum(w.*abs(mu_foam2(sw(2:end-1,i),x(1:3),x(4:5), muf_exp_all(2:end-1,i)) ...
%         -muf_exp_all(2:end-1,i)));
%     x_guess = x_new_ls(:,i);
% %     x_guess = [x_guess5 1e-4 0.5];
%     [x_new_src(:,i), fval_src(i)]=patternsearch(Fun, x_guess, [],[],[],[],...
%         [10 10 swc 1e-6 -2], [500000 100000 swc+0.6 1e-3 5]);
% end
%
% %% visualization; search based
% mu_app2 = muf_exp_all;
% for i =1:n_data
%     for j = 2:length(fg_exp)-1
%         options = optimset('TolX', 1e-30, 'TolFun', 1e-30, 'TolCon', 1e-30, ...
%             'TolGradCon', 1e-30);
% %         mu_app2(j,i) = fzero(@(x)(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), [eps*100, 1]);
%         mu_app2(j,i) = fsolve(@(x)(x-mu_foam2(sw(j,i), x_new_src(1:3,i), x_new_src(4:5,i), x)), 0.5, options);
% %         mu_app2(j,i) = fmincon(@(x)abs(x-mu_foam2(sw(j,i), x_new(1:3,i), x_new(4:5,i), x)), 0.1, ...
% %             [],[],[],[], 0, 10, [], options);
%     end
%
%     figure(4);subplot(2,4,i);plot(fg2(sw(:,i), x_new_src(1:3,i), x_new_src(4:5,i), mu_app2(:,i)), ...
%         mu_foam2(sw(:,i), x_new_src(1:3,i), x_new_src(4:5,i), mu_app2(:,i)), ...
%         fg_exp_all(:,i), muf_exp_all(:,i), 'o');
%     title(['k = ' num2str(k_all(i))]);
% end
