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
% analyze the data
muf_max=max(muf_exp);

% fit a line to the high quality data; note: y=a+b*x
w1=ones(length(fg_high)+1) % weight factor
w1(end)=100.0
w1(end-3:end-1)=20.0
a_high=linreg((fg_high;1.0), (muf_high;0.0), w1)
%a_high(1)=-a_high(2) % force it to go through point (1,0)
% fit a power law to the low quality data
% y=a*x^b; log(y)=log(a)+b*log(x)
w2=ones(length(fg_low)) % weight factor
w2(end-3:end)=1.0 % a larger weight to the points near the transition
a_low=linreg(log(fg_low),log(muf_low), w2)
fg_low_range=linspace(0.0, maximum(fg_low)+0.2, 50)
fg_high_range=linspace(maximum(fg_low),1.0, 10)
plot(fg_low, muf_low, 'vr', fg_low_range, exp(a_low(1))*fg_low_range.^a_low(2), 'r-')
plot(fg_high, muf_high, 'b^', fg_high_range, a_high(1)+a_high(2)*fg_high_range, 'b-')
% find the transition quality
fun1(x)=exp(a_low(1))*x^a_low(2)-(a_high(1)+a_high(2)*x)
fg_trans = fzero(fun1, maximum(fg_low))
muf_trans = a_high(1)+a_high(2)*fg_trans
plot(fg_trans, muf_trans, '*y', markersize=15, (fg_trans, fg_trans), (0.0, 1.1*muf_trans), '--y', linewidth=1)
xlabel('foam quality')
ylabel('apparent viscosity (Pa.s)')
% find the transition saturation (or fmdry)
sw_trans = (1-swc-sgr)*(muw(ind_mu_max).*(1-fg_trans)./(krw0*muf_trans)).^(1/nw)+swc
fmdry_br= sw_trans
% find fmmob
FM_trans = mug(ind_mu_max)/muw(ind_mu_max)*fg_trans/(1-fg_trans)*(krw(sw_trans)/kr(sw_trans))
fmmob_br= 1.0/FM_trans-1.0
% find epdry
epdry_br= 100000.0 % very convenient
% visualize the three-parameter fit
x_br=(fmmob_br, epdry_br, fmdry_br)
sw_val = (linspace(0.0, minimum(sw_exp), 100); linspace(minimum(sw_exp)+eps(), maximum(sw_exp), 100); 
    linspace(maximum(sw_exp)+eps(), 1.0, 100))
fg_opt = fg(sw_val, x_br)
muf_opt = mu_foam(sw_val, x_br)
plot(fg_opt, muf_opt, '-.g', linewidth=2)
% calculate the epcap
dryout(x, sw)=0.5+atan(x(1)*(sw-x(2)))/π
fg_plus=0.7*fg_trans
muf_plus= exp(a_low(1))*fg_plus^a_low(2)
plot(fg_plus, muf_plus, 'k+', markersize=20)
sw_plus=(1-swc-sgr)*(muw(ind_mu_max).*(1-fg_plus)./(krw0*muf_plus)).^(1/nw)+swc
%FM_plus = mug(ind_mu_max)/muw(ind_mu_max)*fg_plus/(1-fg_plus)*(krw(sw_plus)/kr(sw_plus))
%epcap_br=log((1.0/FM_trans-1.0)/(1.0/FM_plus-1.0)*
%dryout((epdry_br,fmdry_br),sw_plus)/dryout((epdry_br,fmdry_br),sw_trans))/log(muf_plus/muf_trans)
epcap_br=log((kr(sw_trans)*muf_plus-fg_plus*mug(ind_mu_max))/
(fmmob_br*mug(ind_mu_max)*fg_plus))/log(muf_trans/muf_plus)
% calculate fmcap
ind_muf_min = indmin(muf_exp)
fmcap_br= muf_exp(ind_muf_min)*u(ind_muf_min)/sigma_wg
% correct for fmmob
fmmob_br_new = fmmob_br*(muf_trans/muf_exp(ind_muf_min))^epcap_br
% visualize the final results
x_br2= (fmmob_br_new, epdry_br, fmdry_br, fmcap_br, epcap_br)
sw_p = (linspace(0.0, minimum(sw_exp), 100); linspace(minimum(sw_exp)+eps(), maximum(sw_exp), 100); 
    linspace(maximum(sw_exp)+eps(), 1.0, 100))
n1 = length(sw_p)
muf_opt = zeros(n1)
u_min = minimum(u)
u_max = maximum(u)
u_ave = mean(u)
for i=1:n1
    fmu(x)= x-mu_foam2(sw_p(i), x_br2(1:3), x_br2(4:5), x, u_ave)
    muf_opt(i)=fzero(fmu, 0.9)
end
plot(fg2(sw_p, x_br2(1:3), x_br2(4:5), muf_opt, u_ave), mu_foam2(sw_p, x_br2(1:3), 
x_br2(4:5), muf_opt, u_ave), 'k--', linewidth=2)
axis((0.0,1.0,0.0,0.12))
fit_br=(fmmob_br_new, epdry_br, fmdry_br, fmcap_br, epcap_br)
