# Viscosity data [Pa.s]
muw = 0.6e-3
mug = 0.05e-3
# IFT value [N/m]
sigma_wg = 0.025
# Foam quality
fg_exp = [20, 30, 40, 50, 60, 70, 80, 90]/100.0
# Foam apparent viscosity [Pa.s]
muf_exp = [103.574871571, 133.536045851, 177.393262285,
    188.494578926, 211.520938156, 231.554481465, 224.811153811,
    150.598233034]/1000.0
# velocity [m/s]
u = 10*0.3048/(24*3600)
# Relperm data
swc = 0.10
sgr = 0.05
krg0 = 0.94
ng = 1.8
krw0 = 0.22
nw = 4.0
sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+
(sw.>=1-sgr).*ones(size(sw)))
kr(sw)=((sw.>=swc).*krg0.*(1-sws(sw)).^ng+(sw.<swc).*(1+(krg0-1)/swc*sw))
krw(sw)=((sw.<=1-sgr).*krw0.*sws(sw).^nw+(sw.>1-sgr).*
(-(1-krw0)/sgr.*(1.0-sw)+1.0))
sw_plot = collect(linspace(0.0,1.0, 100))
figure(1)
plot(sw_plot, krw(sw_plot), sw_plot, kr(sw_plot))
xlabel(L"S_w")
ylabel("RelPerms")
# foam model terms
fm(sw, F2)=1+F2[1]*(0.5+atan(F2[2].*(sw-F2[3]))/pi)
krg(sw, F2)=(kr(sw)./fm(sw, F2))

# Define the main functions
fg(sw, F2)=((krg(sw,F2)/mug)./(krw(sw)/muw+krg(sw,F2)/mug))
mu_foam(sw, F2)=(1./(krw(sw)/muw+krg(sw, F2)/mug))
fm2(sw, F2, F5, mu_f)=(1+F2[1]*(0.5+atan(F2[2].*(sw-F2[3]))/pi).*
(F5[1]./(mu_f*u/sigma_wg)).^F5[2])
krg2(sw, F2, F5, mu_f)=(kr(sw)./fm2(sw, F2, F5, mu_f))
fg2(sw, F2, F5, mu_f)=((krg2(sw,F2, F5, mu_f)/mug)./(krw(sw)/muw+
krg2(sw,F2, F5, mu_f)/mug));
mu_foam2(sw, F2, F5, mu_f)=(1./(krw(sw)/muw+krg2(sw, F2, F5, mu_f)/mug))
sw_exp = (1-swc-sgr)*(muw*(1-fg_exp)./(krw0*muf_exp)).^(1/nw)+swc
figure(2)
plot(fg_exp, sw_exp)
xlabel("Foam quality")
ylabel("Liquid saturation")
FM = mug/muw*fg_exp./(1-fg_exp).*(krw(sw_exp)./kr(sw_exp))
f_foam = 1./FM-1
m = Model()
# define the variables of the foam model
@defVar(m, 10<=fmmob<=1000000)
@defVar(m, 10<=epdry<=100000)
@defVar(m, swc<=fmdry<=swc+0.3)
@defVar(m, 1.0e-6<=fmcap<=1.0e-4)
@defVar(m, -2.0<=epcap<=4.0)
w=ones(length(muf_exp))
w[end-3:end-1]=10
@setNLObjective(m, Min, sum{w[i]*(f_foam[i]-fmmob*(0.5+atan(epdry*(sw_exp[i]-fmdry))/pi))^2, i=1:8})
solve(m)
# get the result
x1 = [getValue(fmmob), getValue(epdry), getValue(fmdry)]
@setNLObjective(m, Min, sum{w[i]*(f_foam[i]-
    fmmob*(0.5+atan(epdry*(sw_exp[i]-fmdry))/pi)
    *(fmcap/(u*muf_exp[i]/sigma_wg))^epcap)^2, i=1:8})
# solve it
solve(m)
# get the results
x2 = [getValue(fmmob), getValue(epdry), getValue(fmdry),
    getValue(fmcap), getValue(epcap)]
    n_data = length(fg_exp)
sw_plot = collect(linspace(0.0,1.0, 100))
n_plot = length(sw_plot)
muf_opt = zeros(n_plot)
for i=1:n_plot
	fmu(x)= x-mu_foam2(sw_plot[i], x2[1:3], x2[4:5], x)
	muf_opt[i]=fzero(fmu, 0.1)
end
figure(3)
plot(fg2(sw_plot, x2[1:3], x2[4:5], muf_opt), mu_foam2(sw_plot, x2[1:3],
x2[4:5], muf_opt), fg_exp, muf_exp, "o")
axis([0,1,0,0.25])
