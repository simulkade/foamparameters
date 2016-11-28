using JuMP, Roots, PyPlot, Ipopt, DataFrames, NLopt, PyCall, Optim, JLD
@pyimport scipy.optimize as so
PyPlot.rc("font", family="sans")

# read the data table
data1 = readtable("AOS0.03_N2_Benth_CT_mu_foam.csv")
deleterows!(data1,7)

# trim the data
# Foam quality
ug = float(data1[:ug]) # [m/s]
uw = float(data1[:uw]) # [m/s]
# Viscosity data [Pa.s]
muw = float(data1[:mu_liq])
mug = float(data1[:mu_gas])
u = ug+uw # [m/s]
fg_exp = ug./u # [-]
# Foam apparent viscosity [Pa.s]
muf_exp = float(data1[:muf_tot]) # [Pa.s]
sw_ct= float(data1[:sw_ave])
n_data = length(ug)
# filter the low quality data for the low quality regime
ind_mu_max= indmax(muf_exp)
fg_trans= fg_exp[ind_mu_max]
ind_fg_low= find(fg_exp.<=fg_trans)
fg_low= fg_exp[ind_fg_low]
sw_low= sw_ct[ind_fg_low]
muf_low= muf_exp[ind_fg_low]
muw_low= muw[ind_fg_low]
# filter the high quality data for the high quality regime
ind_mu_max= indmax(muf_exp)
fg_trans= fg_exp[ind_mu_max]
ind_fg_high= find(fg_exp.>=fg_trans)
fg_high= fg_exp[ind_fg_high]
sw_high= sw_ct[ind_fg_high]
muf_high= muf_exp[ind_fg_high]
println("Data is loaded.")

# IFT value [N/m]
sigma_wg = 0.03;

# I only use the less-scattered data in the low quality region
#filter for measured saturation
ind1= ~isnan(sw_low)
krw_exp= (1.0-fg_low[ind1]).*muw_low[ind1]./muf_low[ind1]
sw_low_ct=sw_low[ind1]
n_low=length(sw_low_ct)
semilogy(sw_low_ct, krw_exp, "o")
m1 = Model(solver=IpoptSolver(print_level=1))
#m = Model()
# define the variables of the foam model
@defVar(m1, 0.05<=nw_fit<=6.0)
@defVar(m1, 0.05<=krw0_fit<=1.0)
@defVar(m1, 0.001<=swc_fit<=0.2)
@defVar(m1, 0.001<=sgr_fit<=0.1)
@setNLObjective(m1, Min, sum{(log(krw_exp[i])-log(krw0_fit)-nw_fit*log((sw_low_ct[i]-swc_fit)/(1.0-swc_fit-sgr_fit)))^2, i=1:n_low})
#@setNLObjective(m1, Min, sum{(krw_exp[i]-krw0_fit*((sw_low_ct[i]-swc_fit)/(1.0-swc_fit-sgr_fit))^nw_fit)^2, i=1:n_low})
solve(m1)
# get the result
liq_relperm=[getValue(krw0_fit), getValue(nw_fit), getValue(swc_fit), getValue(sgr_fit)]
sw1=[linspace(liq_relperm[3], 0.3, 50); linspace(0.3,1.0-liq_relperm[4], 50)]
semilogy(sw1, liq_relperm[1]*((sw1-liq_relperm[3])/(1-liq_relperm[4]-liq_relperm[3])).^liq_relperm[2], "--")
xlabel("Liquid saturation")
ylabel("Liquid relperm")
#df_krw=DataFrame(sw=sw_low_ct, krw=krw_exp)
#writetable("krw_low_AOS.csv", df_krw)
println(liq_relperm)

# Relperm data
# swc = 0.207
swc=0.07
sgr = 0.03
krg0 = 0.587
ng = 0.938
krw0 = 0.713
nw = 2.460
sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+
(sw.>=1-sgr).*ones(size(sw)))
kr(sw)=((sw.>=swc).*krg0.*(1-sws(sw)).^ng+(sw.<swc).*(1+(krg0-1)/swc*sw))
krw(sw)=((sw.<=1-sgr).*krw0.*sws(sw).^nw+(sw.>1-sgr).*
(-(1-krw0)/sgr.*(1.0-sw)+1.0))
