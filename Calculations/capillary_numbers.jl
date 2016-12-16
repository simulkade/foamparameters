# this scrits calculate the capillary numbers in our core flooding experiments
# it is calculated by k grad(p)/sigma_wg
# or
# Nc=H/L dP/dPc
using DataFrames, Plots
# constants
dpc = 1e3 # Pa
L=0.17 # m core length
D_core =0.038 # core diameter
H_over_L=0.2 # simensionless

# data files
data_files=["AOS0.03_N2_Benth_CT_mu_foam.csv",
            "AOS0.1_N2_Benth_CT_mu_foam.csv",
            "AOS0.5_N2_Benth_CT_mu_foam.csv"]
            # "Amph0.5_N2_Benth_CT_mu_foam.csv"
labels=["0.03 AOS", "0.1 AOS", "0.5 AOS", "0.5 Amph"]
markers=[:circle, :square, :utriangle, :diamond]


# read the data files
i=0
plot(size=(500,400), xtickfont = Plots.font(10, "Courier"), ytickfont=Plots.font(10, "Courier"),
  legendfont=Plots.font(10, "Courier"),  guidefont=Plots.font(12, "Courier"))
for data_file in data_files
  data1=readtable(data_file)
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

  i+=1
  scatter!(data1[:dp_tot][ind_fg_high]*1e5/dpc*H_over_L, data1[:sw_ave][ind_fg_high], marker=markers[i]
  , markersize=7, label="", color=:lightgray)
  scatter!(data1[:dp_tot][ind_fg_low]*1e5/dpc*H_over_L, data1[:sw_ave][ind_fg_low], marker=markers[i]
  , markersize=7, label=labels[i])
end
plot!(ylabel="Liquid saturation", xlabel="Capillary number", ylims=(0.1,0.3))

savefig("capillary_number_AOS.png")
# plot pressure drops
i=0
plot(size=(500,400), xtickfont = Plots.font(10, "Courier"), ytickfont=Plots.font(10, "Courier"),
  legendfont=Plots.font(10, "Courier"),  guidefont=Plots.font(12, "Courier"))
for data_file in data_files
  data1=readtable(data_file)
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

  i+=1
  scatter!(fg_high, data1[:dp_tot][ind_fg_high], marker=markers[i]
  , markersize=7, label="", color=:lightgray)
  scatter!(fg_low, data1[:dp_tot][ind_fg_low], marker=markers[i]
  , markersize=7, label=labels[i])
end
plot!(ylabel="Pressure drop [bar]", xlabel="Gas fractional flow")
savefig("pressure_drop_AOS.png")
