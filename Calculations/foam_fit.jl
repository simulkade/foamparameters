# first try to use JuMP or other optimization methods 
#for optimizing foam model parameters in Julia
# written by AAE, TU Delft, October 1st, 2014
# using JuMP
# using PyPlot

# Experimental data
muf_exp = [103.574871571	133.536045851	177.393262285	188.494578926	211.520938156	231.554481465	224.811153811	150.598233034;
105.559077163	142.464971013	162.508418842	192.472894629	213.495239227	234.534091359	219.834132298	146;
50.9868203848	74.0065766016	96.0309285158	120.042787528	148.02635923	161.108579956	148.415937034	89.0878596992;
68.8479722145	90.8756256355	108.931566367	128.97171269	138.102029766	132.347503401	100.795002839	50.3958506662]/1000.0

fg_exp = [20	30	40	50	60	70	80	90]/100

# Viscosity data
muw = 0.6e-3
mug = 0.05e-3
# velocity
u = 10*0.3048/(24*3600)
# Relperm data
swc = 0.10
sgr = 0.05
krg0 = 0.94
ng = 1.8
krw0 = 0.22
nw = 4
# IFT value
sigma_wg = 0.025;
# Relperm functions (pay attention to all the `dots` and two different deffinitions
# for Real and Array types
sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+(sw.>=1-sgr).*ones(size(sw)))
kr(sw)=((sw.>=swc).*krg0.*(1-sws(sw)).^ng+(sw.<swc).*(1+(krg0-1)/swc*sw))
krw(sw)=((sw.<=1-sgr).*krw0.*sws(sw).^nw+(sw.>1-sgr).*(-(1-krw0)/sgr.*(1.0-sw)+1.0))

fm(sw, F2)=1+F2[1]*(0.5+atan(F2[2].*(sw-F2[3]))/pi)
krg(sw, F2)=(kr(sw)./fm(sw, F2))

# Define the main functions
fg(sw, F2)=((krg(sw,F2)/mug)./(krw(sw)/muw+krg(sw,F2)/mug))
mu_foam(sw, F2)=(1./(krw(sw)/muw+krg(sw, F2)/mug))
fm2(sw, F2, F5, mu_f)=(1+F2[1]*(0.5+atan(F2[2].*(sw-F2[3]))/pi).*(F5[1]./(mu_f*u/sigma_wg)).^F5[2])
krg2(sw, F2, F5, mu_f)=(kr(sw)./fm2(sw, F2, F5, mu_f))
fg2(sw, F2, F5, mu_f)=((krg2(sw,F2, F5, mu_f)/mug)./(krw(sw)/muw+krg2(sw,F2, F5, mu_f)/mug));
mu_foam2(sw, F2, F5, mu_f)=(1./(krw(sw)/muw+krg2(sw, F2, F5, mu_f)/mug))

# calculate liquid saturations
sw_exp = (1-swc-sgr)*(muw*(1-fg_exp)./(krw0*muf_exp)).^(1/nw)+swc

# plot saturations and gas fractional flow
plot(fg_exp', sw_exp')

# Now we calculate the experimental foam parameters
FM = mug/muw*fg_exp./(1-fg_exp).*(krw(sw_exp)./kr(sw_exp))
f_foam = 1./FM-1

# Objective function
obj_fun(x) = sum(abs(fm(sw_exp[1,:], x)-f_foam[1,:]))
x_guess = [100.0, 500.0, 0.3]
res1 = optimize(obj_fun, x_guess, method = :gradient_descent, grtol=1e-12, autodiff=true)

lb = [10, 10, swc]
ub = [500000, 100000, swc+0.6]

f1 = DifferentiableFunction(obj_fun)
res2 = fminbox(f1, x_guess, lb, ub)

obj_fun2(x) = sum((mu_foam(sw_exp[1,:],x)-muf_exp[1,:]).^2)
res3 = optimize(obj_fun2, x_guess, method = :gradient_descent, grtol=1e-12, autodiff=true)


f2 = DifferentiableFunction(obj_fun2)
res4 = fminbox(f2, x_guess, lb, ub)

# It is not possible to use a defined objective function in JuMP
# The objective function is defined in one or more lines using `sum{f[i,j], i=1:n, j=1:m}
# It is apparently not yet mature enough.
# define the model

m = Model()
# define the variables of the foam model
@defVar(m, 10<=fmmob<=1000000)
@defVar(m, 10<=epdry<=100000)
@defVar(m, swc<=fmdry<=swc+0.3)
@defVar(m, 1.0e-6<=fmcap<=1.0e-4)
@defVar(m, -2.0<=epcap<=4.0)
three_par = zeros(3,4)
five_par = zeros(5,4)
for k=1:4
	# define the objective function
	@setNLObjective(m, Min, sum{(f_foam[1,i]-fmmob*(0.5+atan(epdry*(sw_exp[k,i]-fmdry))/pi))^2, i=1:8})
	# solve it using ipopt
	solve(m)
	# get the result
	x1 = [getValue(fmmob), getValue(epdry), getValue(fmdry)]
	three_par[:,k]=x1[:]
	# and visualize
	sw_val = linspace(0,1.0, 100)
	fg_opt = fg(sw_val, x1)
	muf_opt = mu_foam(sw_val, x1)
	figure(2)
	plot(fg_opt, muf_opt, fg_exp', muf_exp[k,:]', "o")
	# give the previous results as initial estimates to the new 5 par objective function
	setValue(fmmob, x1[1])
	setValue(epdry, x1[2])
	setValue(fmdry, x1[3])
	# define the objective function
	@setNLObjective(m, Min, sum{(f_foam[k,i]-fmmob*(0.5+atan(epdry*(sw_exp[k,i]-fmdry))/pi)
		*(fmcap/(u*muf_exp[k,i]/sigma_wg))^epcap)^2, i=1:8})
	# solve it
	solve(m)
	# get the results
	x2 = [getValue(fmmob), getValue(epdry), getValue(fmdry), getValue(fmcap), getValue(epcap)]
	five_par[:,k]=x2[:]
	# visualize the results
	n_data = length(fg_exp)
	muf_opt = zeros(n_data)
	for i=1:n_data
		fmu(x)= x-mu_foam2(sw_exp[k,i], x2[1:3], x2[4:5], x)
		muf_opt[i]=fzero(fmu, 0.1)
	end
	figure(3)
	plot(fg2(sw_exp[k,:][:], x2[1:3], x2[4:5], muf_opt), mu_foam2(sw_exp[k,:][:], x2[1:3], x2[4:5], muf_opt),
		fg_exp', muf_exp[k,:]', "o")
	axis([0,1,0,0.25])
end
