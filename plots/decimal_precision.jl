using SoftPosit, PyPlot, Lorenz63, StatsBase, Statistics

include("representable_numbers.jl")
include("lorenz_data.jl")

N = 100_000
L63data_s0 = lorenz_data(N,1/10)
L63data_s1 = lorenz_data(N,1.)
L63data_s2 = lorenz_data(N,10.)

##
float8 = representable_floats(8,3)
float16 = representable_floats(16,5)
bfloat16 = representable_floats(16,8)

f8_am, f8_wda = wcdp_float(float8)
f_am, f_wda = wcdp_float(float16)
bf_am, bf_wda = wcdp_float(bfloat16)

posit8 = Float64.(Posit8.(UInt8.(collect(1:127))))
posit16 = Float64.(Posit16.(UInt16.(collect(1:32767))))
posit16_2 = Float64.(Posit16_2.(UInt16.(collect(1:32767))))

p1_am, p1_wda = wcdp_posit(posit16)
p2_am, p2_wda = wcdp_posit(posit16_2)
p8_am, p8_wda = wcdp_posit(posit8)

i_am = vcat(0.5000001,1.49999999:1:(2^15-1))
i_wda = decprec(i_am,round.(i_am))
i_am[1] = 1

## histograms
bins = 10.0.^(-5:0.05:5)
H_s0 = fit(Histogram,abs.(L63data_s0[:]),bins)
w_s0 = H_s0.weights

H_s1 = fit(Histogram,abs.(L63data_s1[:]),bins)
w_s1 = H_s1.weights

H_s2 = fit(Histogram,abs.(L63data_s2[:]),bins)
w_s2 = H_s2.weights

## PLOTTING
ioff()
fig,(ax1,ax2) = subplots(2,1,figsize=(7,6),sharex=false)

ax1.plot(f_am,f_wda,"k",lw=2)
ax1.plot(bf_am,bf_wda,"0.7",lw=1.4)
ax1.plot(p1_am,p1_wda,"#50C070",lw=1.2)
ax1.plot(p2_am,p2_wda,"#900000",lw=0.8)
ax1.plot(p8_am,p8_wda,"C4",lw=1.8)
ax1.plot(f8_am,f8_wda,"C5",lw=1.3)
ax1.plot(i_am,i_wda,"C0",lw=2)
ax1.plot(i_am/2^10,i_wda,"C1",lw=2)

ax1.fill_between(f_am,-0.1,f_wda,edgecolor="k",facecolor="none",linestyle="--")
ax1.fill_between(f8_am,-0.1,f8_wda,edgecolor="C5",facecolor="none",linestyle="--")
ax1.fill_between(i_am,-0.1,i_wda,edgecolor="C0",facecolor="none",linestyle="--")
ax1.fill_between(i_am/2^10,-0.1,i_wda,edgecolor="C1",facecolor="none",linestyle="--")
ax1.fill_between(p1_am,-0.1,p1_wda,where=((p1_am .>= posit16[1]).*(p1_am .<= posit16[end])),edgecolor="C2",facecolor="none",linestyle="--")
ax1.fill_between(p8_am,-0.1,p8_wda,where=((p8_am .>= posit8[1]).*(p8_am .<= posit8[end])),edgecolor="C4",facecolor="none",linestyle="--")

# ax1.legend(loc=2,fontsize=9)

x0,x1 = 1e-16,1e16

ax1.set_xlim(x0,x1)
ax1.set_xscale("log",basex=10)
ax1.set_ylim(0,6)

ax2.set_xlabel("value")
ax1.set_ylabel("decimal places")
ax2.set_ylabel(L"$N$ [$10^5$]")

# ax1.text(7e-5,2.5,"Posit(16,0)",color="C1",rotation=63)
ax1.text(1e-12,0.5,"Posit(16,1)",color="#50C070")
ax1.text(5e-15,2,"Posit(16,2)",color="#900000")
ax1.text(6e-3,1,"Posit(8,0)",color="C4",rotation=64.5,va="bottom")
ax1.text(3e-2,.6,"Float8",color="C5",rotation=64.5,va="bottom")
ax1.text(5e-7,4.1,"Float16",color="k")
ax1.text(1e-15,3.2,"BFloat16",color="0.3")
ax1.text(6e4,5,"Int16",color="C0")
ax1.text(30,5.3,"Q6.10",color="C1",rotation=0)

# ax11 = ax1[:twiny]()
# ax11[:set_xscale]("log",basex=10)
# ax11[:set_xlim](ax1[:get_xlim]())

xtik = 10.0.^(-16:4:16)
# #ax1[:set_xticks](xtik)

ax2.set_xscale("log",basex=10)
ax2.set_xlim(x0,x1)
#ax2.set_xticks(xtik)
#ax1.set_xticks(xtik)

bins[1] = 1e-16
bins[end-1] = 1e16

ax2.plot(bins[1:end-1],w_s0/1e5,"0.8",label=L"s = 10$^{-1}$",drawstyle="steps-post")
ax2.plot(bins[1:end-1],w_s1/1e5,"0.45",label="s = 1",drawstyle="steps-post")
ax2.plot(bins[1:end-1],w_s2/1e5,"0.0",label=L"s = 10$^1$",drawstyle="steps-post")
ax2.legend(loc=2)

ax1.set_title("Decimal precision",loc="left")
ax2.set_title("Numbers in Lorenz63",loc="left")

ax1.set_title("a",loc="right",fontweight="bold")
ax2.set_title("b",loc="right",fontweight="bold")

tight_layout()
savefig("plots/decimal_precision.png",dpi=300)
savefig("plots/decimal_precision.pdf")
close(fig)
