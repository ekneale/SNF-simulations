from snf_simulations import plotting, flux

#plotting the spectrum for a single dry cask of fuel for either Hartlepool or Sizewell, depending which is set to True/False
spec_single = plotting.plot_single_cask([0.5,1,5,10,20], Sizewell= True) 
#flux.write_spec(spec_single, Sizewell = True) # unhash if you want energy/flux csv for single cask for either Hartlepool or Sizewell, depending which is set to True/False



#multiple casks
spec_multiple_sizewell = plotting.plot_multiple_casks_sizewell([0.5,5,10,20])
spec_multiple_hartlepool = plotting.plot_multiple_casks_hartlepool([3,7,15,19])
flux.write_spec(spec_multiple_sizewell, Sizewell=True) #unhash if you want energy/flux for multiple casks, change to _sizewell/_hartlepool if thats the reactor you want 

#calculation of flux at the detector located 40 m from the dry casks
flux.flux_calc(spec_multiple_sizewell, 40)


lots_of_graphs = plotting.multiple_fluxes(Sizewell= True)

plotting.plot_multiple(lots_of_graphs, Sizewell= True)

for i in range(len(lots_of_graphs)):
    flux.flux_calc(lots_of_graphs[i],40)

plotting.plot_sample(spec_multiple_sizewell, Sizewell= True)
