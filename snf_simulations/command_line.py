from snf_simulations import plotting, flux


#plotting the spectrum for a single dry cask of fuel for either Hartlepool or Sizewell, depending which is set to True
plotting.plot_single_cask([0.5,1,5,10,20], HartlePool = True)


#plotting the spectrum for 10 x 10 tonne x 4 cooling times dry casks of SNF from eithe Sizewell or Hartlepool depending on which is set to True
spec = plotting.plot_multiple_casks(HartlePool=True)

plotting.plot(spec)


#calculation of flux at the detector located 40 m from the dry casks
flux.flux_calc(spec, 40)


lots_of_graphs = plotting.multiple_fluxes()

plotting.plot_multiple(lots_of_graphs)

for i in range(len(lots_of_graphs)):
    flux.flux_calc(lots_of_graphs[i],40)

plotting.plot_sample(spec)
