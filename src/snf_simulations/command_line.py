from . import flux, plotting

# from ROOT import TFile

# creating root file in order to export to simulation
# hfile = TFile("Sizewell_Spectra.root", "RECREATE")

# plotting the spectrum for a single dry cask of fuel for either Hartlepool or Sizewell, depending which is set to True/False
spec_single = plotting.plot_single_cask(
    reactor="sizewell",
    removal_times=[0.5, 1, 5, 10, 20],
)
energy_single, flux_single = flux.write_spec_single(spec_single, reactor="sizewell")


# multiple casks
spec_multiple_sizewell = plotting.plot_multiple_casks_sizewell([0.5, 5, 10, 20])
spec_multiple_hartlepool = plotting.plot_multiple_casks_hartlepool([3, 7, 15, 19])
energy_multiple, flux_multiple = flux.write_spec_multiple(
    spec_multiple_sizewell,
    reactor="sizewell",
)  # change for sizewell or hartlepool


# calculation of flux at the detector located 40 m from the dry casks
flux.flux_calc(spec_multiple_sizewell, 40)
flux.flux_calc(spec_single, 40)


# plot both flux spectra on one graph
flux.multiple_single_plot(
    energy_single,
    flux_single,
    energy_multiple,
    flux_multiple,
)

lots_of_graphs = plotting.multiple_fluxes(reactor="sizewell")

plotting.plot_multiple(lots_of_graphs, reactor="sizewell")

for i in range(len(lots_of_graphs)):
    flux.flux_calc(lots_of_graphs[i], 40)

plotting.plot_sample(spec_multiple_sizewell, reactor="sizewell")

# hfile.Write()
