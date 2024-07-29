import ROOT
import define_proportions
import add_casks

def plot_single_cask(removal_times=[], Sizewell = False, HartlePool=False):
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)

    c.SetLogy()

    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    if Sizewell == True:
        initial = define_proportions.TotSpec(cask_name="main",Kr88_prop=4.972e-6,Rb88_prop=5.3151e-8,Sr90_prop=4.2217e-4, Y90_prop= 1.5087e-7, Zr93_prop=7.0332e-4, Tc99_prop=.0080e-4, Ru106_prop=6.0510e-4, Rh106_prop=1.69e-13, I129_prop=1.7722e-4, Ce144_prop=5.555e-8, Pr144_prop=2.3399e-12, Np239_prop=7.2309e-11, Pu241_prop=8.1328e-4, Am242_prop=5.2489e-12)

        initial.SetLineColor(ROOT.kBlue)

        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.append( define_proportions.TotSpec(cask_name="main",removal_time=removal_times[i],Kr88_prop=4.972e-6,Rb88_prop=5.3151e-8,Sr90_prop=4.2217e-4, Y90_prop= 1.5087e-7, Zr93_prop=7.0332e-4, Tc99_prop=.0080e-4, Ru106_prop=6.0510e-4, Rh106_prop=1.69e-13, I129_prop=1.7722e-4, Ce144_prop=5.555e-8, Pr144_prop=2.3399e-12, Np239_prop=7.2309e-11, Pu241_prop=8.1328e-4, Am242_prop=5.2489e-12))

        legend.AddEntry(initial, "Initial Spectrum after removal from core")

        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], "Spectrum after " + str(removal_times[i]) +" years")

        initial.Draw("hist")
        
        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

    if HartlePool == True:
        initial = define_proportions.TotSpec(cask_name="main",Kr88_prop=2.243e-8,Rb88_prop=2.3878e-9, Sr90_prop=4.2915e-4,Y90_prop=1.1011e-1, Zr93_prop=5.5053e-4, Tc99_prop=6.1165e-4,Ru106_prop=6.1422e-5, Rh106_prop=5.8187e-11, I129_prop=1.2094e-4,Cs135_prop=3.8357e-4,Cs137_prop=8.6122e-4,Pr144_prop=7.3156e-9,Ce144_prop=1.7313e-4, Np239_prop=2.7864e-5,Pu241_prop=5.3082e-4,Am242_prop=3.8686e-8,Cm249_prop=3.1726e-18)

        initial.SetLineColor(ROOT.kBlue)

        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.append( define_proportions.TotSpec(cask_name="main",removal_time=removal_times[i],Kr88_prop=2.243e-8,Rb88_prop=2.3878e-9, Sr90_prop=4.2915e-4,Y90_prop=1.1011e-1, Zr93_prop=5.5053e-4, Tc99_prop=6.1165e-4,Ru106_prop=6.1422e-5, Rh106_prop=5.8187e-11, I129_prop=1.2094e-4,Cs135_prop=3.8357e-4,Cs137_prop=8.6122e-4,Pr144_prop=7.3156e-9,Ce144_prop=1.7313e-4, Np239_prop=2.7864e-5,Pu241_prop=5.3082e-4,Am242_prop=3.8686e-8,Cm249_prop=3.1726e-18))

        legend.AddEntry(initial, "Initial Spectrum after removal from core")

        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], "Spectrum after " + str(removal_times[i]) +" years")
        
        initial.Draw("hist")
        
        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

    #input proportions are for 1 metric ton

    initial.GetXaxis().SetTitle("Energy [keV]")
    initial.GetXaxis().SetLabelSize(0.05)
    initial.GetXaxis().SetTitleSize(0.05)
    initial.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    initial.GetYaxis().SetLabelSize(0.05)
    initial.GetXaxis().SetTitleSize(0.05)

    c.Update()

    legend.Draw()

    input("exit")




def plot_multiple_casks():
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)

    c.SetLogy()


    #Sizewell
    cask1 = define_proportions.TotSpec(cask_name="cask1",removal_time=10,Kr88_prop=4.972e-6,Rb88_prop=5.3151e-8,Sr90_prop=4.2217e-4, Y90_prop= 1.5087e-7, Zr93_prop=7.0332e-4, Tc99_prop=.0080e-4, Ru106_prop=6.0510e-4, Rh106_prop=1.69e-13, I129_prop=1.7722e-4,Cs135_prop=3.1232e-4, Cs137_prop=1.2121e-3, Ce144_prop=5.555e-8, Pr144_prop=2.3399e-12, Np239_prop=7.2309e-11, Pu241_prop=8.1328e-4, Am242_prop=5.2489e-12)
    cask1.SetLineColor(ROOT.kGreen)

    #HartlePool
    cask2 = define_proportions.TotSpec(cask_name="cask2",removal_time=5,Kr88_prop=2.243e-8,Rb88_prop=2.3878e-9, Sr90_prop=4.2915e-4,Y90_prop=1.1011e-1, Zr93_prop=5.5053e-4, Tc99_prop=6.1165e-4,Ru106_prop=6.1422e-5, Rh106_prop=5.8187e-11, I129_prop=1.2094e-4,Cs135_prop=3.8357e-4,Cs137_prop=8.6122e-4,Pr144_prop=7.3156e-9,Ce144_prop=1.7313e-4, Np239_prop=2.7864e-5,Pu241_prop=5.3082e-4,Am242_prop=3.8686e-8,Cm249_prop=3.1726e-18 )
    cask2.SetLineColor(ROOT.kRed)

    casks = ROOT.TList()
    casks.Add(cask1)
    casks.Add(cask2)

    sum = add_casks.add_casks(casks)

    sum.SetLineColor(ROOT.kBlue)

    sum.Draw("hist")

    cask1.Draw("hist same")
    cask2.Draw("hist same")



    spectra = ROOT.TList()
        
    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    legend.AddEntry(sum, "Sum of Dry Cask Emissions")
    legend.AddEntry(cask1, "Sizewell Reactor 10 years after removal")
    legend.AddEntry(cask2, "HartlePool Reactor 5 years after removal")

    legend.Draw()

    sum.GetXaxis().SetTitle("Energy [keV]")
    sum.GetXaxis().SetLabelSize(0.05)
    sum.GetXaxis().SetTitleSize(0.05)
    sum.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    sum.GetYaxis().SetLabelSize(0.05)
    sum.GetXaxis().SetTitleSize(0.05)


    input("exit")

