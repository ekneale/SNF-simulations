import ROOT
import define_proportions
import add_casks


def plot_single_cask(removal_times=[], Sizewell = False, HartlePool=False):
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)

    c.SetLogy()

    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    if Sizewell == True:
        initial = define_proportions.TotSpec(cask_name="main",Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 )

        initial.SetLineColor(ROOT.kBlue)

        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.append( define_proportions.TotSpec(cask_name="main",removal_time=removal_times[i],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8))

        legend.AddEntry(initial, "1 Day since removal from core")

        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], str(removal_times[i]) +" Years since removal from core")

        initial.Draw("hist")
        
        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

    if HartlePool == True:
        initial = define_proportions.TotSpec(cask_name="main",Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8)

        initial.SetLineColor(ROOT.kBlue)

        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.append( define_proportions.TotSpec(cask_name="main",removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        legend.AddEntry(initial, "1 Day since removal from core")

        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], str(removal_times[i]) +" Years since removal from core")
        
        initial.Draw("hist")
        
        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

    #input proportions are for 1 metric ton
    initial.SetTitle("")
    initial.GetXaxis().SetTitle("Energy [keV]")
    initial.GetXaxis().SetLabelSize(0.05)
    initial.GetXaxis().SetTitleSize(0.05)
    initial.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    initial.GetYaxis().SetLabelSize(0.05)
    initial.GetXaxis().SetTitleSize(0.05)

    c.Update()

    legend.Draw()

    input("exit")

    c.SaveAs("HartlepoolSpectra.pdf")



def plot_multiple_casks(Sizewell = False, HartlePool = False):
    if Sizewell == True:
        removal_times = [0.5,4,10,15,20,1,12,2.5,18,8,0.01, 10.5, 11, 3,7]

        casks = ROOT.TList()

        for i in range(len(removal_times)):
            casks.append(define_proportions.TotSpec(cask_name="Sizewell",removal_time=removal_times[i],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8))

        total = add_casks.add_casks(casks)

        return total

    if HartlePool == True:
        removal_times = [0.5, 10, 17, 9, 0.01, 11, 15, 20, 12, 17.5, 13,2,1, 6, 19]

        casks = ROOT.TList()

        for i in range(len(removal_times)):
            casks.Add(define_proportions.TotSpec(cask_name="Hartlepool",removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        total = add_casks.add_casks(casks)

        return total
    

def plot(spectrum):
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)
    c.SetLogy()
    spectrum.Draw("histS")

    spectrum.SetTitle("")
    spectrum.GetXaxis().SetTitle("Energy [keV]")
    spectrum.GetXaxis().SetLabelSize(0.05)
    spectrum.GetXaxis().SetTitleSize(0.05)
    spectrum.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    spectrum.GetYaxis().SetLabelSize(0.05)
    spectrum.GetXaxis().SetTitleSize(0.05)


    c.Update()
    input("exit")

    c.SaveAs("HartlepoolCasks.pdf")

