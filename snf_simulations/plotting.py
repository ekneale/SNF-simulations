import ROOT
from snf_simulations import define_proportions, add_casks, sample


def plot_single_cask(removal_times=[], Sizewell = False, HartlePool=False):
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)

    c.SetLogy()

    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    #plotting a spectrum of antineutrino energies for the Sizewell PWR with varying cooling times

    if Sizewell == True:
        initial = define_proportions.TotSpec(cask_name="main",total_m = 100000, Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 )

        initial.SetLineColor(ROOT.kBlue)

        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.append( define_proportions.TotSpec(cask_name="main",total_m = 100000,removal_time=removal_times[i],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8))

        legend.AddEntry(initial, "1 Day since removal from core")

        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], str(removal_times[i]) +" Years since removal from core")

        initial.Draw("hist")
        
        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

    #plotting a spectrum of antineutrino energies for the Hartlepool AGR with varying cooling times

    if HartlePool == True:
        initial = define_proportions.TotSpec(cask_name="main",total_m = 100000,Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8)

        initial.SetLineColor(ROOT.kBlue)

        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.append( define_proportions.TotSpec(cask_name="main",total_m = 100000,removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        legend.AddEntry(initial, "1 Day since removal from core")

        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], str(removal_times[i]) +" Years since removal from core")
        
        initial.Draw("hist")
        
        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

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

    #plotting 10 dry casks x 10tonnes  x 4 cooling times for Sizewell PWR

    if Sizewell == True:
        removal_times = [0.5,5,10,20]

        casks = ROOT.TList()

        for i in range(len(removal_times)):
            casks.append(define_proportions.TotSpec(cask_name="Sizewell",total_m = 100000, removal_time=removal_times[i],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8))

        total = add_casks.add_casks(casks)

        return total

    #plotting 10 dry casks x 10tonnes  x 4 cooling times for Hartlepool AGR

    if HartlePool == True:
        removal_times = [3,7,15,19]

        casks = ROOT.TList()

        for i in range(len(removal_times)):
            casks.Add(define_proportions.TotSpec(cask_name="Hartlepool",total_m = 100000,removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        total = add_casks.add_casks(casks)

        return total
    

#just using the plot_multiple_casks results in a graphics error due to the amount spectra being added so plot has to be saved separately and then plotted

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

    c.SaveAs("Hartlepoolcasks.pdf")


#plotting previous spectrum for 10 dry casks x 10 tonned x 4 cooling times for varying cooling times since the inital measurement
def multiple_fluxes():
    removal_times = [3,7,15,19]
     
    extra_times = [1,5,10,20]
    casks0 = ROOT.TList()
    casks1 = ROOT.TList()
    casks5 = ROOT.TList()
    casks10 = ROOT.TList()
    casks20 = ROOT.TList()

    #for Hartlepool

    for i in range(len(removal_times)):
        casks0.Add(define_proportions.TotSpec(cask_name="Hartlepool0",total_m = 100000,removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

    sum0 = add_casks.add_casks(casks0)  

    for i in range(len(removal_times)):
        casks1.Add(define_proportions.TotSpec(cask_name="Hartlepool1",total_m = 100000,removal_time=removal_times[i]+extra_times[0],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

    sum1 = add_casks.add_casks(casks1)  

    for i in range(len(removal_times)):
        casks5.Add(define_proportions.TotSpec(cask_name="Hartlepool2",total_m = 100000,removal_time=removal_times[i]+extra_times[1],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

    sum5 = add_casks.add_casks(casks5)  

    for i in range(len(removal_times)):
        casks10.Add(define_proportions.TotSpec(cask_name="Hartlepool3",total_m = 100000,removal_time=removal_times[i]+extra_times[2],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

    sum10 = add_casks.add_casks(casks10)  

    for i in range(len(removal_times)):
        casks20.Add(define_proportions.TotSpec(cask_name="Hartlepool4",total_m = 100000,removal_time=removal_times[i]+extra_times[3],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))
    
    sum20 = add_casks.add_casks(casks20)  

    sums = ROOT.TList()#
    sums.Add(sum0)
    sums.Add(sum1)
    sums.Add(sum5)
    sums.Add(sum10)
    sums.Add(sum20)

    return sums


#plotting separately to avoid graphics errors 
def plot_multiple(multiple):
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)
    c.SetLogy()

    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    extra_times = [0,1,5,10,20]

    multiple[0].SetLineColor(ROOT.kBlue)

    multiple[0].Draw("hist")

    legend.AddEntry(multiple[0], "Total spectrum")

    colors = [ROOT.kYellow, ROOT.kGreen, ROOT.kViolet, ROOT.kBlack, ROOT.kRed]

    for i in range (1,len(multiple)):
        multiple[i].SetLineColor(colors[i])
        multiple[i].Draw("hist same")
        legend.AddEntry(multiple[i], "Spectrum after " + str(extra_times[i]) +" years" )

    multiple[0].SetTitle("")
    multiple[0].GetXaxis().SetTitle("Energy [keV]")
    multiple[0].GetXaxis().SetLabelSize(0.05)
    multiple[0].GetXaxis().SetTitleSize(0.05)
    multiple[0].GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    multiple[0].GetYaxis().SetLabelSize(0.05)
    multiple[0].GetXaxis().SetTitleSize(0.05)

    legend.Draw()

    c.Update()
    input("exit")

    c.SaveAs("HartlepoolMultipleCasks.pdf")

def plot_sample(totSpec):
    h = ROOT.TH1D("sample", "", 6000,0,6000)

    sampled = sample.sample(totSpec, N=1000000)

    for i in range(1000000):
        h.Fill(sampled[i])
        
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)

    c.SetLogy()

    h.Draw("hist")

    h.SetStats(0)
    h.SetTitle("")
    h.GetXaxis().SetTitle("Energy [keV]")
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitle("Relative Flux [keV ^{-1} s^{-1}]")
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)

    c.Update()

    input("exit")

    c.SaveAs("HartlepoolSampled.pdf")

        
        
