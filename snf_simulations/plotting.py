import ROOT
from snf_simulations import define_proportions, add_casks, sample

 
# the things that have been changed so far: 
# for the single casks, the energy and flux will save to a csv file. You can pick which cooling time by changing the number in the spectra.At(). 
# You can pick the reactor in command line, the name of the reactor is now a variable so you don't have to manually change the pdf name .
# The multiple flux function has been updated so that sizewell is defined as well as Hartlepool. The plot spectrum function also has the reactor variable.
# Two seperate functions are defined for plotting multiple casks for either sizewell or hartlepool as they have different removal times so energy/flux csv can be found per time
#At the moment you have to manually change the csv file name for what reactor, cooling time etc but eventually this would be done automatically.
 

def plot_single_cask(removal_times=[], Sizewell = False, HartlePool=False):
    c = ROOT.TCanvas("c","Total Spectrum", 1200,600)
    c.SetLogy()
    legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    #plotting a spectrum of antineutrino energies for the Sizewell PWR with varying cooling times
    if Sizewell == True:
        reactor = "Sizewell"
        initial = define_proportions.TotSpec(cask_name="main",total_m = 100000, Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 )
        initial.SetLineColor(ROOT.kBlue)
        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.Add( define_proportions.TotSpec(cask_name="main",total_m = 100000,removal_time=removal_times[i],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8))

        legend.AddEntry(initial, "1 Day since removal from core")
        colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kViolet, ROOT.kOrange]

        for i in range(len(spectra)):
            spectra[i].SetLineColor(colors[i])
            legend.AddEntry(spectra[i], str(removal_times[i]) +" Years since removal from core")

        initial.Draw("hist")

        for i in range(len(spectra)):
            spectra[i].Draw("hist same")

    #plotting a spectrum of antineutrino energies for the Hartlepool AGR with varying cooling times
    if Sizewell == False:  
        reactor = "Hartlepool"
        initial = define_proportions.TotSpec(cask_name="main",total_m = 100000,Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8)
        initial.SetLineColor(ROOT.kBlue)
        spectra = ROOT.TList()

        for i in range(len(removal_times)):
            spectra.Add( define_proportions.TotSpec(cask_name="main",total_m = 100000,removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))


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
    c.SaveAs(f"{reactor}_Spectra_0.5.C") #change no for cooling time

    return (spectra.At(0))  #needs .At() as its a Tlist not a python array



 
def plot_multiple_casks_sizewell(removal_times):

    #plotting 10 dry casks x 10tonnes  x 4 cooling times for Sizewell PWR
    removal_times = [0.5,5,10,20]
    casks = ROOT.TList()
    
    for i in range(len(removal_times)):
        casks.Add(define_proportions.TotSpec(cask_name="Sizewell",total_m = 100000, removal_time=removal_times[i],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8))

    total_sizewell = add_casks.add_casks(casks)

    #c.SaveAs("Sizewell_MultipleCasks.C")

    return total_sizewell 
 
def plot_multiple_casks_hartlepool(removal_times):
      
        removal_times = [3,7,15,19]
        casks = ROOT.TList()

        for i in range(len(removal_times)):
            casks.Add(define_proportions.TotSpec(cask_name="Hartlepool",total_m = 100000,removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        total_hartlepool = add_casks.add_casks(casks)

        #c.SaveAs("Hartlepool_MultipleCasks.C")

        return total_hartlepool 
 



 

#just using the plot_multiple_casks results in a graphics error due to the amount spectra being added so plot has to be saved separately and then plotted

def plot(spectrum,Hartlepool = False, Sizewell = False):
    if Hartlepool == True:
         reactor = "Hartlepool"
    if Sizewell == True:
         reactor = "Sizewell"

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
    c.SaveAs(f"{reactor}_casks.C")

 

 

#plotting previous spectrum for 10 dry casks x 10 tonned x 4 cooling times for varying cooling times since the initial measurement
def multiple_fluxes(Sizewell = False, Hartlepool = False):
    extra_times = [1,5,10,20]
    
        #for Hartlepool
    if Hartlepool == True:
        removal_times = [3,7,15,19]
        casks0_h = ROOT.TList()
        casks1_h = ROOT.TList()
        casks5_h = ROOT.TList()
        casks10_h = ROOT.TList()
        casks20_h= ROOT.TList()
        for i in range(len(removal_times)):
                casks0_h.Add(define_proportions.TotSpec(cask_name="Hartlepool0",total_m = 100000,removal_time=removal_times[i],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        sum0 = add_casks.add_casks(casks0_h)  

        for i in range(len(removal_times)):
                casks1_h.Add(define_proportions.TotSpec(cask_name="Hartlepool1",total_m = 100000,removal_time=removal_times[i]+extra_times[0],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        sum1 = add_casks.add_casks(casks1_h)  

        for i in range(len(removal_times)):
                casks5_h.Add(define_proportions.TotSpec(cask_name="Hartlepool2",total_m = 100000,removal_time=removal_times[i]+extra_times[1],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        sum5 = add_casks.add_casks(casks5_h)  

        for i in range(len(removal_times)):
                casks10_h.Add(define_proportions.TotSpec(cask_name="Hartlepool3",total_m = 100000,removal_time=removal_times[i]+extra_times[2],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))

        sum10 = add_casks.add_casks(casks10_h)  

        for i in range(len(removal_times)):
                casks20_h.Add(define_proportions.TotSpec(cask_name="Hartlepool4",total_m = 100000,removal_time=removal_times[i]+extra_times[3],Kr88_prop=6.4374e-11, Rb88_prop=7.5089e-11, Sr90_prop=4.2912e-4, Y90_prop=1.0953e-7, Zr93_prop=5.5068e-4, Tc99_prop=6.12e-4, Ru106_prop=6.1306e-5, Rh106_prop=5.7428e-11, I129_prop=1.2097e-4, Cs135_prop=3.8379e-4, Cs137_prop=8.6117e-4,Ce144_prop=1.7271e-4, Pr144_prop=7.2749e-9, Np239_prop=2.0904e-5, Pu241_prop=5.3075e-4, Am242_prop=1.3708e-8))
            
        sum20 = add_casks.add_casks(casks20_h)  


        sums_h = ROOT.TList()#
        sums_h.Add(sum0)
        sums_h.Add(sum1)
        sums_h.Add(sum5)
        sums_h.Add(sum10)
        sums_h.Add(sum20)
        return sums_h

#for sizewell
    if Sizewell == True:
        removal_times = [0.5,5,10,20]
        casks0_s = ROOT.TList()
        casks1_s = ROOT.TList()
        casks5_s = ROOT.TList()
        casks10_s = ROOT.TList()
        casks20_s = ROOT.TList()

        for i in range(len(removal_times)):
            casks0_s.Add(define_proportions.TotSpec(cask_name="Sizewell0",total_m = 100000,removal_time=removal_times[i], Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 ))

        sum0_s = add_casks.add_casks(casks0_s)  

        for i in range(len(removal_times)):
            casks1_s.Add(define_proportions.TotSpec(cask_name="Sizewell1",total_m = 100000,removal_time=removal_times[i]+extra_times[0], Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 ))


        sum1_s = add_casks.add_casks(casks1_s)  

        for i in range(len(removal_times)):
            casks5_s.Add(define_proportions.TotSpec(cask_name="Sizewell2",total_m = 100000,removal_time=removal_times[i]+extra_times[1], Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 ))


        sum5_s = add_casks.add_casks(casks5_s)  

        for i in range(len(removal_times)):
            casks10_s.Add(define_proportions.TotSpec(cask_name="Sizewell3",total_m = 100000,removal_time=removal_times[i]+extra_times[2], Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 ))


        sum10_s = add_casks.add_casks(casks10_s)  

        for i in range(len(removal_times)):
            casks20_s.Add(define_proportions.TotSpec(cask_name="Sizewell4",total_m = 100000,removal_time=removal_times[i]+extra_times[3],Kr88_prop=1.427e-10, Rb88_prop=1.6645e-11, Sr90_prop=5.356e-4,Y90_prop=1.3922e-7, Zr93_prop=1.7681e-6,Tc99_prop=7.9742e-4, Ru106_prop=1.7496e-4, Rh106_prop=1.6389e-10, I129_prop=1.7535e-4, Cs135_prop=3.1282e-4, Cs137_prop=1.212e-3, Ce144_prop=4.0111e-4, Pr144_prop=1.6896e-8, Np239_prop=7.5852e-5, Pu241_prop=1.316e-3,Am242_prop=3.554e-8 ))
        
        sum20_s = add_casks.add_casks(casks20_s)

        sums_s = ROOT.TList()#
        sums_s.Add(sum0_s)
        sums_s.Add(sum1_s)
        sums_s.Add(sum5_s)
        sums_s.Add(sum10_s)
        sums_s.Add(sum20_s)
        return sums_s
    

 
#plotting separately to avoid graphics errors

def plot_multiple(multiple):
   # if Hartlepool == True:
    #     reactor = "Hartlepool"
    #if Sizewell == True:
    #     reactor == "Sizewell"

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

    c.SaveAs("Sizewell_MultipleCasks.pdf")

 

def plot_sample(totSpec):
   # if Hartlepool == True:
    #     reactor = "Hartlepool"
   # if Sizewell == True:
    #     reactor == "Sizewell"
         
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
    c.SaveAs("Sizewell_Sampled.pdf")