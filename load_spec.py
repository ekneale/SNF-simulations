import ROOT
import numpy as np


def load_spec(Energy,dN,isotope):
	c = ROOT.TCanvas()

	h = ROOT.TH1D(isotope, isotope, len(Energy)-1, np.array(Energy))
		
	for i in range(len(Energy)):
		hist.Fill(Energy[i],dN[i])
	
	hist.SetStats(0)

	hist.GetXaxis().SetTitle("Energy [keV]")
	hist.GetYaxis().SetTitle("dN/dE")

	c.Update()

	hist.Draw("hist")
	input("exit")

	return hist


def load_equal(name, isotope,E,dN, max_E, min_E=0,):
	c = ROOT.TCanvas()

	h = ROOT.TH1D("h","", len(E)-1, np.array(E))

	for i in range(len(E)):
		h.Fill(E[i],dN[i])
	
	new_edges = np.linspace(min_E, max_E, max_E +1)

	new_centres = np.linspace(0.5, max_E-0.5, max_E)

	new_content = []
	
	for i in range(len(new_centres)):
		new_content.append(h.Interpolate(new_centres[i]))

	hnew = ROOT.TH1D(str(name) , isotope + "equal bin widths", len(new_centres), new_edges)

	for i in range(len(new_centres)):
		hnew.Fill(new_centres[i], new_content[i])
	
	hnew.SetStats(0)

	hnew.GetXaxis().SetTitle("Energy [keV]")
	hnew.GetYaxis().SetTitle("dN/dE")
	
	c.Update()
	

	return hnew
		





