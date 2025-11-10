import ROOT
import numpy as np


def load_spec(Energy,dN,errors,isotope):
	c = ROOT.TCanvas()

	#the energies are the bin edges and the bin centres are between these, however these graphs have unequal bin widths
	h = ROOT.TH1D(isotope, isotope, len(Energy)-1, np.array(Energy))
		
	for i in range(len(Energy)):
		h.Fill(Energy[i],dN[i])

	# Fill the bins (note that ROOT histograms are 1-indexed)
	for i in range(1, len(Energy)+1):
		h.SetBinError(i, errors[i-1])

	h.SetStats(0)

	h.GetXaxis().SetTitle("Energy [keV]")
	h.GetYaxis().SetTitle("dN/dE")

	c.Update()

	h.Draw("hist E")
	input("exit")

	return h


def load_equal(name, isotope,E,dN, error, max_E, min_E=0,):
	h = ROOT.TH1D("h","", len(E)-1, np.array(E))
	
	#loading in original spectra

	for i in range(len(E)):
		h.Fill(E[i],dN[i])
	for i in range(1, len(E)+1):
		h.SetBinError(i, error[i-1])

	#creating new bin edges and new bin centres to ensure equal bin widths so that the histograms can later be added together

	new_edges = np.linspace(min_E, max_E, max_E +1)

	new_centres = np.linspace(0.5, max_E-0.5, max_E)

	new_content = []
	
	for i in range(len(new_centres)):
		new_content.append(h.Interpolate(new_centres[i]))

	hnew = ROOT.TH1D(str(name) , isotope + "equal bin widths", len(new_centres), new_edges)

	for i in range(len(new_centres)):
		hnew.Fill(new_centres[i], new_content[i])
	
	#calculation of new errors for interpolated graph

	new_errors = []

	#seperate calculation of error on first interpolated point as iteration has to start from 1

	point1 = [h.GetBinCenter(0), h.GetBinError(0)]
	point2 = [h.GetBinCenter(1), h.GetBinError(1)]
	d0 = point2[0] - point1[1]
	d01 = abs(new_centres[0] - point1[0])
	d02 = abs(new_centres[0] - point2[0])

	first_err = np.sqrt(((d01/d0)*(point1[1])**2) + ((d02/d0)*(point2[1]**2)))
	new_errors.append(first_err)

	#calcualtion of rest of interpolated errors

	for i in range(1,len(new_centres)):
		comparison1 =[h.GetBinCenter(i-1), h.GetBinError(i-1)]
		comparison2 =[h.GetBinCenter(i+1), h.GetBinError(i+1)]

		d = comparison2[0] - comparison1[0]
		d1 = abs(new_centres[i] - comparison1[0])
		d2 = abs(comparison2[0] - new_centres[i])

		errors = np.sqrt(((d1/d)*(comparison1[1])**2) + ((d2/d)*(comparison2[1]**2)))
		new_errors.append(errors)

	for i in range(len(new_centres)):
		hnew.SetBinError(i,new_errors[i])

	hnew.SetStats(0)

	hnew.GetXaxis().SetTitle("Energy [keV]")
	hnew.GetYaxis().SetTitle("#frac{dN}{dE} [keV^{-1}]" )
	
	hnew.GetXaxis().SetLabelSize(0.05)
	hnew.GetYaxis().SetLabelSize(0.05)

	hnew.GetXaxis().SetTitleSize(0.047)
	hnew.GetYaxis().SetTitleSize(0.047)


	return hnew
