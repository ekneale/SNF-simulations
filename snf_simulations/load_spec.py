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

	new_edges = np.linspace(min_E, max_E, (max_E-min_E) + 1)
	new_centres = new_edges[:-1] + 0.5

	new_content = []
	
	for i in range(len(new_centres)):
		new_content.append(h.Interpolate(new_centres[i]))

	hnew = ROOT.TH1D(str(name) , isotope + "equal bin widths", len(new_centres), new_edges)

	for i in range(len(new_centres)):
		hnew.Fill(new_centres[i], new_content[i])
	
	#calculation of new errors for interpolated graph

	new_errors = []

	for i in range(0,len(new_centres)):
		# The TH1D::Interpolate function looks for the two closest bins surrounding the new point.
		# If it's before the first bin centre or above the final bin centre it just returns
		# the content of the first/last bin.
		# So here we do the same with the errors.
		if new_centres[i] < h.GetBinCenter(1):
			new_errors.append(h.GetBinError(1))
			continue
		if new_centres[i] >= h.GetBinCenter(h.GetNbinsX()):
			new_errors.append(h.GetBinError(h.GetNbinsX()))
			continue

		# find which of the old bins the new centre would be in
		idx = h.FindBin(new_centres[i])

		# find which of the surrounding bins is lower and which is upper
		# need to check which side of the centre of the found bin the new centre is
		if new_centres[i] < h.GetBinCenter(idx):
			comparison1 =[h.GetBinCenter(idx - 1), h.GetBinError(idx - 1)]
			comparison2 =[h.GetBinCenter(idx), h.GetBinError(idx)]
		else:
			comparison1 =[h.GetBinCenter(idx), h.GetBinError(idx)]
			comparison2 =[h.GetBinCenter(idx + 1), h.GetBinError(idx + 1)]

		d = comparison2[0] - comparison1[0]
		d1 = abs(new_centres[i] - comparison1[0])
		d2 = abs(comparison2[0] - new_centres[i])

		errors = np.sqrt(((d2/d)**2 * comparison1[1]**2) + ((d1/d)**2 * comparison2[1]**2))
		new_errors.append(errors)

	for i in range(1, len(new_centres)+1):
		hnew.SetBinError(i,new_errors[i-1])


	hnew.SetStats(0)

	hnew.GetXaxis().SetTitle("Energy [keV]")
	hnew.GetYaxis().SetTitle("#frac{dN}{dE} [keV^{-1}]" )
	
	hnew.GetXaxis().SetLabelSize(0.05)
	hnew.GetYaxis().SetLabelSize(0.05)

	hnew.GetXaxis().SetTitleSize(0.047)
	hnew.GetYaxis().SetTitleSize(0.047)


	return hnew
