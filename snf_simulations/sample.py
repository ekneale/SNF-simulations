import numpy as np
from ROOT import TH1D, TCanvas
import csv

#sampling the total spectrum to simulate what a detector could observe
#samples can be saved as csv for further analysis

def sample(total_spec, N=100, save = True):
    sampled = []

    for i in range(N):
        sampled.append(total_spec.GetRandom())

    if save == True:
        with open('sampled_spectrum.csv', mode = 'w', newline = '') as file:
            writer = csv.writer(file)
            for item in sampled:
                writer.writerow([item])

        print('saved as csv')

    return sampled
