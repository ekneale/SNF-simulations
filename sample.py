import numpy as np
from ROOT import TH1D, TCanvas
import csv



def sample(total_spec, N=100):
    sampled = []
    sampled_e = []

    for i in range(N):
        sampled.append(total_spec.GetRandom())

    with open('sampled_spectrum.csv', mode = 'w', newline = '') as file:
        writer = csv.writer(file)

        for item in sampled:
            writer.writerow([item])


    print('saved as csv')
