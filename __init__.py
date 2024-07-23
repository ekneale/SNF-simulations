import numpy as np


def load_data():
    Sr90 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Sr90_an.txt",skip_header=1)[:,[7,10]]

    Y90 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Y90_an.txt",skip_header=1)
    end =int(np.where(Y90[:,7] == 2278.5)[0])
    Y90 = Y90[:end,[7,10]]

    Pu241 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Pu241_an.txt",skip_header=1)[:,[7,10]]
    Cs137 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Cs137_an.txt",skip_header=1)[:,[7,10]]
    Am242 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Am242_an.txt",skip_header=1)[:,[7,10]]
    Cm249 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Cm249_an.txt",skip_header=1)[:,[7,10]]
    Cs135 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Cs135_an.txt",skip_header=1)[:,[7,10]]
    I129 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/I129_an.txt",skip_header=1)[:,[7,10]]
    Np239 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Np239_an.txt",skip_header=1)[:,[7,10]]

    Tc99 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Tc99_an.txt",skip_header=1)
    end0 = int(np.where(Tc99[:,7] == 297.5)[0])
    Tc99 = Tc99[:end0,[7,10]]

    Zr93 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Zr93_an.txt",skip_header=1)[:,[7,10]]

    Ce144 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Ce144_an.txt",skip_header=1)[:,[7,10]]
    Kr88 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Kr88_an.txt",skip_header=1)[:,[7,10]]
    Pr144 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Pr144_an.txt",skip_header=1)[:,[7,10]]
    Rb88 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Rb88_an.txt",skip_header=1)[:,[7,10]]
    Rh106 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Rh106_an.txt",skip_header=1)[:,[7,10]]
    Ru106 = np.genfromtxt("/home/zuzanna/SURE/snf_sims/plot_spectra/data/Ru106_an.txt",skip_header=1)[:,[7,10]]




    return(Sr90, Y90,Pu241, Cs137, Am242, Cm249, Cs135, I129, Np239, Tc99, Zr93, Ce144, Kr88, Pr144, Rb88, Rh106, Ru106)