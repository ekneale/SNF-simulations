# About

SNF-simulations generates antineutrino spectra that can be used to simulate the antineutrino emission from spent nuclear fuel (SNF).

As the world moves away from reliance on fossil fuels, nuclear power generation will increase, leading to unprecedented quantities of highly radioactive nuclear waste. After initial cooling, this SNF is stored in above-ground dry casks at interim storage facilities, where it is monitored for decades for safety and nonproliferation. Effective monitoring technologies are of great importance, as interim facilities store hundreds of tonnes of highly radioactive spent fuel, and increasing reliance on civil nuclear power has led to more interest in developing improved methods of monitoring SNF. One such method is the detection of electron antineutrinos (antineutrinos) emitted by isotopes undergoing beta decay. Predicting the emission from a SNF storage facility is vital in understanding the expected signal and capabilities of a neutrino detector for SNF monitoring.

This software package has been created to predict the antineutrino emission spectrum from any SNF storage facility of any number of dry storage casks with any initial isotopic composition and cooling time after removal from the core. It takes the isotopic composition for the most relevant isotopes from FISPIN[1] calculations of nuclear fuel after burn up in a reactor core, and outputs the total antineutrino spectrum, for input into detector simulations. The user can define the number of casks and cooling time for each cask.

The decay rate for each isotope after the cooling period is calculated from the initial activity and the decay constant of the isotope. Where the decay chain contains further beta-emitting isotopes, the activity of the beta-decaying daughter isotope also takes into account the amount of the parent isotope that has decayed during the cooling period.

Calculations of the antineutrino emission from each isotope are based on antineutrino spectrum data for each decay, taken from IAEA Nuclear Data Services Live Chart of Nuclides[2], which are based on Evaluated Nuclear Structure Data File (ENSDF)[3] and generated using Betashape[4].

The total calculated spectra are finally sampled to generate antineutrino avents that can be used in your detector simulation.

The software can also be used to provide an approximate prediction for the sensitivity of any given neutrino detector to SNF neutrino emission, according to user-defined values for the detector efficiency, size and distance from the SNF source.

[1] R. F. Burstall and UKAEA Risley Nuclear Power Development Establishment (Oct. 1979). FISPIN - a computer code for nuclide inventory calculations. Northern Division Report, United Kingdom Atomic Energy Authority Northern Division.

[2] IAEA Nuclear Data Section, 'Live Chart of Nuclides', https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html. (Accessed August 2024).

[3] M. R. Bhat (1992). Evaluated Nuclear Structure Data File (ENSDF). In: Qaim, S.M. (eds) Nuclear Data for Science and Technology. Research Reports in Physics. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-58113-7_227

[4] X. Mougeot (2017). BetaShape: A new code for improved analytical calculations of beta spectra. EPJ Web Conf., 146 12015
