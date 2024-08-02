from setuptools import setup

setup(name = "snf_simulations",
      version = "0.1",
      description = "simulates antineutrino spectra for different compositions of spent nuclear fuel",
      packages = ["simulation"],
      install_requires = ['numpy', 'ROOT'],
      author='Zuzanna Leliwa',
      author_email='zleliwa1@sheffield.ac.uk', 
      zip_safe=False
)
   
