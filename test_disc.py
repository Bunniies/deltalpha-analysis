import numpy 

pp  = '/Users/alessandroconigli/Lattice/data/HVP/disc/disc/2pt/88.npz'
data = numpy.load(pp, allow_pickle=True)
numpy.__version__
tt = data['VV']