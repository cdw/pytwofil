# $filaments_and_heads
#     $Used to instantiate thick and thin filaments
###############################################################################
# Revisions
#    v 0.01/cdw     20090201  (Re)Created
# 
###############################################################################

import random as rn
from math import *

class XB:
	"""A crossbridge that will be instantiated by the thick filament"""
	def __init__(self, head_id, thick_fil, thin_fil):
		# The first thing to do is define our default parameters.
		self.Cs = pi/3 # rest angle of converter domain
		self.Ck = 200  # torsional spring const of converter domain
		self.Gs = 10.5 # rest length of globular domain
		self.Gk = 5    # spring constant of globular domain
		self.Fm = 0    # mean of forces exerted on myosin heads
		self.Fv = 10.0 # variance of forces exerted on myosin heads
		self.Bd = 0.55 # dist at which binding becomes likely
		# Thick and thin fil connectivity
		self.thick = thick_fil
		self.thin = thin_fil
		# Current state and identity of XB
		self.id = head_id
		self.loc = self.thick.mln + self.thick.uda + self.id * self.thick.s
		self.bound = False
	
	def __repr__(self):
		pass

	def bop(self):
		"""Knock an unbound XB around, return a new position """
		Fmag = rn.normalvariate(self.Fm,self.Fv)
		Fdir = rn.uniform(0, 2*pi)
		Fperp = Fmag * sin(self.Cs - Fdir)
		Fpara = Fmag * cos(self.Cs - Fdir)
		Rad = (Fpara / self.Gk) + self.Gs
		Phi = (Fperp * Rad / self.Ck) + self.Cs
		return (Rad * cos(Phi), Rad * sin(Phi))
	
	def head_loc(self, use_rest_values=False):
		"""Return the location of myosin head"""
		# Unpack component values (or their rest positions)
		if use_rest_values == False:
			[T,N,C,G]=[self.T, self.N, self.C, self.G]
		else:
			[T,N,C,G]=[self.R_T, self.R_N, self.R_C, self.R_G]
		# Calculate the x and y location of the head
		# y_T = 0 # this is assumed
		x_C = cos(T)*N
		y_C = sin(T)*N
		x_H = x_C + cos(C)*G
		y_H = y_C + sin(C)*G
		# Give back a tuple
		return (x_H, y_H)
	

class ThickFil:
	"""An instantiation of the 2D, 2fil system thick filament"""
	def __init__(self, mln=0, sep=12, thin_fil=None):
		self.k = 2020  # spring const between thick fil sites
		self.s = 43    # rest len between thick fil sites
		self.n = 20    # number of thick fil sites
		self.uda = 40  # dist from the m-line to first myosin head
		self.mln = mln # location of the m-line
		self.sep = sep # distance to thin fil
		self.thin = thin_fil # opposing filament
		if thin_fil is not None:
			self.myo = [XB(i, self, self.thin) for i in range(self.n)]
		else:
			print('Warning: ThickFil not given a thin_fil, XBs not created')


class ThinFil:
	def __init__(self, zln=1000, thick_fil=None):
		self.k = 1743  # spring const between thin fil sites
		self.s = 37.3  # rest length betwen thin fil sites
		self.n = 30    # number of thin fil sites
		self.zln = zln # location of the z-line
		self.thick = thick_fil # thick fil facing this thin fil

	def link_thick(self, thick_fil):
		"""An easy way to remind oneself to link to the thick fil"""
		self.thick = thick_fil


class FilPair:
	def _init_(self):
		thin = ThinFil()
		thick = ThickFil(thin_fil=thin)
		thin.link_thick(thick)
		
