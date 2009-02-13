# $filaments_and_heads
#     $Used to instantiate thick and thin filaments
###############################################################################
# Revisions
#    v 0.01/cdw     20090201  (Re)Created
#    v 0.02/cdw     20090209  Essential structures in, bop (re)written
#
# ToDo
#    u cdw          20090209  (Re)write the binding check function
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
		self.rest_head_loc() # Sets head_loc to rest location
		self.bound = False
		self.state = 1 # 1 is unbound, 2 is loosely, 3 is strongly
	
	def __repr__(self):
		pass

	def bop(self):
		"""Knock an unbound XB around, update and return head position """
		Fmag = rn.normalvariate(self.Fm,self.Fv)
		Fdir = rn.uniform(0, 2*pi)
		Fperp = Fmag * sin(self.Cs - Fdir)
		Fpara = Fmag * cos(self.Cs - Fdir)
		Rad = (Fpara / self.Gk) + self.Gs
		Phi = (Fperp * Rad / self.Ck) + self.Cs
		Xoffset = Rad * cos(Phi)
		Yoffset = Rad * sin(Phi)
		self.head_loc = (self.loc + Xoffset, Yoffset)
		return (Xoffset, Yoffset)

	def bind(self):
		"""Check if the XB binds, update if it does"""
		thin_loc = self.thin.closest_binding_site(self)
		dist = hypot(thin_loc - self.head_loc[0],
					 self.thick.sep - self.head_loc[1])
		if 1-(exp(-pow(dist,2))) < rn.random() and Af.bst(MinInd)==0:
			pass
		
		# Unfinished
	
	def rest_head_loc(self):
		"""Set the head loc to its rest location"""
		self.head_loc = (self.loc + self.Gs * cos(self.Cs),
						 self.Gs * sin(self.Cs))


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
		self.loc = ([self.mln + self.uda + (x * self.s)
					 for x in range(self.n)])
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
		self.loc = ([self.zln - (self.n - x) * self.s
					 for x in range(self.n)])
		self.bound = [False for x in range(self.n)]

	def link_thick(self, thick_fil):
		"""An easy way to remind oneself to link to the thick fil"""
		self.thick = thick_fil

	def closest_binding_site(self,XB):
		"""Return the closest binding site, but only if it is free"""
		# Make the above line accurate
		closest_val = min([act_loc - XB.head_loc for act_loc in self.loc])
		closest_ind = self.loc.index(closest_val)
		if self.bound[closest_ind] is False:
			return (0, 
		else:
			return (closest_val, closest_ind)
		

class FilPair:
	def __init__(self, hsl=1200):
		"""Creates and returns a set of two fils with the specified half-sl"""
		self.thin = ThinFil(hsl)
		self.thick = ThickFil(thin_fil=self.thin)
		self.thin.link_thick(self.thick)


# Test by creation of a FilPair
fp = FilPair(1100)
