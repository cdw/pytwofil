# $filaments_and_heads
#     $Used to instantiate thick and thin filaments
###############################################################################
# Revisions
#    v 0.01/cdw     20090201  (Re)Created
#    v 0.02/cdw     20090209  Essential structures in, bop (re)written
#
# ToDo
#    u cdw          20090209  (Re)write the binding check function
#    u cdw          20090224  Fix the 'dist' in trans_loosely etc.
#
###############################################################################

import random as rn
from math import *
from scipy.optimize import fmin

class XB:
	"""A crossbridge that will be instantiated by the thick filament"""
	def __init__(self, head_id, thick_fil, thin_fil):
		# The first thing to do is define our default parameters.
		self.Cs = pi/3 # rest angle of converter domain
		self.Ck = 200  # torsional spring const of converter domain
		self.Cv = (pi/3, pi/3, 1.2*pi/3) # normal and rigor values of Cs
		self.Gs = 10.5 # rest length of globular domain
		self.Gk = 5    # spring constant of globular domain
		self.Gv = (5, 5, 5) # normal and rigor values of Gs
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
		self.state = 0 # 0 is unbound, 1 is loosely, 2 is strongly

	def __repr__(self):
		return ('XB #'+str(self.id)+' Info\n===========\nLoc:   '+str(self.loc)+
			'\nBound: '+str(self.bound)+'\nState: '+str(self.state))

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
		thin_loc = self.thin.closest_binding_site(self) # get closest thin site
		# thin_loc = (closest_val, closest_ind) or (False, False)
		if thin_loc is False: # if that site is alrady occupied
			return # there's no binding
		dist = dist_to_thin(thin_loc)
		if 1 - (exp(-pow(dist, 1/self.Bd))) < rn.random():
			# Where we are close enough, bind
			self.bound = thin_loc  # Where XB binds to 
			self.state = 1 # Now in the loosely bound state
			self.thin.bound[thin_loc] = self.id # Thin fil link back
	
	def trans_loosely(self):
		"""When in the loosely bound state: unbind, strongly bind, or pass"""
		dist = glob_len()
		angl = conv_ang()
		random = rn.random()
		# possibly bind tightly or unbind
		if 1-.001*100/sqrt(2*0.2515)*(1-tanh(1*sqrt(2*0.2515)*(dist-6))) < random:
			self.state = 2
			self.Cs = self.Cv[2]
		elif (pow(exp(15 * self.Gk * 0.2515 * 
				 pow(dist - self.Gs,2) + 
				 3 * self.Cs * 0.2515 / dist *
				 (angl - self.Cs)),2)  < random):
			self.thin.bound[self.bound] = False
			self.bound = False
			self.state = 0
	
	def trans_tightly(self):
		"""When in the tightly bound state: unbind, loosely bind, or pass"""
		dist = glob_len()
		angl = conv_ang()		
		random = rn.random()
		# possibly unbind or become loosely bound
		if (0.005 * 1 / exp(self.Gk * 0.2515 *pow(dist-self.Gv[1],2) +
				    self.Ck * 0.2515 / dist *pow(angl-self.Cv[1],2) -
				    self.Gk * 0.2515 *pow(dist-self.Gv[2],2) -
				    self.Ck * 0.2515 / dist *pow(angl-self.Cv[2],2)
				)) < random:
			self.state = 1
			self.Cs = self.Cv[1]
		elif (.001*(sqrt(self.Gk*0.2515)*
			    (sqrt(3600*pow(dist,2)) - 40*dist)+ 20)) > random:
			self.thin.bound[self.bound] = False
			self.bound = False
			self.state = 0
			self.Cs = self.Cv[0]
	
	def dist_to_thin(self, thin_point):
		return hypot(self.thin.loc[thin_point] - self.head_loc[0],
					 self.thick.sep - self.head_loc[1])
	
	def glob_len(self):
		"""Return the globular length for a bound XB"""
		return hypot(self.thin.loc[self.bound] - self.loc, self.thick.sep)
	
	def conv_ang(self):
		"""Return the converter angle for a bound XB"""
		return atan2(self.thin.loc[self.bound] - self.loc,self.thick.sep)
	
	def rest_head_loc(self):
		"""Set the head loc to its rest location"""
		self.head_loc = (self.loc + self.Gs * cos(self.Cs),
						 self.Gs * sin(self.Cs))


	def transition(self):
		"""Transition to a new state, maybe"""
		if True:
			pass


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
	
	def __repr__(self):
		o_loc = ([self.mln + self.uda + (x * self.s) for x in range(self.n)])
		return ('Thick Fil Info\n==============\nLoc Offset: ' + 
			str([self.loc[i] - o_loc[i] for i in range(len(o_loc))])+
			'\nBound: '+str([m.bound for m in self.myo])+
			'\nState: '+str([m.state for m in self.myo]))


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
	
	def __repr__(self):
		o_loc = ([self.zln - (self.n - x) * self.s for x in range(self.n)])
		return ('Thin Fil Info\n=============\nLoc Offset: '+
			str([self.loc[i] - o_loc[i] for i in range(len(o_loc))])+
			'\nBound: '+str(self.bound))

	def link_thick(self, thick_fil):
		"""An easy way to remind oneself to link to the thick fil"""
		self.thick = thick_fil

	def closest_binding_site(self,XB):
		"""Return the closest binding site, but only if it is free"""
		closest_val = min([act_loc - XB.head_loc for act_loc in self.loc])
		closest_ind = self.loc.index(closest_val)
		if self.bound[closest_ind] is False:
			return closest_ind
		else:
			return False


class FilPair:
	def __init__(self, hsl=1200):
		"""Creates and returns a set of two fils with the specified half-sl"""
		self.thin = ThinFil(hsl)
		self.thick = ThickFil(thin_fil=self.thin)
		self.thin.link_thick(self.thick)
		
	def settle(self):
		"""Balance the forces felt by the two filaments"""
		# Create our initial guess, which is just the current
		# location of the nodes along the thick and thin filaments
		x0 = self.thick.loc + self.thin.loc
		# Unpack some variables for ease of writing the force function
		Mk  = self.thick.k
		Ms  = self.thick.s
		Mu = self.thick.uda;
		Ml = [m.bound for m in self.thick.myo]
		Mb = [m.state for m in self.thick.myo]
		Gk  = [m.Gk for m in self.thick.myo]
		Gs  = [m.Gs for m in self.thick.myo]
		Ck = [m.Ck for m in self.thick.myo]
		Cs = [m.Cs for m in self.thick.myo]
		Ak  = self.thin.k
		As  = self.thin.s
		Al = self.thin.bound
		
		def  force(self, x):
			""" Return a matrix of the forces on all points of a two 
			filament system based on the locations fed in. It is
			good to note that the location inputs are in the form:
			[ThickXLoc1, ThickXLoc2,...,ThinXLoc1,ThinXLoc2...]"""
			
		# Settle and subfunctions are unfinished


# Test by creation of a FilPair
fp = FilPair(1100)
