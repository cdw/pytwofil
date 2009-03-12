# $filaments_and_heads
#     $Used to instantiate thick and thin filaments
###############################################################################
# Revisions
#    v 0.01/cdw     20090201  (Re)Created
#    v 0.02/cdw     20090209  Essential structures in, bop (re)written
#    v 0.03/cdw     20090306  Believed to be complete recreation of Matlab machinery
#
# ToDo
#    u cdw          20090209  (Re)write the binding check function
#    u cdw          20090224  Fix the 'dist' in trans_loosely etc.
#
###############################################################################

import random as rn
import numpy as np
from numpy import array,cos,sin,arctan2,hypot,pi
from scipy.optimize import fsolve

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
		return arctan2(self.thin.loc[self.bound] - self.loc,self.thick.sep)
	
	def rest_head_loc(self):
		"""Set the head loc to its rest location"""
		self.head_loc = (self.loc + self.Gs * cos(self.Cs),
						 self.Gs * sin(self.Cs))


	def transition(self, give_detail=0):
		"""Transition to a new state, maybe
		
		This checks the XB for its current state and 
		checks if it transitions to a new state based on the 
		rate constants. It directly sets the new state of the 
		XB and returns a parameter describing what occured.
		
		:Input:
		
		give_detail: int
			0 - just return the transition flag (default)
			1 - give back all the optional output information
		
		:Output:
		
		flag: boolean
			False - no transition occured
			True - some transition occured
		"""
		state = self.state
		if state == 0:
			self.bind()
		elif state == 1:
			self.trans_loosely
		elif state ==2:
			self.trans_tightly()
		flag = state != self.state
		if give_detail == 0:
			return flag
		else:
			# TODO create a returned data structure for use 
			# when more detail is needed 
			return flag


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
		self.loc = array([self.mln + self.uda + (x * self.s)
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
	def __init__(self, zln=1200, thick_fil=None):
		self.k = 1743  # spring const between thin fil sites
		self.s = 37.3  # rest length betwen thin fil sites
		self.n = 30    # number of thin fil sites
		self.zln = zln # location of the z-line
		self.thick = thick_fil # thick fil facing this thin fil
		self.loc = array([self.zln - (self.n - x) * self.s
					 for x in range(self.n)])
		if self.loc[0] < 0:
			print('Warning: Thin Fil error, fil is too short')
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
		closest_val = np.argmin(np.absolute(self.loc - XB.head_loc[0] ))
		# FIXME I think this works but it needs checking
		closest_ind = np.searchsorted(self.loc, closest_val)
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
		
	def settle(self, give_detail=0):
		""" Relax positions to balance forces.
		
		Balance the forces felt by the two filaments by solving
		for roots of the forces at each node. Optional parameter,
		'give_detail' returns info about the optimization if set 
		to other than zero.
		
		:Input:
		
		give_detail : number
			0 - return only the new location vector (default)
			1 -  return all output arguments (x1, infodict, ier, mesg)
			
		:Output:
		
		x1: numpy array
			the solution (or the result of last iteration)
		infodict : a dictionary of optional outputs with the keys
			see documentation of scipy.optimize.fsolve
		ier: an integer flag equal to 1 when a solution is found
		mesg: a string message about the cause of failure, should
			fsolve fail to produce a solution
		
		"""
		# Create our initial guess, which is just the current
		# location of the nodes along the thick and thin filaments
		x0 = np.hstack((self.thick.loc, self.thin.loc))
		# Unpack some variables for ease of writing the force function
		Mk = self.thick.k # myosin spring const
		Ms = self.thick.s # myosin spring length
		Mu = self.thick.uda #  thick fil undecorated length
		Ml = array([m.bound for m in self.thick.myo]) # linkers when bound
		Mb = array([m.state for m in self.thick.myo]) # bound state
		Mn = len(self.thick.myo) # number of myosins
		Gk  = array([m.Gk for m in self.thick.myo]) # glob spring consts
		Gs  = array([m.Gs for m in self.thick.myo]) # glob spring lengths
		Ck = array([m.Ck for m in self.thick.myo]) # converter spring const
		Cs = array([m.Cs for m in self.thick.myo]) # converter rest angle
		Ak  = self.thin.k # actin spring constant
		As  = self.thin.s # actin spring rest length
		Al = self.thin.bound # actin linkers
		An = len(self.thin.loc) # number of actin sites
		Sep = self.thick.sep # Separation between thick and thin fils
		Zln = self.thin.zln # The end of the  1/2 sarc, the z-line
		
		def  force(x):
			""" Return a matrix of the forces on all points of a two 
			filament system based on the locations fed in. It is
			good to note that the location inputs are in the form:
			[ThickXLoc0, ThickXLoc1,...,ThinXLoc0,ThinXLoc1...]
			and that the outputs are in the form:
			[FThickXLoc0, FThickXLoc1,... and so on]"""
			# Initialize the force return array
			f = np.ones_like(x)
			## Begin with the forces on the thick filaments
			# First thick filament site
			f[0] = Mk*(x[1]-x[0]-Ms) - Mk*(x[0]-0-Mu) # force from adj site
			if Ml[0] != False:
				G = hypot(x[Mn+Ml[0]] - x[0], Sep)
				C = arctan2(Sep, x[Mn+Ml[0]] - x[0])
				f[0] = (f[0] +  
					Gk[0] * (G-Gs[0]) * cos(C) - 
					(1 / G) * Ck[0] * (C-Cs[0]) * sin(C))
			# Most thick filament sites
			for i in range(1, Mn-1):
				f[i]=Mk*(x[i+1]-x[i]-Ms) - Mk*(x[i]-x[i-1]-Ms)
				if Ml[i] != False:
					G = hypot(x[Mn+Ml[i]] - x[i], Sep)
					C = arctan2(Sep, x[Mn+Ml[i]] - x[i])
					f[i] = (f[i] +  
						Gk[i] * (G-Gs[i]) * cos(C) - 
						(1 / G) * Ck[i] * (C-Cs[i]) * sin(C))
			# Last thick filament site
			f[Mn-1] = - Mk*(x[Mn-1]-x[Mn-1-1]-Ms)
			if Ml[Mn-1] != False:
				G = hypot(x[Mn+Ml[Mn-1]] - x[Mn-1], Sep)
				C = arctan2(Sep, x[Mn+Ml[Mn-1]] - x[Mn-1])
				f[Mn-1] = (f[Mn-1] +  
					Gk[Mn-1] * (G-Gs[Mn-1]) * cos(C) - 
					(1 / G) * Ck[Mn-1] * (C-Cs[Mn-1]) * sin(C))
			## Continue to the forces on the thin filament
			# First thin filament site
			f[Mn] = Ak*(x[Mn+1]-x[Mn]-As)
			if Al[0] != False:
				G = hypot(x[Mn+0]-x[Al[0]], Sep)
				C = arctan2(Sep, x[Mn+0]-x[Al[0]])
				f[Mn] = (f[Mn] - 
					Gk[Al[0]] * (G-Gs[Al[0]]) * cos(C) + 
					(1 / G) * Ck[Al[0]] * (C-Cs[Al[0]]) * sin(C))
			# Most thin filament sites
			for i in range(Mn+1, Mn+An-1):
				f[i] = Ak*(x[i+1]-x[i]-As) - Ak*(x[i]-x[i-1]-As)
				if Al[i-Mn] != False:
					G= hypot(x[i]-x[Al[i-Mn]], Sep)
					C = arctan2(Sep, x[i]-x[Al[i-Mn]])
					f[i] = (f[i] - 
						Gk[Al[i-Mn]] * (G-Gs[Al[i-Mn]]) * cos(C) + 
						(1 / G) * Ck[Al[i-Mn]] * (C-Cs[Al[i-Mn]]) * sin(C))
			# Last thin filament site
			f[Mn+An-1] = Ak*(Zln-x[Mn+An-1]-As) - Ak*(x[Mn+An-1]-x[Mn+An-2]-As)
			if Al[An-1] != False:
				G = hypot(x[Mn+An-1]-x[Al[An-1]], Sep)
				C = arctan2(Sep, x[Mn+An-1]-x[Al[An-1]])
				f[Mn+An-1] = (f[Mn+An-1] - 
					Gk[Al[An-1]] * (G-Gs[Al[An-1]]) * cos(C) + 
					(1 / G) * Ck[Al[An-1]] * (C-Cs[Al[An-1]]) * sin(C))
			return f
			
		## Optimize with fsolve and update filament  with new locations
		x1 = fsolve(force, x0,full_output=give_detail)
		self.thick.loc = x1[0:Mn]
		self.thin.loc = x1[Mn:Mn+An]
		# Return our information
		return (x1)
	
	def step(self):
		"""Take a single time step.
		
		In other words, this applies all the processes that
		need to be gone through in order to update the 
		model through one time step; including knocking 
		the XBs around and balancing the forces felt by 
		the nodes along the filaments.
		"""
		trans = [m.transition() for m in self.thick.myo]
		
