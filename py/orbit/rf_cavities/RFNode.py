"""
Module. Includes classes for RF accelerator nodes.
"""

import sys
import os
import math
from math import pi, sin
#from scipy.interpolate import interp1d

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode,\
AccActionsContainer, AccNodeBunchTracker

# import teapot drift class
from orbit.teapot import DriftTEAPOT

#import RF cavity classes
from rfcavities import Frequency_Cav
from rfcavities import Harmonic_Cav
from rfcavities import Barrier_Cav

class Base_RFNode(DriftTEAPOT):

	def __init__(self, length, name = "base_rfnode"):
		"""
			Constructor. Creates Base RF Cavity TEAPOT element.
			It will never be called.
		"""
		DriftTEAPOT.__init__(self, name)
		self.setType("base rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		#self..trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		#self..trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class Frequency_RFNode(Base_RFNode):

	def __init__(self, RFFreq, RFE0TL, RFPhase,\
		length, name = "frequency_rfnode"):
		"""
			Constructor. Creates Frequency
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.frequencynode = Frequency_Cav(RFFreq, RFE0TL, RFPhase)
		self.setType("frequency rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.frequencynode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.frequencynode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase,\
		length, name = "harmonic_rfnode"):
		"""
			Constructor. Creates Harmonic
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			RFHNum, RFVoltage, RFPhase)
		self.setType("harmonic rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class Dual_Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, RFHNum, RatioRFHNum, RFVoltage, RatioVoltage, RFPhase, RFPhase2,\
		length, name = "harmonic_rfnode"):
		"""
			Constructor. Creates Harmonic
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.dualharmonicnode = Dual_Harmonic_Cav(ZtoPhi, \
			RFHNum, RatioRFHNum, RFVoltage, RatioVoltage, RFPhase, RFPhase2)
		self.setType("dual harmonic rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.dualharmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.dualharmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class BRhoDep_Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, accelDict, bunch,\
	        name = "brho_timedep_harmonic_rfnode"):
		"""
			Constructor. Creates BRho Time Dependent
			Harmonic RF Cavity TEAPOT element
		"""
                Base_RFNode.__init__(self, length=0.0, name=name)
		self.Z2Phi     = ZtoPhi
		self.localDict = accelDict

                gammaTrans      = self.localDict["gammaTrans"]
                RFHNum          = self.localDict["RFHNum"]

                BRho_vs_t       = self.localDict["BRho_vs_t"]
                RFVoltage_vs_t  = self.localDict["RFVoltage_vs_t"]
                RFPhase_vs_t    = self.localDict["RFPhase_vs_t"]
                phaseCor        = self.localDict["RFJumpCor"]                

                #print("phaseCor = {:g}".format(phaseCor))

                time = bunch.getSyncParticle().time()
		mass = bunch.mass()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()

		BRho      = BRho_vs_t(time)
                RFVoltage = RFVoltage_vs_t(time)
		RFPhase   = RFPhase_vs_t(time)

                
                pcnew = 0.299792458 * charge * BRho
		bunch.getSyncParticle().momentum(pcnew)
                dESync  = charge*RFVoltage*sin(RFPhase*pi/180.0) # is this needed at all ? 

		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			                         RFHNum, RFVoltage, RFPhase, gammaTrans, phaseCor)

		self.setType("harmonic rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())

                gammaTrans      = self.localDict["gammaTrans"]
                RFHNum          = self.localDict["RFHNum"]
                BRho_vs_t       = self.localDict["BRho_vs_t"]
                RFVoltage_vs_t  = self.localDict["RFVoltage_vs_t"]
                RFPhase_vs_t    = self.localDict["RFPhase_vs_t"]
                
		time = bunch.getSyncParticle().time()
		mass = bunch.mass()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()

		BRho      = BRho_vs_t(time)
                RFVoltage = RFVoltage_vs_t(time)
		RFPhase   = RFPhase_vs_t(time)

                pcnew = 0.299792458 * charge * BRho
		bunch.getSyncParticle().momentum(pcnew)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)

		#put the track method here:
                self.harmonicnode.trackBunch(bunch)

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length     = self.getLength(self.getActivePartIndex())
		bunch      = paramsDict["bunch"]

                gammaTrans      = self.localDict["gammaTrans"]
		RFHNum          = self.localDict["RFHNum"]
                BRho_vs_t       = self.localDict["BRho_vs_t"]        # callable - interpolation function 
                RFVoltage_vs_t  = self.localDict["RFVoltage_vs_t"]   # callable - interpolation function
                RFPhase_vs_t    = self.localDict["RFPhase_vs_t"]     # callable - interpolation function 

                time = bunch.getSyncParticle().time()
		mass = bunch.mass()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
	

                BRho      = BRho_vs_t(time)
                RFVoltage = RFVoltage_vs_t(time)
		RFPhase   = RFPhase_vs_t(time)

                pcnew = 0.299792458 * charge * BRho
                bunch.getSyncParticle().momentum(pcnew) 

                ZtoPhi = self.Z2Phi
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)

		#put the track method here:
		self.harmonicnode.trackBunch(bunch)



class SyncPhaseDep_Harmonic_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, accelDict, bunch,\
		length, name = "syncphase_timedep_harmonic_rfnode"):
		"""
			Constructor. Creates SyncPhase Time Dependent
			Harmonic RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.Z2Phi = ZtoPhi
		self.localDict = accelDict
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode = Harmonic_Cav(ZtoPhi, dESync,\
			RFHNum, RFVoltage, RFPhase)
		self.setType("harmonic rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		ZtoPhi = self.Z2Phi
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		gammaTrans = self.localDict["gammaTrans"]
		RFHNum = self.localDict["RFHNum"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		SyncPhase_tuple = self.localDict["SyncPhase"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhase_tuple = self.localDict["RFPhase"]
		time = bunch.getSyncParticle().time()
		charge = bunch.charge()
		gamma = bunch.getSyncParticle().gamma()
		keold = bunch.getSyncParticle().kinEnergy()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhase = interp(time, n_tuple,\
			time_tuple, RFPhase_tuple)
		SyncPhase = interp(time, n_tuple,\
			time_tuple, SyncPhase_tuple)
		dESync = charge * RFVoltage *\
			math.sin(math.pi *\
			(RFHNum * SyncPhase + RFPhase) / 180.0)
		kenew = keold + dESync
		bunch.getSyncParticle().kinEnergy(kenew)
		ZtoPhi = self.Z2Phi
		Zsync = - math.pi * SyncPhase / (180.0 * ZtoPhi)
		bunch.getSyncParticle().z(Zsync)
		self.harmonicnode.dESync(dESync)
		self.harmonicnode.RFVoltage(RFVoltage)
		self.harmonicnode.RFPhase(RFPhase)
		#put the track method here:
		self.harmonicnode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class Barrier_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, RFVoltage, RFPhasep,\
		RFPhasem, dRFPhasep, dRFPhasem,\
		length, name = "barrier_rfnode"):
		"""
			Constructor. Creates Barrier
			RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.barriernode = Barrier_Cav(ZtoPhi, RFVoltage,\
			RFPhasep, RFPhasem, dRFPhasep, dRFPhasem)
		self.setType("barrier rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		#put the track method here:
		self.barriernode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		#put the track method here:
		self.barriernode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

class TimeDep_Barrier_RFNode(Base_RFNode):

	def __init__(self, ZtoPhi, accelDict, bunch,\
		length, name = "syncphase_timedep_harmonic_rfnode"):
		"""
			Constructor. Creates SyncPhase Time Dependent
			Harmonic RF Cavity TEAPOT element
		"""
		Base_RFNode.__init__(self, length, name)
		self.localDict = accelDict
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhasep_tuple  = self.localDict["RFPhasep"]
		RFPhasem_tuple  = self.localDict["RFPhasem"]
		dRFPhasep_tuple = self.localDict["dRFPhasep"]
		dRFPhasem_tuple = self.localDict["dRFPhasem"]
		time = bunch.getSyncParticle().time()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhasep = interp(time, n_tuple,\
			time_tuple, RFPhasep_tuple)
		RFPhasem = interp(time, n_tuple,\
			time_tuple, RFPhasem_tuple)
		dRFPhasep = interp(time, n_tuple,\
			time_tuple, dRFPhasep_tuple)
		dRFPhasem = interp(time, n_tuple,\
			time_tuple, dRFPhasem_tuple)
		self.barriernode = Barrier_Cav(ZtoPhi, RFVoltage,\
			RFPhasep, RFPhasem, dRFPhasep, dRFPhasem)
		self.setType("barrier rf node")
		self.setLength(0.0)

	def trackBunch(self, bunch):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhasep_tuple  = self.localDict["RFPhasep"]
		RFPhasem_tuple  = self.localDict["RFPhasem"]
		dRFPhasep_tuple = self.localDict["dRFPhasep"]
		dRFPhasem_tuple = self.localDict["dRFPhasem"]
		time = bunch.getSyncParticle().time()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhasep = interp(time, n_tuple,\
			time_tuple, RFPhasep_tuple)
		RFPhasem = interp(time, n_tuple,\
			time_tuple, RFPhasem_tuple)
		dRFPhasep = interp(time, n_tuple,\
			time_tuple, dRFPhasep_tuple)
		dRFPhasem = interp(time, n_tuple,\
			time_tuple, dRFPhasem_tuple)
		self.barriernode.RFVoltage(RFVoltage)
		self.barriernode.RFPhasep(RFPhasep)
		self.barriernode.RFPhasem(RFPhasem)
		self.barriernode.dRFPhasep(dRFPhasep)
		self.barriernode.dRFPhasem(dRFPhasem)
		#put the track method here:
		self.barriernode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length

	def track(self, paramsDict):
		"""
			The rfcavity-teapot class implementation of the
			AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		n_tuple = self.localDict["n_tuple"]
		time_tuple = self.localDict["time"]
		RFVoltage_tuple = self.localDict["RFVoltage"]
		RFPhasep_tuple  = self.localDict["RFPhasep"]
		RFPhasem_tuple  = self.localDict["RFPhasem"]
		dRFPhasep_tuple = self.localDict["dRFPhasep"]
		dRFPhasem_tuple = self.localDict["dRFPhasem"]
		time = bunch.getSyncParticle().time()
		RFVoltage = interp(time, n_tuple,\
			time_tuple, RFVoltage_tuple)
		RFPhasep = interp(time, n_tuple,\
			time_tuple, RFPhasep_tuple)
		RFPhasem = interp(time, n_tuple,\
			time_tuple, RFPhasem_tuple)
		dRFPhasep = interp(time, n_tuple,\
			time_tuple, dRFPhasep_tuple)
		dRFPhasem = interp(time, n_tuple,\
			time_tuple, dRFPhasem_tuple)
		self.barriernode.RFVoltage(RFVoltage)
		self.barriernode.RFPhasep(RFPhasep)
		self.barriernode.RFPhasem(RFPhasem)
		self.barriernode.dRFPhasep(dRFPhasep)
		self.barriernode.dRFPhasem(dRFPhasem)
		#put the track method here:
		self.barriernode.trackBunch(bunch)
		#print "debug tracking the bunch through the rf node = ",\
		self.getName(), " part ind = ", self.getActivePartIndex(),\
		" length = ", length
		#print "time = ", time
		#print "RFVoltage = ", RFVoltage
		#print "RFPhaseplus = ", RFPhasep + dRFPhasep, RFPhasep - dRFPhasep
		#print "RFPhasemnus = ", RFPhasem - dRFPhasem, RFPhasem + dRFPhasem

def interp(x, n_tuple, x_tuple, y_tuple):
	"""
	Linear interpolation: Given n-tuple + 1 points,
	x_tuple and y_tuple, routine finds y = y_tuple
	at x in x_tuple. Assumes x_tuple is increasing array.
	"""
	if x <= x_tuple[0]:
		y = y_tuple[0]
		return y
	if x >= x_tuple[n_tuple]:
		y = y_tuple[n_tuple]
		return y
	dxp = x - x_tuple[0]
	for n in range(n_tuple):
		dxm = dxp
		dxp = x - x_tuple[n + 1]
		dxmp = dxm * dxp
		if dxmp <= 0:
			break
	y = (-dxp * y_tuple[n] + dxm * y_tuple[n + 1]) /\
		(dxm - dxp)
	return y

def syncZ(ZtoPhi, gammaTrans, gamma, charge,\

        # JFO: as far as I can tell, this function serves no purpose except perhaps as a diagnostic
        # for the sync particle, unless the voltage/rf phase are inconsistent for synchronism   
        # syncZ always return a very small number, nearly zero. If the synchronization is perfect,
        # this would be exactly 0. ???? The z coordinate of the syncParticle does not appear to be
        # used anywhere else in the code.
          
        dESync, RFHNum, RFVoltage, RFPhase):
	"""
	Calculates position of synchronous particle.
	"""
        Zsync = 0
	if abs(dESync) > abs(charge * RFVoltage): # this cannot happen: -1 < sin(phi) < 1 so one must have abs(dE/V) < 1    
		return Zsync
	PhaseTot = math.asin(dESync / (charge * RFVoltage))
	if gamma > gammaTrans and gammaTrans > 0:
		PhaseTot = math.pi - PhaseTot
	Zsync = -(PhaseTot - math.pi * RFPhase / 180.0)/ (RFHNum * ZtoPhi)
        return Zsync

