from Vec3_t import *
class Boundary:
    def __init__(self):
        self.inDensity	= 0.#/< Apply a certain density to inflow particles
        self.outDensity	= 0.#/< Apply a certain density to outflow particles
        self.allDensity	= 0.#/< Apply a certain density to outflow particles

        self.inv =Vec3_t(0.,0.,0.)		#/< Apply a certain velocity to inflow particles
        self.outv=Vec3_t(0.,0.,0.)		#/< Apply a certain velocity to outflow particles
        self.allv=Vec3_t(0.,0.,0.)		#/< Apply a certain velocity to all particle

        #TODO:UNCOMMENT
        # Periodic[3]	= False#/< Considering periodic in all directions => 0=X, 1=Y, 2=Z

        self.InOutFlow=0	#/< Considering inflow in all directions  by adding and deleting particles=> [0]=X, [1]=Y, [2]=Z and 0=none, 1=-
        # InFlowLoc1
        # InFlowLoc2
        # InFlowLoc3
        # OutFlowLoc
        # cellfac
        # inoutcounter
        # MassConservation

        # Array <>	OutPart
        # Array <>	InPart
