# .pyx for BargerPropagator.h
cdef extern from "/Users/cjesus/Documents/Prob3plusplus/BargerPropagator.h":
    cdef cppclass BargerPropagator:
        BargerPropagator()
        BargerPropagator( const char * )
        void propagate( int )
        void propagateLinear( int, double, double )
        void DefinePath( double, double, bool )
        void SetMNS( double, double, double, double, double, double, double, bool, int )
        void ResetMNS(double, int)
        void SetDensityConversion( double )
        double GetProb( int, int )
        double GetPathLength()
        void SetPathLength( double )
        void SetEnergy( double )
        void SetMatterPathLength( )
        void SetAirPathLength( double )
        void UseMassEigenstates( bool )
        void SetWarningSuppression( bool )
        void SetOneMassScaleMode( bool )
        double GetPathAveragedDensity( )
        void SetDefaultOctant( int, int )
cdef class PyBargerPropagator:
    cdef BargerPropagator* thisptr
    def __cinit__(self):    # Constructor to initialize the c++ BargerPropagator object
        self.thisptr = new BargerPropagator()
    def __init__(self, const_char_ptr=None):
        if const_char_ptr is not None:
            self.thisptr = new BargerPropagator(const_char_ptr.encode('utf-8'))
    def __dealloc__(self):  # Destructor to delete the c++ BargerPropagator object
        del self.thisptr
    def definePath(self, cosZ, prodHeight, kSetProfile:bool): # Define Path
        self.thisptr.DefinePath(cosZ, prodHeight, kSetProfile)
    def getProb(self, nuIn: int, nuOut: int):  # Get Osc Prob
        return self.thisptr.GetProb(nuIn, nuOut)
    def propagate(self, kNuType:int): # Propagate function
        self.thisptr.propagate(kNuType)
    def propagateLinear(self, kNuFlavor:int, pathLengh, density): # Propagate function
        self.thisptr.propagateLinear(kNuFlavor, pathLengh, density)
    def setMNS(self, theta12, theta13, theta23,
               dm21, dm32, deltacp,
               energy_, kSquared:bool, kNuType:int): # Set MNS
        self.thisptr.SetMNS(theta12, theta13, theta23, dm21, dm32, deltacp, energy_, kSquared, kNuType)