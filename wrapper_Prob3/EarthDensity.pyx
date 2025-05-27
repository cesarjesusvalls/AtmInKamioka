#cython .pyx for EarthDensity.h
# Import necessary declarations
from libc.string cimport memcpy
from libc.stdlib cimport malloc, free

cdef extern from "/Users/cjesus/Documents/Prob3plusplus/EarthDensity.h":
    cdef cppclass EarthDensity:
        EarthDensity()
        EarthDensity( const char * )
        void LoadDensityProfile( const char * )
        void SetDensityProfile( double, double, double )
        int get_LayersTraversed( )
        double get_DistanceAcrossLayer( int )
        double get_DensityInLayer( int )
        double get_YpInLayer( int )

cdef class PyEarthDensity:
    cdef EarthDensity* thisptr
    cdef object owner # keep track of Python-side owner of this object
    def __cinit__(self):    # Constructor to initialize the c++ EarthDensity object
        self.thisptr = new EarthDensity()
    def __init__(self, const_char_ptr=None):
        if const_char_ptr is not None:
            self.thisptr = new EarthDensity(const_char_ptr.encode('utf-8'))
    def __dealloc__(self):  # Destructor to delete the c++ EarthDensity object
        del self.thisptr
    #def set_earth(self, EarthDensity* ptr):    # set the c++ EarthDensity object
    #    self.thisptr = ptr
    #def get_earth(self):    # get the c++ EarthDensity object
    #    return PyEarthDensity().set_earth(self.thisptr.GetSelf())
    def loadDensityProfile(self, const_char_ptr:str):  # load density profile
        self.thisptr.LoadDensityProfile(const_char_ptr.encode('utf-8'))
    def setDensityProfile(self, cosZ, pathlength, prodHeight):  # set density profile
        self.thisptr.SetDensityProfile(cosZ, pathlength, prodHeight)
    def getLayersTraversed(self):  # get number of layers traversed
        return self.thisptr.get_LayersTraversed()
    def getDistanceAcrossLayer(self, layer:int):  # get distance across layer
        return self.thisptr.get_DistanceAcrossLayer(layer)
    def getDensityInLayer(self, layer:int):  # get density in layer
        return self.thisptr.get_DensityInLayer(layer)
    def getYpInLayer(self, layer:int):  # get Yp in layer
        return self.thisptr.get_YpInLayer(layer)