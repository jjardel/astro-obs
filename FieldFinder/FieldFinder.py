import numpy as np
from scipy.spatial.distance import cdist
import astropy
from astropy.coordinates import ICRS
from astropy import units as u
import pysao

# Program to determine optimal positioning of the IFU when observing a
# resolved stellar system so as to put as many bright stars as possible on fibers.
#
# Edit the arguments in kwargs when called from main() to customize.
# =====================================================================
# set imageFile as path to image to display over
# set photFile as path to file with photometry data
#                (see getBrightStars method for format)
# set centerCoords to file containing RA and Dec in decimal degrees as
#                152.117081 12.306500
# set posThreshold to control strictness of proximity of stars to fibers
# set magThreshold to specify cut in stars below a luminosity
# =====================================================================
#
#    Usage: python FieldFinder.py
#
#
#    Copyright (C) John R. Jardel 2013
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. <http://www.gnu.org/licenses/>.



 

class FieldFinder:
    """
    Determines optimal IFU positioning.  Takes its arguments from kwargs in main.
    Most of the work is done in optimize, visualizations are created in displayFit.
    """

    def __init__( self, **kwargs ): 

        self.IFUoffsets = self.getIFUoffsets( **kwargs )
        self.starCoords = self.getBrightStarCoords( **kwargs )
        self.bestOffset = self.optimize( **kwargs )
        
        self.displayFit( **kwargs )

    def getIFUoffsets( self, **kwargs ):
        # read IFUfile for IFU offsets.  Should work with VP or VW
        
        filename = kwargs[ 'IFUfile' ]

        xIFU = []
        yIFU = []
        with open( filename ) as f:
            for line in f:
                xIFU.append( float( line.split()[ 1 ] ) )
                yIFU.append( float( line.split()[ 2 ] ) )

        return np.array( zip( xIFU, yIFU ) )

    def getBrightStarCoords( self, **kwargs ):
        # open file with photometry (photFile)
        # formatted as FLAG MAG RA DEC where FLAG is nonzero if not point source

        photFilename = kwargs[ 'photFile' ]
        centerCoordsFile = open( kwargs[ 'centerCoords' ] )
        threshold = kwargs[ 'magThreshold' ]

        line = centerCoordsFile.readline()
        self.centerCoords = np.array( [ float( line.split()[ 0 ] ),
                                        float( line.split()[ 1 ] ) ] )
        
        centerCoordsFile.close()

        starCoords = []
        with open( photFilename ) as f:
            for line in f:
                flag = int( line.split()[ 0 ] )
                mag = float( line.split()[ 1 ] )
                if flag == 0 and mag < threshold:
                    sCoords = line.split()[ 2 ] + ' ' + line.split()[ 3 ]
                    coords = ICRS( sCoords, unit = ( u.hour, u.degree ) )
                    starCoords.append( coords )
                    

        return starCoords # this is an ICRS object

    def countStars( self, testCoords, starCoords, **kwargs ):

        posThreshold = kwargs[ 'posThreshold' ]
        # factor to multiply a fiber diameter by to control how strictly
        # we determine whether a star is "on" a fiber
        # e.g. set to .5 for 50% of a fiber's diameter

        nStars = 0

        # calculate distances of all stars to all fibers
        dists = cdist( testCoords, starCoords )
        # count stars whose position is within a fiber diameter
    
        closeDists = np.where( dists < kwargs[ 'fiberDiam' ] / 3600. * posThreshold )
        nStars = len( set( closeDists[ 0 ] ) )

        return nStars
            
    def optimize( self, **kwargs ):

        yStepRange = kwargs[ 'yStepRange' ]
        xStepRange = kwargs[ 'xStepRange' ]
        nSteps = kwargs[ 'nSteps' ]

        xCoords = [ x.ra.degree for x in self.starCoords ]
        yCoords = [ x.dec.degree for x in self.starCoords ]

        nStarsMax = -1

        # main loop to shuffle the IFU position
        for i in np.linspace( xStepRange[ 0 ], xStepRange[ 1 ], num = nSteps ):
            for j in np.linspace( yStepRange[ 0 ], yStepRange[ 1 ], num = nSteps ):
                testCoords = np.copy( self.IFUoffsets )  / 3600. + self.centerCoords
                testCoords[ :, 0 ] += i / 3600.
                testCoords[ :, 1 ] += j / 3600.
                nStars = self.countStars( testCoords,
                                          zip( xCoords, yCoords ), **kwargs )
                if nStars > nStarsMax:
                    nStarsMax = nStars
                    xOffset = i
                    yOffset = j
        
        print nStarsMax, ' stars on fibers'
        print 'best offset is ', [ xOffset, yOffset ]

        # print final coordinates too
        finalRA = self.centerCoords[ 0 ] + xOffset / 3600.
        finalDec = self.centerCoords[ 1 ] + yOffset / 3600.

        print 'optimal pointing is ', [ finalRA, finalDec ]
        return np.array( [ xOffset, yOffset ] )

    def displayFit( self, **kwargs ):
        # get ready to show the IFU overlaid on the image provided by imageFile
        
        self.drawIFU( **kwargs )
        self.drawStars( **kwargs )

        # load region files with IFU and star coordinates
        ds9 = pysao.ds9()
        ds9.load_fits( kwargs[ 'imageFile' ] )
        ds9.set( 'regions load ifu.reg' )
        ds9.set( 'regions load stars.reg' )

        # pan the ds9 frame to center on the IFU
        xPan = self.bestOffset[ 0 ] / 3600. + self.centerCoords[ 0 ]
        yPan = self.bestOffset[ 1 ] / 3600. + self.centerCoords[ 1 ]

        coords = ICRS( ra = xPan, dec = yPan, unit = ( u.degree, u.degree ) )
        ds9.set( 'pan to %(x)s %(y)s wcs fk5' % { "x": coords.ra.to_string(),
                                                  "y": coords.dec.to_string() } )
        ds9.set( 'scale zscale' ) # set to zscale
        raw_input( "Press ENTER when finished" )


    
    def drawIFU( self, **kwargs ):
        # make region file for IFU

        # does this need a cos( dec ) correction?
        IFUcoords = ( self.IFUoffsets + self.bestOffset) / 3600. + self.centerCoords
        fRegion = open( 'ifu.reg', 'w' )
        for fiber in IFUcoords:
            line = 'fk5; circle( %(x)s,%(y)s,%(rad)s")\n' % {"x": fiber[ 0 ],
                                                           "y": fiber[ 1 ],
                                                           "rad": kwargs[ 'fiberDiam' ] / 2. }
            fRegion.write( line )
        fRegion.close()
        
    def drawStars( self, **kwargs ):
        # make region file for bright stars
        
        starCoords = self.starCoords
        fRegion = open( 'stars.reg', 'w' )
        for star in starCoords:
            line = 'fk5; diamond point %(x)s %(y)s # color=red\n' % { "x": star.ra.degree,
                                                                      "y": star.dec.degree }
            fRegion.write( line )
        fRegion.close()
        

def main( **kwargs ):

    pointing = FieldFinder( **kwargs )
    

if __name__ == '__main__':
    kwargs = { 'IFUfile': '/Users/jardel/research/lowdisp/leo1/IFUcen.txt', 
               'xStepRange': [ -100, 100 ], # max offset in arcsec
               'yStepRange': [ -100, 100 ],
               'nSteps': 200, # max number of steps
               'fiberDiam': 3.14, # in arcsec on sky
               'magThreshold': 20.5, # target stars less than this magnitude
               'posThreshold': 0.5, # proximity threshold
               'photFile': '/Users/jardel/research/lowdisp/leo1/stars.phot',
               'imageFile': '/Users/jardel/research/lowdisp/leo1/frame-g-003631-3-0238.fits.bz2',
               'centerCoords': '/Users/jardel/research/lowdisp/leo1/center.coords'
               }
    main( **kwargs )

    
    
