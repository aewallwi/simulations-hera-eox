from farfieldData import FarFieldData
import numpy as np

beamDir='../data/beams/varyRiggingHeight_short/beamVaryRiggingHeight_run000003/'
beamDirFeedOnly='../data/beams/beamCosmicTwilightPolarimeter_diskConnect_short/'
beamFile=beamDir+'farfield (f=%s) [1].txt'%('100')
beamFileFeedOnly=beamDirFeedOnly+'farfield (f=%s) [1].txt'%('100')
rigging_height=4.90
farfield=FarFieldData([beamFileFeedOnly],[beamFile],
                      ['100'],64,pols=['X'],rotateY=True,
                      invertFeedOnly=False,
                      rotatexz=True,
                      dFocus=rigging_height)
farfield.beamFeedAndDish.export_fits_prisim('../data/simulations/test.fits',['X'],[100e6])
