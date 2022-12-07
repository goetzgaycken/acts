#!/usr/bin/env python3
import os
import pathlib, acts, acts.examples, acts.examples.itk

u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir,
                                                                            logLevel=acts.logging.WARNING,
)
print ("buildITkGeometry done")

def readAthenaMeasurements(trackingGeometry, inputFile, outputDir, s=None):
    # hepmc_dir = os.path.join(outputDir, "hepmc3")
    # if not os.path.exists(hepmc_dir):
    #     os.mkdir(hepmc_dir)

    # s = acts.examples.Sequencer(events=1, numThreads=-1, outputDir=str(outputDir))
    s = acts.examples.Sequencer(
        events=int(os.environ.get("NEVENTS", 5)), numThreads=-11,
        outputDir=str(outputDir),
    )

    # field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
    digiReader = acts.examples.RootDigiReader (
        level=acts.logging.DEBUG,
        # name of the particles container which is populated
        particles = 'particles_final',
        #  name of the sourceLink collection (Reader)
        sourceLinks = 'sourcelinks',
        # name of the measurements collection which is populated
        measurements = 'measurements',
        #  name for the map between measurements and truth particles which is populated
        measurementParticlesMap = 'measurement_particles_map',
        # Path to the input file, which is a "special" AOD with xAOD prep raw data obejcts and has a split-level of 99
        filePath = inputFile,
        #  Name of the tree within the output file.
        treeName = 'CollectionTree',

        # truth particle branch prefix
        truthParticlePrefix = 'TruthParticlesAux.',
        # truth vertex branch prefix
        truthVertexPrefix = 'TruthVerticesAux.',
        # truth vertex branch prefix
        eventInfoPrefix = 'EventInfoAux.',
        # truth  branch prefix
        measurementPrefix = ['ITkPixelClustersAux.','ITkStripClustersAux.'],
        # for each measurement prefix the type of the measurement : 0 pixel, 1 strips 
        measurementType = [0,1],

        trackingGeometry = trackingGeometry,

        outputFilePath = outputDir +'/geoMatching_EC_mu200_fullTruth.root',

        surfaceCenterFile = '/tmp/surface_centers.txt',
        # when set to true, prrogram will be stopped upon construction to allow for a debugger to be attached
        stop = False
    )

    s.addReader(digiReader)

    return s


if "__main__" == __name__:

    # trackingGeometry=None

#    readAthenaMeasurements(trackingGeometry, inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.pool.root', outputDir=os.getcwd()).run()
    readAthenaMeasurements(trackingGeometry, inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.mu200_fullTruth.pool.root', outputDir=os.getcwd()).run()
