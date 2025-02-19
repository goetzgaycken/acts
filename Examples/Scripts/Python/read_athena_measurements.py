#!/usr/bin/env python3
import os
import pathlib, acts, acts.examples, acts.examples.itk

from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    TruthSeedRanges,
    addCKFTracks,
    CKFPerformanceConfig,
    TrackSelectorRanges,
)

ttbar_pu200 = True

u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"

if False :
    detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir,
                                                                            logLevel=acts.logging.WARNING,
                                                                            )


config = acts.MaterialMapJsonConverter.Config()
customLogLevel = acts.examples.defaultLogging(logLevel=acts.logging.INFO)
# mdecorator = acts.JsonMaterialDecorator(
#     rConfig=config,
#    level=customLogLevel(minLevel=acts.logging.WARNING),
#    jFileName="acts-ws/build-sft/geometry-maps-volbounds_2023_material.json",
# )

jsonTGReaderCfg =     acts.examples.TrackingGeometryJsonReader.Config(detectorName="ITK",
                                                    toolLogLevel = acts.logging.VERBOSE,
                                                    logLevel = acts.logging.VERBOSE,
                                                   )
jsonTGReader=acts.examples.TrackingGeometryJsonReader(jsonTGReaderCfg)
# trackingGeometry = jsonTGReader.read("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps-volbounds_2023.json", mdecorator)
# trackingGeometry = jsonTGReader.read("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps_track_reco_default.json", mdecorator)
# trackingGeometry = jsonTGReader.read("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps_track_reco.json",mdecorator)
# material provided by the same json file
trackingGeometry = jsonTGReader.read("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps_track_reco.json", None)
# trackingGeometry = jsonTGReader.read("/data/goetz/src/ACTS/source/tg_dup.json", None)
# trackingGeometry = jsonTGReader.read("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps_track_reco.json::/data/goetz/src/ACTS/source/tg_dup.json", None)

jsonTGWriterCfg = acts.examples.TrackingGeometryJsonWriter.Config(
    logLevel = acts.logging.VERBOSE,
    processSensitives=True,
    processApproaches=True,
    processRepresenting=True,
    processBoundaries=True,
    processVolumes=True,
    processNonMaterial = True)

jsonTGWriter=acts.examples.TrackingGeometryJsonWriter(jsonTGWriterCfg)
jsonTGWriter.writeNominal(trackingGeometry,"tg_dup2.json")

print ("buildITkGeometry done")

def readAthenaMeasurements(trackingGeometry, inputFile, outputDir, s=None):
    # hepmc_dir = os.path.join(outputDir, "hepmc3")
    # if not os.path.exists(hepmc_dir):
    #     os.mkdir(hepmc_dir)

    # s = acts.examples.Sequencer(events=1, numThreads=-1, outputDir=str(outputDir))
    s = acts.examples.Sequencer(
        events=int(os.environ.get("NEVENTS", 10)), numThreads=1,
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
        stop = False,

        singleParticle = True
    )

    s.addReader(digiReader)

    field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))

    addSeeding(
        s,
        trackingGeometry,
        field,
        TruthSeedRanges(pt=(.1 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
        if ttbar_pu200
        else TruthSeedRanges(),
        *acts.examples.itk.itkSeedingAlgConfig("PixelSpacePoints"),
        #    seedingAlgorithm = SeedingAlgorithm.Orthogonal,
        inputParticles = "particles_final",
        seedingAlgorithm = SeedingAlgorithm.Default,
        # geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
        # geoSelectionConfigFile="geoSelection-athenaITk.json",
        geoSelectionConfigFile="geoSelection-athenaITkHGTD.json",
        outputDirRoot=outputDir,
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFPerformanceConfig(ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=6),
        TrackSelectorRanges(pt=(1.0 * u.GeV, None), absEta=(None, 4.0), removeNeutral=True),
        outputDirRoot=outputDir,
    )


    return s


if "__main__" == __name__:

    # trackingGeometry=None

#    readAthenaMeasurements(trackingGeometry, inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.pool.root', outputDir=os.getcwd()).run()
    readAthenaMeasurements(trackingGeometry,
                           # inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.mu200_fullTruth.pool.root',
                           # inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.ttbar.pool.root',
                           # inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.minbias.pool.root',
                           # inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.ttbar.pool.root',
                           inputFile='/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/AOD.single_muon.pool.root',
                           outputDir=os.getcwd()).run()
