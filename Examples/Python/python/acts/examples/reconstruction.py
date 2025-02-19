from pathlib import Path
from typing import Optional, Union
from enum import Enum
from collections import namedtuple
import warnings

import acts
import acts.examples

u = acts.UnitConstants

SeedingAlgorithm = Enum(
    "SeedingAlgorithm", "Default TruthSmeared TruthEstimated Orthogonal"
)

TruthSeedRanges = namedtuple(
    "TruthSeedRanges",
    ["rho", "z", "phi", "eta", "absEta", "pt", "nHits"],
    defaults=[(None, None)] * 7,
)

ParticleSmearingSigmas = namedtuple(
    "ParticleSmearingSigmas",
    ["d0", "d0PtA", "d0PtB", "z0", "z0PtA", "z0PtB", "t0", "phi", "theta", "pRel"],
    defaults=[None] * 10,
)

SeedFinderConfigArg = namedtuple(
    "SeedFinderConfig",
    [
        "maxSeedsPerSpM",
        "cotThetaMax",
        "sigmaScattering",
        "radLengthPerSeed",
        "minPt",
        "impactMax",
        "interactionPointCut",
        "arithmeticAverageCotTheta",
        "deltaZMax",
        "maxPtScattering",
        "zBinEdges",
        "skipPreviousTopSP",
        "zBinsCustomLooping",
        "rRangeMiddleSP",
        "useVariableMiddleSPRange",
        "binSizeR",
        "forceRadialSorting",
        "seedConfirmation",
        "centralSeedConfirmationRange",
        "forwardSeedConfirmationRange",
        "deltaR",  # (min,max)
        "deltaRBottomSP",  # (min,max)
        "deltaRTopSP",  # (min,max)
        "deltaRMiddleSPRange",  # (min,max)
        "collisionRegion",  # (min,max)
        "r",  # (min,max)
        "z",  # (min,max)
    ],
    defaults=[None] * 20 + [(None, None)] * 7,
)
SeedFinderOptionsArg = namedtuple(
    "SeedFinderOptions", ["beamPos", "bFieldInZ"], defaults=[(None, None), None]
)

SeedFilterConfigArg = namedtuple(
    "SeedFilterConfig",
    [
        "impactWeightFactor",
        "zOriginWeightFactor",
        "compatSeedWeight",
        "compatSeedLimit",
        "numSeedIncrement",
        "seedWeightIncrement",
        "seedConfirmation",
        "curvatureSortingInFilter",
        "maxSeedsPerSpMConf",
        "maxQualitySeedsPerSpMConf",
        "useDeltaRorTopRadius",
        "deltaRMin",
    ],
    defaults=[None] * 12,
)

SpacePointGridConfigArg = namedtuple(
    "SeedGridConfig",
    [
        "rMax",
        "zBinEdges",
        "phiBinDeflectionCoverage",
        "impactMax",
        "deltaRMax",
        "phi",  # (min,max)
    ],
    defaults=[None] * 5 + [(None, None)] * 1,
)

SeedingAlgorithmConfigArg = namedtuple(
    "SeedingAlgorithmConfig",
    [
        "allowSeparateRMax",
        "zBinNeighborsTop",
        "zBinNeighborsBottom",
        "numPhiNeighbors",
    ],
    defaults=[None] * 4,
)

TrackParamsEstimationConfig = namedtuple(
    "TrackParamsEstimationConfig",
    [
        "deltaR",  # (min,max)
    ],
    defaults=[(None, None)],
)

TrackSelectorRanges = namedtuple(
    "TrackSelectorRanges",
    [
        "loc0",
        "loc1",
        "time",
        "eta",
        "absEta",
        "pt",
        "phi",
        "removeNeutral",
        "removeCharged",
    ],
    defaults=[(None, None)] * 7 + [None] * 2,
)

CKFPerformanceConfig = namedtuple(
    "CKFPerformanceConfig",
    ["truthMatchProbMin", "nMeasurementsMin", "ptMin"],
    defaults=[None] * 3,
)

AmbiguityResolutionConfig = namedtuple(
    "AmbiguityResolutionConfig",
    ["maximumSharedHits"],
    defaults=[None] * 1,
)


class VertexFinder(Enum):
    Truth = (1,)
    AMVF = (2,)
    Iterative = (3,)


@acts.examples.NamedTypeArgs(
    seedingAlgorithm=SeedingAlgorithm,
    truthSeedRanges=TruthSeedRanges,
    particleSmearingSigmas=ParticleSmearingSigmas,
    seedFinderConfigArg=SeedFinderConfigArg,
    seedFinderOptionsArg=SeedFinderOptionsArg,
    seedFilterConfigArg=SeedFilterConfigArg,
    spacePointGridConfigArg=SpacePointGridConfigArg,
    seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg,
    trackParamsEstimationConfig=TrackParamsEstimationConfig,
    logLevel=acts.logging.Level,
)
def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.Default,
    truthSeedRanges: Optional[TruthSeedRanges] = TruthSeedRanges(),
    particleSmearingSigmas: ParticleSmearingSigmas = ParticleSmearingSigmas(),
    initialVarInflation: Optional[list] = None,
    seedFinderConfigArg: SeedFinderConfigArg = SeedFinderConfigArg(),
    seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(),
    seedFilterConfigArg: SeedFilterConfigArg = SeedFilterConfigArg(),
    spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
    trackParamsEstimationConfig: TrackParamsEstimationConfig = TrackParamsEstimationConfig(),
    inputParticles: str = "particles_initial",
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
) -> None:
    """This function steers the seeding
    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    truthSeedRanges : TruthSeedRanges(rho, z, phi, eta, absEta, pt, nHits)
        TruthSeedSelector configuration. Each range is specified as a tuple of (min,max).
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TruthSeedSelector.hpp
        If specified as None, don't run ParticleSmearing at all (and use addCKFTracks(selectedParticles="particles_initial"))
    particleSmearingSigmas : ParticleSmearingSigmas(d0, d0PtA, d0PtB, z0, z0PtA, z0PtB, t0, phi, theta, pRel)
        ParticleSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSmearing.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSmearing.hpp
    seedFinderConfigArg : SeedFinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, impactMax, interactionPointCut, arithmeticAverageCotTheta, deltaZMax, maxPtScattering, zBinEdges, skipPreviousTopSP, zBinsCustomLooping, rRangeMiddleSP, useVariableMiddleSPRange, binSizeR, forceRadialSorting, seedConfirmation, centralSeedConfirmationRange, forwardSeedConfirmationRange, deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z)
        SeedFinderConfig settings. deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z are ranges specified as a tuple of (min,max). beamPos is specified as (x,y).
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFinderOptionsArg :  SeedFinderOptionsArg(bFieldInZ, beamPos)
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFilterConfigArg : SeedFilterConfigArg(compatSeedWeight, compatSeedLimit, numSeedIncrement, seedWeightIncrement, seedConfirmation, curvatureSortingInFilter, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf, useDeltaRorTopRadius)
                                Defaults specified in Core/include/Acts/Seeding/SeedFilterConfig.hpp
    spacePointGridConfigArg : SpacePointGridConfigArg(rMax, zBinEdges, phiBinDeflectionCoverage, phi, impactMax)
                                SpacePointGridConfigArg settings. phi is specified as a tuple of (min,max).
        Defaults specified in Core/include/Acts/Seeding/SpacePointGrid.hpp
    seedingAlgorithmConfigArg : SeedingAlgorithmConfigArg(allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors)
                                Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/SeedingAlgorithm.hpp
    trackParamsEstimationConfig : TrackParamsEstimationConfig(deltaR)
        TrackParamsEstimationAlgorithm configuration. Currently only deltaR=(min,max) range specified here.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    inputParticles : str, "particles_initial"
        input particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    logLevel = acts.examples.defaultLogging(s, logLevel)()
    logger = acts.logging.getLogger("addSeeding")

    if truthSeedRanges is not None:
        selectedParticles = addSeedingTruthSelection(
            s, inputParticles, truthSeedRanges, logLevel
        )
    else:
        selectedParticles = inputParticles

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        logger.info("Using smeared truth particles for seeding")
        addTruthSmearing(
            s,
            rnd,
            selectedParticles,
            particleSmearingSigmas,
            initialVarInflation,
            logLevel,
        )
    else:
        spacePoints = addSpacePointsMaking(
            s, trackingGeometry, geoSelectionConfigFile, acts.logging.VERBOSE # logLevel
        )
        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            logger.info("Using truth track finding from space points for seeding")
            # Use truth tracking
            truthTrackFinder = acts.examples.TruthTrackFinder(
                level=logLevel,
                inputParticles=selectedParticles,
                inputMeasurementParticlesMap="measurement_particles_map",
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(truthTrackFinder)
            inputProtoTracks = truthTrackFinder.config.outputProtoTracks
            inputSeeds = ""
        elif seedingAlgorithm == SeedingAlgorithm.Default:
            logger.info("Using default seeding")
            inputProtoTracks, inputSeeds = addStandardSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                acts.logging.VERBOSE # logLevel,
                , trackingGeometry=trackingGeometry
            )
        elif seedingAlgorithm == SeedingAlgorithm.Orthogonal:
            logger.info("Using orthogonal seeding")
            inputProtoTracks, inputSeeds = addOrthogonalSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                logLevel,
            )
        else:
            logger.fatal("unknown seedingAlgorithm %s", seedingAlgorithm)

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=logLevel,
            inputSeeds=inputSeeds,
            inputProtoTracks=inputProtoTracks,
            inputSpacePoints=[spacePoints],
            inputSourceLinks="sourcelinks",
            outputTrackParameters="estimatedparameters",
            outputProtoTracks="prototracks_estimated",
             trackingGeometry=trackingGeometry,
            magneticField=field,
            **acts.examples.defaultKWArgs(
                initialVarInflation=initialVarInflation,
                deltaRMin=trackParamsEstimationConfig.deltaR[0],
                deltaRMax=trackParamsEstimationConfig.deltaR[1],
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        if outputDirRoot is not None:
            addSeedPerformanceWriters(
                s,
                outputDirRoot,
                inputProtoTracks,
                selectedParticles,
                inputParticles,
                parEstimateAlg.config.outputTrackParameters,
                parEstimateAlg.config.outputProtoTracks,
                logLevel,
            )
    return s


def addSeedingTruthSelection(
    s: acts.examples.Sequencer,
    inputParticles: str,
    truthSeedRanges: str,
    logLevel: acts.logging.Level = None,
):
    """adds truth particles filtering before filtering
    For parameters description see addSeeding
    """
    selAlg = acts.examples.TruthSeedSelector(
        **acts.examples.defaultKWArgs(
            ptMin=truthSeedRanges.pt[0],
            ptMax=truthSeedRanges.pt[1],
            etaMin=truthSeedRanges.eta[0],
            etaMax=truthSeedRanges.eta[1],
            nHitsMin=truthSeedRanges.nHits[0],
            nHitsMax=truthSeedRanges.nHits[1],
            rhoMin=truthSeedRanges.rho[0],
            rhoMax=truthSeedRanges.rho[1],
            zMin=truthSeedRanges.z[0],
            zMax=truthSeedRanges.z[1],
            phiMin=truthSeedRanges.phi[0],
            phiMax=truthSeedRanges.phi[1],
            absEtaMin=truthSeedRanges.absEta[0],
            absEtaMax=truthSeedRanges.absEta[1],
        ),
        level=logLevel,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap="measurement_particles_map",
        outputParticles="truth_seeds_selected",
    )
    s.addAlgorithm(selAlg)
    return selAlg.config.outputParticles


def addTruthSmearing(
    sequence: acts.examples.Sequencer,
    rnd: acts.examples.RandomNumbers,
    selectedParticles: str,
    particleSmearingSigmas: ParticleSmearingSigmas,
    initialVarInflation: float,
    logLevel: acts.logging.Level = None,
):
    """adds algorithm that would mimic detector response uncertainties for truth seeding
    For parameters description see addSeeding
    """

    rnd = rnd or acts.examples.RandomNumbers(seed=42)
    # Run particle smearing
    ptclSmear = acts.examples.ParticleSmearing(
        level=logLevel,
        inputParticles=selectedParticles,
        outputTrackParameters="estimatedparameters",
        randomNumbers=rnd,
        # gaussian sigmas to smear particle parameters
        **acts.examples.defaultKWArgs(
            sigmaD0=particleSmearingSigmas.d0,
            sigmaD0PtA=particleSmearingSigmas.d0PtA,
            sigmaD0PtB=particleSmearingSigmas.d0PtB,
            sigmaZ0=particleSmearingSigmas.z0,
            sigmaZ0PtA=particleSmearingSigmas.z0PtA,
            sigmaZ0PtB=particleSmearingSigmas.z0PtB,
            sigmaT0=particleSmearingSigmas.t0,
            sigmaPhi=particleSmearingSigmas.phi,
            sigmaTheta=particleSmearingSigmas.theta,
            sigmaPRel=particleSmearingSigmas.pRel,
            initialVarInflation=initialVarInflation,
        ),
    )
    sequence.addAlgorithm(ptclSmear)


def addSpacePointsMaking(
    sequence: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geoSelectionConfigFile: Union[Path, str],
    logLevel: acts.logging.Level = None,
):
    """adds space points making
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputSourceLinks="sourcelinks",
        inputMeasurements="measurements",
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
    )
    sequence.addAlgorithm(spAlg)
    return spAlg.config.outputSpacePoints


def addStandardSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFilterConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
    trackingGeometry = None,
):
    """adds standard seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    seedFinderConfig = acts.SeedFinderConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            arithmeticAverageCotTheta=seedFinderConfigArg.arithmeticAverageCotTheta,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            skipPreviousTopSP=seedFinderConfigArg.skipPreviousTopSP,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            forceRadialSorting=seedFinderConfigArg.forceRadialSorting,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=acts.Vector2(0.0, 0.0)
            if seedFinderOptionsArg.beamPos == (None, None)
            else acts.Vector2(
                seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfig.deltaRMin
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
            curvatureSortingInFilter=seedFilterConfigArg.curvatureSortingInFilter,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    gridConfig = acts.SpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfig.minPt,
            rMax=(
                seedFinderConfig.rMax
                if spacePointGridConfigArg.rMax is None
                else spacePointGridConfigArg.rMax
            ),
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            deltaRMax=(
                seedFinderConfig.deltaRMax
                if spacePointGridConfigArg.deltaRMax is None
                else spacePointGridConfigArg.deltaRMax
            ),
            cotThetaMax=seedFinderConfig.cotThetaMax,
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            impactMax=spacePointGridConfigArg.impactMax,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
        )
    )

    gridOptions = acts.SpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            numPhiNeighbors=seedingAlgorithmConfigArg.numPhiNeighbors,
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        trackingGeometry=trackingGeometry,
        stop=True,
    )
    sequence.addAlgorithm(seedingAlg)
    return seedingAlg.config.outputProtoTracks, seedingAlg.config.outputSeeds


def addOrthogonalSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds orthogonal seeding algorithm
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    seedFinderConfig = acts.SeedFinderOrthogonalConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=acts.Vector2(0.0, 0.0)
            if seedFinderOptionsArg.beamPos == (None, None)
            else acts.Vector2(
                seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfigArg.deltaR[0]
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            curvatureSortingInFilter=seedFilterConfigArg.curvatureSortingInFilter,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )
    seedingAlg = acts.examples.SeedingOrthogonalAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
    )
    sequence.addAlgorithm(seedingAlg)
    return seedingAlg.config.outputProtoTracks, seedingAlg.config.outputSeeds


def addSeedPerformanceWriters(
    sequence: acts.examples.Sequencer,
    outputDirRoot: Union[Path, str],
    inputProtoTracks: str,
    selectedParticles: str,
    inputParticles: str,
    outputTrackParameters: str,
    outputProtoTracks: str,
    logLevel: acts.logging.Level = None,
):
    """Writes seeding related performance output"""
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir()
    sequence.addWriter(
        acts.examples.TrackFinderPerformanceWriter(
            level=customLogLevel(),
            inputProtoTracks=inputProtoTracks,
            inputParticles=selectedParticles,  # the original selected particles after digitization
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDirRoot / "performance_seeding_trees.root"),
        )
    )

    sequence.addWriter(
        acts.examples.SeedingPerformanceWriter(
            level=customLogLevel(minLevel=acts.logging.DEBUG),
            inputProtoTracks=inputProtoTracks,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDirRoot / "performance_seeding_hists.root"),
        )
    )

    if False :
        sequence.addWriter(
            acts.examples.RootTrackParameterWriter(
                level=customLogLevel(),
                inputTrackParameters=outputTrackParameters,
                inputProtoTracks=outputProtoTracks,
                inputParticles=inputParticles,
                inputSimHits="simhits",
                inputMeasurementParticlesMap="measurement_particles_map",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / "estimatedparams.root"),
                treeName="estimatedparams",
            )
        )


def addKalmanTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    directNavigation=False,
    reverseFilteringMomThreshold=0 * u.GeV,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=customLogLevel(),
        inputParticles="truth_seeds_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputProtoTracks="prototracks",
    )
    s.addAlgorithm(truthTrkFndAlg)

    if directNavigation:
        srfSortAlg = acts.examples.SurfaceSortingAlgorithm(
            level=customLogLevel(),
            inputProtoTracks="prototracks",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            outputProtoTracks="sortedprototracks",
        )
        s.addAlgorithm(srfSortAlg)
        inputProtoTracks = srfSortAlg.config.outputProtoTracks
    else:
        inputProtoTracks = "prototracks"

    kalmanOptions = {
        "multipleScattering": True,
        "energyLoss": True,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        outputTrajectories="kfTrajectories",
        directNavigation=directNavigation,
        pickTrack=-1,
        trackingGeometry=trackingGeometry,
        fit=acts.examples.makeKalmanFitterFunction(
            trackingGeometry, field, **kalmanOptions
        ),
    )
    s.addAlgorithm(fitAlg)

    s.addWhiteboardAlias("trajectories", fitAlg.config.outputTrajectories)

    return s


def addTruthTrackingGsf(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    gsfOptions = {
        "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
        "maxComponents": 12,
        "abortOnError": False,
        "disableAllMaterialHandling": False,
        "finalReductionMethod": acts.examples.FinalReductionMethod.maxWeight,
        "weightCutoff": 1.0e-4,
        "level": customLogLevel(),
    }

    gsfAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputProtoTracks="prototracks",
        inputInitialTrackParameters="estimatedparameters",
        outputTrajectories="gsf_trajectories",
        directNavigation=False,
        pickTrack=-1,
        trackingGeometry=trackingGeometry,
        fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions),
    )

    s.addAlgorithm(gsfAlg)

    return s


@acts.examples.NamedTypeArgs(
    ckfPerformanceConfig=CKFPerformanceConfig,
    trackSelectorRanges=TrackSelectorRanges,
)
def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorRanges: Optional[TrackSelectorRanges] = None,
    selectedParticles: str = "truth_seeds_selected",
    ckfPerformanceConfig: CKFPerformanceConfig = CKFPerformanceConfig(),
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    ckfPerformanceConfigArg : CKFPerformanceConfig(truthMatchProbMin, nMeasurementsMin, ptMin)
        CKFPerformanceWriter configuration.
        Defaults specified in Examples/Io/Performance/ActsExamples/Io/Performance/CKFPerformanceWriter.hpp
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    selectedParticles : str, "truth_seeds_selected"
        CKFPerformanceWriter truth input
    trackSelectorRanges : TrackSelectorRanges(loc0, loc1, time, eta, absEta, pt, phi, removeNeutral, removeCharged)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrajectories : bool, True
        write trackstates_ckf.root and tracksummary_ckf.root ntuples? These can be quite large.
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [(acts.GeometryIdentifier(), ([], [15.0], [10]))]
        ),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputInitialTrackParameters="estimatedparameters",
        outputTrajectories="ckfTrajectories",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field
        ),
    )
    s.addAlgorithm(trackFinder)

    s.addWhiteboardAlias("trajectories", trackFinder.config.outputTrajectories)

    if trackSelectorRanges is not None:
        trackSelector = addTrackSelection(
            s,
            trackSelectorRanges,
            inputTrajectories=trackFinder.config.outputTrajectories,
            outputTrajectories="selectedTrajectories",
            logLevel=customLogLevel(),
        )

        s.addWhiteboardAlias("trajectories", trackSelector.config.outputTrajectories)

    addTrajectoryWriters(
        s,
        name="ckf",
        trajectories=trackFinder.config.outputTrajectories,
        ckfPerformanceConfig=ckfPerformanceConfig,
        trackingGeometry=trackingGeometry,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        writeFinderPerformance=False,
        writeFitterPerformance=False,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    ckfPerformanceConfig=CKFPerformanceConfig,
    trackSelectorRanges=TrackSelectorRanges,
)
def addTrajectoryWriters(
    s: acts.examples.Sequencer,
    name: str,
    trajectories: str = "trajectories",
    ckfPerformanceConfig: CKFPerformanceConfig = CKFPerformanceConfig(),
    trackingGeometry: acts.TrackingGeometry = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeStates: bool = False,
    writeSummary: bool = True,
    writeCKFperformance: bool = True,
    writeFinderPerformance: bool = True,
    writeFitterPerformance: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        if writeStates and False:
            # write track states from CKF
            trackStatesWriter = acts.examples.RootTrajectoryStatesWriter(
                level=customLogLevel(),
                inputTrajectories=trajectories,
                # @note The full particles collection is used here to avoid lots of warnings
                # since the unselected CKF track might have a majority particle not in the
                # filtered particle collection. This could be avoided when a seperate track
                # selection algorithm is used.
                inputParticles="particles_selected",
                inputSimHits="simhits",
                inputMeasurementParticlesMap="measurement_particles_map",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeSummary:
            # write track summary from CKF
            trackSummaryWriter = acts.examples.RootTrajectorySummaryWriter(
                level=customLogLevel(),
                inputTrajectories=trajectories,
                # @note The full particles collection is used here to avoid lots of warnings
                # since the unselected CKF track might have a majority particle not in the
                # filtered particle collection. This could be avoided when a seperate track
                # selection algorithm is used.
                inputParticles="particles_final", #"particles_selected",
                inputMeasurementParticlesMap="measurement_particles_map",
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
            )
            s.addWriter(trackSummaryWriter)

        if writeCKFperformance:
            # Write CKF performance data
            ckfPerfWriter = acts.examples.CKFPerformanceWriter(
                level=customLogLevel(),
                inputParticles="truth_seeds_selected",
                inputTrajectories=trajectories,
                inputMeasurementParticlesMap="measurement_particles_map",
                trackingGeometry=trackingGeometry,
                dumpDuplicates=True,
                **acts.examples.defaultKWArgs(
                    # The bottom seed could be the first, second or third hits on the truth track
                    nMeasurementsMin=ckfPerformanceConfig.nMeasurementsMin,
                    ptMin=ckfPerformanceConfig.ptMin,
                    truthMatchProbMin=ckfPerformanceConfig.truthMatchProbMin,
                ),
                filePath=str(outputDirRoot / f"performance_{name}.root"),
            )
            s.addWriter(ckfPerfWriter)

        if writeFinderPerformance:
            s.addWriter(
                acts.examples.TrackFinderPerformanceWriter(
                    level=acts.logging.INFO,
                    inputProtoTracks="prototracks",
                    inputParticles="truth_seeds_selected",
                    inputMeasurementParticlesMap="measurement_particles_map",
                    filePath=str(
                        outputDirRoot / f"performance_track_finder_{name}.root"
                    ),
                )
            )

        if writeFitterPerformance:
            s.addWriter(
                acts.examples.TrackFitterPerformanceWriter(
                    level=acts.logging.INFO,
                    inputTrajectories="trajectories",
                    inputParticles="truth_seeds_selected",
                    inputMeasurementParticlesMap="measurement_particles_map",
                    filePath=str(
                        outputDirRoot / f"performance_track_fitter_{name}.root"
                    ),
                )
            )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        if writeSummary:
            csvMTJWriter = acts.examples.CsvMultiTrajectoryWriter(
                level=customLogLevel(),
                inputTrajectories=trajectories,
                inputMeasurementParticlesMap="measurement_particles_map",
                outputDir=str(outputDirCsv),
            )
            s.addWriter(csvMTJWriter)


@acts.examples.NamedTypeArgs(
    trackSelectorRanges=TrackSelectorRanges,
)
def addTrackSelection(
    s: acts.examples.Sequencer,
    trackSelectorRanges: TrackSelectorRanges,
    inputTrackParameters: Optional[str] = None,
    inputTrajectories: Optional[str] = None,
    outputTrackParameters: Optional[str] = None,
    outputTrajectories: Optional[str] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.TrackSelector:

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    trackSelector = acts.examples.TrackSelector(
        level=customLogLevel(),
        inputTrackParameters=inputTrackParameters
        if inputTrackParameters is not None
        else "",
        inputTrajectories=inputTrajectories if inputTrajectories is not None else "",
        outputTrackParameters=outputTrackParameters
        if outputTrackParameters is not None
        else "",
        outputTrajectories=outputTrajectories if outputTrajectories is not None else "",
        **acts.examples.defaultKWArgs(
            loc0Min=trackSelectorRanges.loc0[0],
            loc0Max=trackSelectorRanges.loc0[1],
            loc1Min=trackSelectorRanges.loc1[0],
            loc1Max=trackSelectorRanges.loc1[1],
            timeMin=trackSelectorRanges.time[0],
            timeMax=trackSelectorRanges.time[1],
            phiMin=trackSelectorRanges.phi[0],
            phiMax=trackSelectorRanges.phi[1],
            etaMin=trackSelectorRanges.eta[0],
            etaMax=trackSelectorRanges.eta[1],
            absEtaMin=trackSelectorRanges.absEta[0],
            absEtaMax=trackSelectorRanges.absEta[1],
            ptMin=trackSelectorRanges.pt[0],
            ptMax=trackSelectorRanges.pt[1],
            removeCharged=trackSelectorRanges.removeCharged,
            removeNeutral=trackSelectorRanges.removeNeutral,
        ),
    )

    s.addAlgorithm(trackSelector)

    return trackSelector


ExaTrkXBackend = Enum("ExaTrkXBackend", "Torch Onnx")


def addExaTrkX(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geometrySelection: Union[Path, str],
    modelDir: Union[Path, str],
    outputDirRoot: Optional[Union[Path, str]] = None,
    backend: Optional[ExaTrkXBackend] = ExaTrkXBackend.Torch,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
    s.addAlgorithm(
        acts.examples.TruthSeedSelector(
            level=customLogLevel(),
            ptMin=500 * u.MeV,
            nHitsMin=9,
            inputParticles="particles_initial",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputParticles="particles_seed_selected",
        )
    )

    # Create space points
    s.addAlgorithm(
        acts.examples.SpacePointMaker(
            level=customLogLevel(),
            inputSourceLinks="sourcelinks",
            inputMeasurements="measurements",
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
    )

    # For now we don't configure only the common options so this works
    exaTrkxModule = (
        acts.examples.ExaTrkXTrackFindingTorch
        if backend == ExaTrkXBackend.Torch
        else acts.examples.ExaTrkXTrackFindingOnnx
    )

    exaTrkxFinding = exaTrkxModule(
        modelDir=str(modelDir),
        spacepointFeatures=3,
        embeddingDim=8,
        rVal=1.6,
        knnVal=500,
        filterCut=0.21,
    )

    s.addAlgorithm(
        acts.examples.TrackFindingAlgorithmExaTrkX(
            level=customLogLevel(),
            inputSpacePoints="spacepoints",
            outputProtoTracks="protoTracks",
            trackFinderML=exaTrkxFinding,
        )
    )

    # Write truth track finding / seeding performance
    if outputDirRoot is not None:
        s.addWriter(
            acts.examples.TrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputProtoTracks="protoTracks",
                inputParticles="particles_initial",  # the original selected particles after digitization
                inputMeasurementParticlesMap="measurement_particles_map",
                filePath=str(Path(outputDirRoot) / "performance_seeding_trees.root"),
            )
        )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionConfig,
    ckfPerformanceConfig=CKFPerformanceConfig,
)
def addAmbiguityResolution(
    s,
    config: AmbiguityResolutionConfig = AmbiguityResolutionConfig(),
    selectedParticles: str = "truth_seeds_selected",
    ckfPerformanceConfig: CKFPerformanceConfig = CKFPerformanceConfig(),
    trackingGeometry: acts.TrackingGeometry = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:

    from acts.examples import AmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    alg = AmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputSourceLinks="sourcelinks",
        inputTrajectories="trajectories",
        outputTrajectories="filteredTrajectories",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
        ),
    )
    s.addAlgorithm(alg)

    s.addWhiteboardAlias("trajectories", alg.config.outputTrajectories)

    addTrajectoryWriters(
        s,
        name="ambi",
        trajectories=alg.config.outputTrajectories,
        ckfPerformanceConfig=ckfPerformanceConfig,
        trackingGeometry=trackingGeometry,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        writeFinderPerformance=False,
        writeFitterPerformance=False,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    trackSelectorRanges=TrackSelectorRanges,
)
def addVertexFitting(
    s,
    field,
    outputDirRoot: Optional[Union[Path, str]] = None,
    trajectories: Optional[str] = "trajectories",
    trackParameters: Optional[str] = None,
    associatedParticles: Optional[str] = None,
    vertexFinder: VertexFinder = VertexFinder.Truth,
    trackSelectorRanges: Optional[TrackSelectorRanges] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:

    """This function steers the vertex fitting

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addVertexFitting)
    field : magnetic field
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    associatedParticles : str, "associatedTruthParticles"
        RootVertexPerformanceWriter.inputAssociatedTruthParticles
    vertexFinder : VertexFinder, Truth
        vertexFinder algorithm: one of Truth, AMVF, Iterative
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    """
    from acts.examples import (
        TruthVertexFinder,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        AdaptiveMultiVertexFinderAlgorithm,
        RootVertexPerformanceWriter,
    )
    print ('DEBUG args',type(field),field)
    print ('DEBUG args',type(outputDirRoot),outputDirRoot)
    print ('DEBUG args',type(associatedParticles),associatedParticles)
    print ('DEBUG args',type(trajectories),trajectories)
    print ('DEBUG args',type(trackParameters),trackParameters)
    print ('DEBUG args',type(trackParametersTips),trackParametersTips)
    print ('DEBUG args',type(vertexFinder),vertexFinder)
    print ('DEBUG args',type(logLevel),logLevel)

    trajectories = trajectories if trajectories is not None else ""
    trackParameters = trackParameters if trackParameters is not None else ""
    associatedParticles = associatedParticles if associatedParticles is not None else ""

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if trackSelectorRanges is not None:
        trackSelector = addTrackSelection(
            s,
            trackSelectorRanges,
            inputTrackParameters=trackParameters,
            inputTrajectories=trajectories,
            outputTrackParameters="selectedTrackParametersVertexing",
            outputTrajectories="selectedTrajectoriesVertexing",
            logLevel=customLogLevel(),
        )

        trajectories = trackSelector.config.outputTrajectories if trajectories else ""
        trackParameters = (
            trackSelector.config.outputTrackParameters if trackParameters else ""
        )

    # inputParticles = "particles_input"
    #selectedParticles = "particles_selected"
    inputParticles = "particles_final"
    selectedParticles = "particles_final"
    outputVertices = "fittedVertices"

    outputTime = ""
    if vertexFinder == VertexFinder.Truth:
        findVertices = TruthVertexFinder(
            level=customLogLevel(),
            inputParticles=selectedParticles,
            outputProtoVertices="protovertices",
            excludeSecondaries=True,
        )
        s.addAlgorithm(findVertices)
        fitVertices = VertexFitterAlgorithm(
            level=customLogLevel(),
            bField=field,
            inputTrajectories=trajectories,
            inputTrackParameters=trackParameters,
            inputProtoVertices=findVertices.config.outputProtoVertices,
            outputVertices=outputVertices,
        )
        s.addAlgorithm(fitVertices)
    elif vertexFinder == VertexFinder.Iterative:
        findVertices = IterativeVertexFinderAlgorithm(
            level=customLogLevel(),
            bField=field,
            inputTrajectories=trajectories,
            inputTrackParameters=trackParameters,
            outputProtoVertices="protovertices",
            outputVertices=outputVertices,
        )
        s.addAlgorithm(findVertices)
    elif vertexFinder == VertexFinder.AMVF:
        outputTime = "outputTime"
        findVertices = AdaptiveMultiVertexFinderAlgorithm(
            level=customLogLevel(),
            bField=field,
            inputTrajectories=trajectories,
            inputTrackParameters=trackParameters,
            outputProtoVertices="protovertices",
            outputVertices=outputVertices,
            outputTime=outputTime,
        )
        s.addAlgorithm(findVertices)
    else:
        raise RuntimeError("Invalid finder argument")

    print('DEBUG output dir',outputDirRoot)
    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        if associatedParticles == selectedParticles:
            warnings.warn(
                "Using RootVertexPerformanceWriter with smeared particles is not necessarily supported. "
                "Please get in touch with us"
            )
        s.addWriter(
            RootVertexPerformanceWriter(
                level=customLogLevel(),
                inputAllTruthParticles=inputParticles,
                inputSelectedTruthParticles=selectedParticles,
                inputMeasurementParticlesMap="measurement_particles_map",
                inputTrajectories=trajectories,
                inputTrackParameters=trackParameters,
                inputAssociatedTruthParticles=associatedParticles,
                inputVertices=outputVertices,
                minTrackVtxMatchFraction=0.5 if associatedParticles else 0.0,
                inputTime=outputTime,
                treeName="vertexing",
                filePath=str(outputDirRoot / "performance_vertexing.root"),
            )
        )

    return s
