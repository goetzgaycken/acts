#!/usr/bin/env python3
import os
# from acts.examples.odd import getOpenDataDetector
#from acts.examples import (
#    TrackingGeometryJsonReader
#    )
from acts.examples import (
#    GenericDetector,
    # AlignedDetector,
    TrackingGeometryJsonReader,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
#    CsvTrackingGeometryWriter,
#    ObjTrackingGeometryWriter,
    SvgTrackingGeometryWriter,
#    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

import acts

from acts import MaterialMapJsonConverter


def runGeometry(
    trackingGeometry,
    outputDir,
    events=1,
    outputSvg=False,
    outputJson=True,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)


        if outputSvg:
            outputDirSvg = os.path.join(outputDir, "svg")
            if not os.path.exists(outputDirSvg):
                os.makedirs(outputDirSvg)
            svgWriterCfg = SvgTrackingGeometryWriter.Config(outputDir=outputDirSvg)
            svgWriter = SvgTrackingGeometryWriter(
                level=acts.logging.INFO, config=svgWriterCfg
            )
            svgWriter.write(context, trackingGeometry)

        if outputJson:
            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=os.path.join(outputDir, "geometry-map"),
                writeFormat=JsonFormat.Json,
            )
            jmw.write(trackingGeometry)


if "__main__" == __name__:
    jsonTGReaderCfg = TrackingGeometryJsonReader.Config(detectorName="ITK",
                                                        toolLogLevel = acts.logging.VERBOSE,
                                                        logLevel = acts.logging.VERBOSE)
    jsonTGReader=TrackingGeometryJsonReader(jsonTGReaderCfg)
    trackingGeometry = jsonTGReader.read("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps-volbounds_2023.json")
    # detector, trackingGeometry, decorators = AlignedDetector.create()
    # detector, trackingGeometry, decorators = GenericDetector.create()
    # detector, trackingGeometry, decorators = getOpenDataDetector(getOpenDataDetectorDirectory() )

    runGeometry(trackingGeometry, outputDir=os.getcwd())
