from Base.Container import Container

class Kinfit4PiContainer(Container):
    def __init__(self, path, variables):
        branches = {
            "KF_4Pi_IsConverged":   variables["PipPimPipPimKinfitIsConverged"],
            "KF_4Pi_Chi2":          variables["PipPimPipPimKinfitChi2"],
            "KF_4Pi_TrackMomenta":  variables["PipPimPipPimKinfitTrackMomenta"],
            "KF_4Pi_TrackThetas":   variables["PipPimPipPimKinfitTrackThetas"],
            "KF_4Pi_TrackPhis":     variables["PipPimPipPimKinfitTrackPhis"],
            "KF_4Pi_TrackEnergies": variables["PipPimPipPimKinfitTrackEnergies"],
            "KF_4Pi_TrackIndices":  variables["PipPimPipPimKinfitTrackIndices"],
        }
        Container.__init__(self, path, 'read', branches)
