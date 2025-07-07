from Base.Container import Container

class Kinfit2K2PiContainer(Container):
    def __init__(self, path, variables):
        branches = {
            "KF_2K2Pi_IsConverged":   variables["KpKmPipPimKinfitIsConverged"],
            "KF_2K2Pi_Chi2":          variables["KpKmPipPimKinfitChi2"],
            "KF_2K2Pi_TrackMomenta":  variables["KpKmPipPimKinfitTrackMomenta"],
            "KF_2K2Pi_TrackThetas":   variables["KpKmPipPimKinfitTrackThetas"],
            "KF_2K2Pi_TrackPhis":     variables["KpKmPipPimKinfitTrackPhis"],
            "KF_2K2Pi_TrackEnergies": variables["KpKmPipPimKinfitTrackEnergies"],
            "KF_2K2Pi_TrackIndices":  variables["KpKmPipPimKinfitTrackIndices"],
        }
        Container.__init__(self, path, 'read', branches)
