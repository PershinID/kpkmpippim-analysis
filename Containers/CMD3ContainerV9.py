from Base.Container import Container

class CMD3ContainerV9(Container):
    def __init__(self, path: str, variables):
        branches = [
            "nt",
            "ntlxe",
            "nph",
            
            "emeas",
            "demeas",
            "xbeam",
            "ybeam",
            "runnum",
            "evnum",
            "ecaltot",
            "ecalneu",
            "psumch",
            "psumnu",
            "nt_total",
            "nv_total",
            "tnhit",
            "tlength",
            "tphi",
            "tth",
            "tptot",
            "tphiv",
            "tthv",
            "tptotv",
            "trho",
            "tdedx",
            "tz",
            "tchi2r",
            "tchi2z",
            "tchi2ndf",
            "tt0",
            "tant",
            "tcharge",
            "ten",
            "tfc",
            "tenlxe",
            "tlengthlxe",
            "tenslxe_layers",
            "tencsi",
            "tenbgo",
            "tclth",
            "tclphi",
            "terr",
            "terr0",
            "tindlxe",
            "txyzatcl",
            "txyzatlxe",
            "tenconv",
            "ntlxe_total",
            "ntlxelayers",
            "tlxenhit",
            "tlxelength",
            "tlxededx",
            "tlxeir",
            "tlxeitheta",
            "tlxeiphi",
            "tlxevtheta",
            "tlxevphi",
            "tlxechi2",
            "tlxesen",
            "tlxesen_layers",
            "finalstate_id",
            "nph_total",
            "phen",
            "phth",
            "phphi",
            "phrho",
            "phen0",
            "phth0",
            "phphi0",
            "phlxe",
            "phslxe_layers",
            "pherr",
            "phcsi",
            "phbgo",
            "phflag",
            "phconv",
            "phfc",
        ]
        Container.__init__(self, path, "read", {branch: variables[branch] for branch in branches})
