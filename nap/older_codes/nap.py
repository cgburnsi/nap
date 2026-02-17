import numpy as np

# Initialize constants and variables
title = ["" for _ in range(8)]
ui = np.zeros(21)
vi = np.zeros(21)
pi = np.zeros(21)
roi = np.zeros(21)

# Common blocks translated as dictionaries
AV = {"iav": 0, "cav": 0, "nst": 0, "smp": 0.0, "lss": 0, "cta": 0.0,
    "xmu": 0.0, "xla": 0.0, "rkmu": 0.0,
    "qut": np.zeros((81, 21)), "qvt": np.zeros((81, 21)), "qpt": np.zeros((81, 21))
}

onesid = {"ud": np.zeros(4), "vd": np.zeros(4), "pd": np.zeros(4), "rod": np.zeros(4)}

solutn = {
    "u": np.zeros((81, 21, 2)),
    "v": np.zeros((81, 21, 2)),
    "p": np.zeros((81, 21, 2)),
    "ro": np.zeros((81, 21, 2))
}

cntrlc = {
    "lmax": 0, "mmax": 0, "nmax": 0, "nprint": 0, "tconv": 0.0, "fdt": 0.0,
    "gamma": 0.0, "rgas": 0.0, "gam1": 0.0, "gam2": 0.0, "l1": 0, "l2": 0,
    "l3": 0, "m1": 0, "m2": 0, "dx": 0.0, "dy": 0.0, "dt": 0.0, "n": 0,
    "n1": 0, "n3": 0, "nasm": 0, "ivel": 0, "ichar": 0, "n1d": 0, "ljet": 0,
    "jflag": 0, "ierr": 0, "iui": 0, "iuo": 0, "dxr": 0.0, "dyr": 0.0,
    "ld": 0, "md": 0, "lmd1": 0, "lmd3": 0, "ib": 0, "rstar": 0.0,
    "rstars": 0.0, "nplot": 0, "g": 0.0, "pc": 0.0, "tc": 0.0, "lc": 0.0,
    "plow": 0.0, "rolow": 0.0
}

geometryc = {
    "nge": 0, "xi": 0.0, "ri": 0.0, "xt": 0.0, "rt": 0.0, "xe": 0.0,
    "re": 0.0, "rci": 0.0, "rct": 0.0, "angi": 0.0, "ange": 0.0,
    "xw": np.zeros(81), "yw": np.zeros(81), "xwi": np.zeros(81),
    "ywi": np.zeros(81), "nxny": np.zeros(81), "nwpts": 0, "iint": 0,
    "idif": 0, "lt": 0, "ndim": 0
}

bcc = {
    "pt": np.zeros(21), "tt": np.zeros(21), "theta": np.zeros(21),
    "pe": 0.0, "masse": 0.0, "massi": 0.0, "masst": 0.0, "thrust": 0.0,
    "nstag": 0
}

# Define default values
tconv = 0.0
fdt = 1.0
tstop = 1.0
nasm = 1
nstag = 0
name = 0
ipunch = 0
ngcb = 0
iintcb = 1
idifcb = 1
nstart = 0
tstart = 0.0
iint = 1
idif = 1
nmax = 0
nprint = 0
gamma = 1.4
rgas = 53.35
n1d = 1
ndim = 1
theta = [0.0] * 21
pe = 14.7
nst = 0
n = 0
iex = 1
nconvi = 1
ierr = 0
jflag = 0
iui = 1
iuo = 1
smp = 0.95
isuper = 0
iav = 0
cav = 4.0
nplot = -1
g = 32.174
pc = 144.0
tc = 460.0
lc = 12.0
iunit = 0
lss = 2
cta = 0.5
xmu = 0.2
xla = 1.0
rkmu = 0.7
plow = 0.01
rolow = 0.0001
rstar = 0.0
rstars = 0.0
