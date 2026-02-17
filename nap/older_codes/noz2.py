import numpy as np

class GlobalData:
    """Class to encapsulate global data."""
    def __init__(self):
        # Variables in AV common block
        self.IAV = 0
        self.NST = 0
        self.SMP = 0
        self.LSS = 0
        self.CAV = ""
        self.CTA = 0.0
        self.XMU = 0.0
        self.XLA = 0.0
        self.RKMU = 0.0
        self.QUT = np.zeros((81, 21))
        self.QVT = np.zeros((81, 21))
        self.QPT = np.zeros((81, 21))

        # Variables in ONESID common block
        self.UD = np.zeros(4)
        self.VD = np.zeros(4)
        self.PD = np.zeros(4)
        self.ROD = np.zeros(4)

        # Variables in SOLUTN common block
        self.U = np.zeros((81, 21, 2))
        self.V = np.zeros((81, 21, 2))
        self.P = np.zeros((81, 21, 2))
        self.RO = np.zeros((81, 21, 2))

        # Variables in CNTRLC common block
        self.LMAX = 0
        self.MMAX = 0
        self.NMAX = 0
        self.NPRINT = 0
        self.NASM = 0
        self.IVEL = 0
        self.ICHAR = 0
        self.N1D = 0
        self.LJET = 0
        self.JFLAG = 0
        self.IERR = 0
        self.IUI = 0
        self.IUO = 0
        self.LD = 0
        self.MD = 0
        self.LMD1 = 0
        self.LMD3 = 0
        self.IB = 0
        self.NPLOT = 0
        self.TCONV = 0.0
        self.FDT = 0.0
        self.GAMMA = 1.4
        self.RGAS = 287.05
        self.GAM1 = 0.0
        self.GAM2 = 0.0
        self.RSTAR = 0.0
        self.RSTARS = 0.0
        self.DX = 0.1
        self.DY = 0.1
        self.DT = 0.01
        self.G = 0.0
        self.PC = 0.0
        self.TC = 0.0
        self.PLOW = 0.0
        self.ROLOW = 0.0

        # Variables in GEMTRYC common block
        self.NGEOM = 0
        self.NXNY = np.zeros(81, dtype=int)
        self.NWPTS = 0
        self.IINT = 0
        self.IDIF = 0
        self.LT = 0
        self.NDIM = 0
        self.XI = 0.0
        self.RI = 0.0
        self.XT = 0.0
        self.RT = 0.0
        self.XE = 0.0
        self.RE = 0.0
        self.RCI = 0.0
        self.RCT = 0.0
        self.ANGI = 0.0
        self.ANGE = 0.0
        self.XW = np.zeros(81)
        self.YW = np.zeros(81)
        self.XWI = np.zeros(81)
        self.YWI = np.zeros(81)

        # Variables in GCB common block
        self.NGCB = 0
        self.NXNYCB = np.zeros(81, dtype=int)
        self.NCBPTS = 0
        self.IINTCB = 0
        self.IDIFCB = 0
        self.LECB = 0
        self.XICB = 0.0
        self.RICB = 0.0
        self.XTCB = 0.0
        self.RTCB = 0.0
        self.XECB = 0.0
        self.RECB = 0.0
        self.RCICB = 0.0
        self.RCTCB = 0.0
        self.ANGICB = 0.0
        self.ANGECB = 0.0
        self.XCB = np.zeros(81)
        self.YCB = np.zeros(81)
        self.XCBI = np.zeros(81)
        self.YCBI = np.zeros(81)

        # Variables in BCC common block
        self.PT = np.zeros(21)
        self.TT = np.zeros(21)
        self.THETA = np.zeros(21)
        self.PE = 0.0
        self.MASSE = 0.0
        self.MASSI = 0.0
        self.MASST = 0.0
        self.THRUST = 0.0
        self.NSTAG = 0

def main():
    # Create an instance of GlobalData
    data = GlobalData()

    # Program initialization
    print("**************************************************************")
    print("NAP: Two-Dimensional Time-Dependent Inviscid Nozzle Flow")
    print("**************************************************************")
    print("Program starting...")

    # Example: Set some initial values and perform calculations
    data.GAM1 = data.GAMMA - 1
    data.GAM2 = 2 / data.GAMMA
    data.QUT[0, 0] = data.GAMMA * data.RGAS

    # Add additional program logic here

    print("Program completed.")

if __name__ == "__main__":
    main()
