# Configuration file for TTP-V-McMD
#
# Written by Satoshi Ono
#
class vmcmd_config:
    def __init__(self):
        self.gmx = "$HOME/gmx2024.4_vmcmd_cpu/bin/GMXRC.bash"
        self.md_num = 128
        # self.n = 10
        # self.vinterval = (20000,20000,10000,10000,5000,5000,5000,5000,5000,5000)
        # self.nvst = 8
        # self.temp = 1500.0 # or 800.0
