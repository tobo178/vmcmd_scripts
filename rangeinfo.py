#
# from omegagene
#

import sys
import re
import os
import kkkit

class RangeInfo(kkkit.FileI):
    """ RangeInfo:
    
    Read TTP-V-MCMD input file
    ri=ReadInfo("temp_s")
    ri.read_ttp_inp()
    len(ri.vs_range.keys())
    len(ri.vs_order)

    """
    def __init__(self, fn):
        super(RangeInfo, self).__init__(fn)
        ## ene_range ... min, max
        self.ene_range = (0.0, 0.0, 100.0)
        self.ene_width = 0.0
        self.vinterval = 20000
        self.temp = 300.0
        ## vs_range[vs_id] = (min_ene, max_ene, tpro1, tpro2)
        self.vs_range = {}
        self.vs_order = {}
        ## vs_params[vs_id] = (coeff[0], ..., coeff[n+1], alphalw, alphaup)
        self.vs_params = {}
    def read_ttp_inp(self):
        self.open()
        line = self.f.readline().strip()
        n_vs = int(line.split()[0])
        line = self.f.readline().strip()
        vinterval = int(line.split()[0])
        self.vinterval = vinterval
        self.vs_range = {}
        for vsid in range(n_vs):
            line = self.f.readline().strip()
            terms = line.split()
            ene_min = float(terms[0])
            ene_max = float(terms[1])
            line = self.f.readline().strip()
            terms = line.split()
            p1 = float(terms[0])
            p2 = float(terms[1])
            self.vs_range[vsid] = (ene_min, ene_max, p1, p2)
        for vsid in range(n_vs):
            self.vs_params[vsid] = []
            self.vs_order[vsid] = []
            order = int(self.f.readline().strip())
            self.vs_order[vsid].append(order)
            for i in range(order+3):
                self.vs_params[vsid].append(self.f.readline().strip())
        self.temp = float(self.f.readline().strip())
        self.close()
        return 0

    def read(self):
        """Range dividing ratios, like
-101700.0 -38000.0 100.0
7
0 0.00 0.10 0.0 1.0
1 0.05 0.20 1.0 1.0
2 0.10 0.30 1.0 1.0
3 0.20 0.45 1.0 1.0
4 0.30 0.60 1.0 1.0
5 0.45 0.80 1.0 1.0
6 0.60 1.00 1.0 0.0

        ri_temp=RangeInfo('test_range.dat')
        ri_temp.read()
        ri_temp.vs_range[0]
        """
        self.open()

        line = self.f.readline().strip()
        terms = line.split()
        self.ene_range = (float(terms[0]), float(terms[1]), float(terms[2]))
        self.ene_width = self.ene_range[1] - self.ene_range[0]
        
        line = self.f.readline().strip()
        n_vs = int(line.split()[0])
        e_width = self.ene_width
        b_width = self.ene_range[2]
        e_min = self.ene_range[0]
        
        self.vs_range = {}
        for i in range(0,n_vs):
            line = self.f.readline().strip()
            terms = line.split()
            vsid = int(terms[0])
            ratio_min = float(terms[1])
            ratio_max = float(terms[2])
            ene_min = int((e_min + e_width * ratio_min)/b_width)*b_width
            ene_max = int((e_min + e_width * ratio_max)/b_width)*b_width
            self.vs_range[vsid] = (ene_min, ene_max, float(terms[3]),
                                   float(terms[4]))
        self.close()
        return 0

    def get_near_state(self, ene):
        # near_vsid = -1
        near_vsid = 0
        for vsid, ene_range in self.vs_range.items():
            if ene_range[0] <= ene:
                near_vsid = vsid
        return near_vsid

