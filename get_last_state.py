#!/usr/bin/env python
#
# Written by Satoshi Ono
#

import argparse
import panedr
import kkkit
import rangeinfo

parser = argparse.ArgumentParser(description='Get last vstate from edr.')
parser.add_argument('-f', dest='edr', nargs='?',
                    required=True, type=str,
                    help='EDR trajectory')
parser.add_argument('-t', dest='temp', nargs='?',
                    required=True, type=str,
                    help='VMCMD paramaeter file')
args = parser.parse_args()

df = panedr.edr_to_df('{}'.format(args.edr), verbose=False)
#last_pot = df['Potential'].tolist()[-1]
last_pot = df['Potential'].tail(1).tolist()[0]

ri = rangeinfo.RangeInfo('{}'.format(args.temp))
ri.read_ttp_inp()

last_vst = ri.get_near_state(last_pot)
print(last_vst)

