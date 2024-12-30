#!/bin/bash

sed -i -z 's/\n\[ moleculetype \]/; Include Position restraint file\n#ifdef POSRES\n#include \"posre.itp\"\n#endif\n\n\[ moleculetype \]/2' $1
