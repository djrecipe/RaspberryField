#! /usr/bin/env python
from subprocess import call
call(['./spectrometer matrix.cfg 2>/dev/null'], shell=True)