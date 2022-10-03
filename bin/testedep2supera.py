#!/usr/bin/env python3
import sys
from ROOT import TChain
from edep2supera import edep2supera

proc = edep2supera.SuperaDriver()

ch = TChain('EDepSimEvents')
print('Adding input:', sys.argv[1])
ch.AddFile(sys.argv[1])
print("Chain has", ch.GetEntries(), "entries")
sys.stdout.flush()

print("step1")
for entry in range(ch.GetEntries()):
    print("considering event:", entry)
    sys.stdout.flush()
    bytes = ch.GetEntry(entry)
    if bytes < 1: break
    ev = ch.Event
    proc.ReadEvent(ev)

    if entry > 4:
        break