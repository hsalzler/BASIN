#!/usr/bin/env python3
"""Regression test for HARMONY downstream geometry command construction.

No BAM is read. This verifies the production command builder used by both the
main HARMONY workflow and regenerate_downstream.py.
"""
from pathlib import Path
import harmony

BAM=Path("sample.bam")
BW=Path("sample.bw")
BASE={"min_mapq":0, "downstream_bin_size":1}

def cmd(center, extend, representation="downstream"):
    cfg = dict(BASE, downstream_center_reads=center, downstream_extend_reads=extend)
    return harmony.bamcoverage_cmd(BAM, BW, cfg, norm="NONE", scale_factor=2.5,
                                   representation=representation)

def has(c, flag): return flag in c

# ordinary aligned-read derivative
c=cmd(False, False)
assert not has(c,"--centerReads") and not has(c,"--extendReads")

# requested read-length footprint centered on PE fragment: deepTools needs BOTH flags
c=cmd(True, False)
assert has(c,"--centerReads") and has(c,"--extendReads")

# requested full fragment-spanning derivative
c=cmd(False, True)
assert not has(c,"--centerReads") and has(c,"--extendReads")

# ambiguous user request must fail
try:
    cmd(True, True)
except SystemExit:
    pass
else:
    raise AssertionError("center=true, extend=true should be rejected")

# canonical path must remain immutable regardless of downstream config
for center, extend in [(False,False),(True,False),(False,True),(True,True)]:
    c=cmd(center, extend, representation="canonical")
    assert not has(c,"--centerReads") and not has(c,"--extendReads")

print("PASS: HARMONY downstream geometry contract")
