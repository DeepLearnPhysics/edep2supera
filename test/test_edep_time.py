import pytest
import ROOT

from edep2supera import edep2supera


def test_make_edeps_interpolates_time_in_microseconds():
    hit = ROOT.TG4HitSegment()
    hit.GetStart().SetXYZT(0.0, 0.0, 0.0, 1000.0)
    hit.GetStop().SetXYZT(0.59, 0.0, 0.0, 4000.0)

    edeps = edep2supera.SuperaDriver().MakeEDeps(hit)

    assert len(edeps) == 2
    assert [edep.t for edep in edeps] == pytest.approx([1.75, 3.25])
