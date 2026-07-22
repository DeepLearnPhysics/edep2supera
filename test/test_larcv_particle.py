from types import SimpleNamespace
from unittest.mock import patch

from edep2supera import utils


class FakeParticle:
    def __init__(self):
        self.values = {}

    def __getattr__(self, name):
        def setter(*values):
            self.values[name] = values[0] if len(values) == 1 else values

        return setter


class FakeEnergy:
    def sum(self):
        return 1.0

    def size(self):
        return 1


def test_larcv_particle_stores_ancestor_track_id_not_instance_id():
    output = FakeParticle()

    vertex = SimpleNamespace(
        pos=SimpleNamespace(x=0.0, y=0.0, z=0.0), time=0.0
    )
    part = SimpleNamespace(
        id=11,
        shape=0,
        trackid=22,
        genid=0,
        pdg=2212,
        px=0.0,
        py=0.0,
        pz=1.0,
        end_px=0.0,
        end_py=0.0,
        end_pz=0.0,
        dist_travel=1.0,
        vtx=vertex,
        end_pt=vertex,
        first_step=vertex,
        last_step=vertex,
        parent_vtx=vertex,
        ancestor_vtx=vertex,
        energy_init=1000.0,
        process="primary",
        parent_trackid=22,
        parent_pdg=2212,
        parent_process="primary",
        parent_id=11,
        children_id=[],
        ancestor_id=11,
        ancestor_trackid=22,
        ancestor_pdg=2212,
        ancestor_process="primary",
        group_id=11,
        interaction_id=0,
    )
    particle = SimpleNamespace(
        part=part,
        energy=FakeEnergy(),
    )

    with patch.object(utils.larcv, "Particle", return_value=output):
        utils.larcv_particle(particle)

    assert output.values["ancestor_track_id"] == 22
    assert output.values["ancestor_track_id"] != part.ancestor_id
