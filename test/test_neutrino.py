import math
from types import SimpleNamespace
from unittest.mock import patch

from edep2supera import utils
from edep2supera.neutrino import neutrinos_from_event


class FourVector:
    def __init__(self, x, y, z, t):
        self.values = (x, y, z, t)

    def X(self):
        return self.values[0]

    def Y(self):
        return self.values[1]

    def Z(self):
        return self.values[2]

    def T(self):
        return self.values[3]

    def E(self):
        return self.values[3]


class FakeNeutrino:
    def __init__(self):
        self.values = {}

    def __getattr__(self, name):
        def setter(*values):
            self.values[name] = values[0] if len(values) == 1 else values

        return setter


def particle(pdg, momentum, track_id=-1):
    return SimpleNamespace(PDGCode=pdg, Momentum=FourVector(*momentum), TrackId=track_id)


def test_extracts_genie_neutrino_and_derived_kinematics():
    incoming = particle(14, (0.0, 0.0, 1000.0, 1000.0))
    target = particle(1000180400, (0.0, 0.0, 0.0, 37224.0))
    lepton = particle(13, (300.0, 0.0, 800.0, 860.0), track_id=4)
    initial = SimpleNamespace(
        GeneratorName="initial-state", Particles=[incoming, target]
    )
    vertex = SimpleNamespace(
        Informational=[initial],
        Particles=[lepton, particle(2212, (0.0, 0.0, 200.0, 960.0))],
        Reaction="nu:14;tgt:1000180400;N:2212;proc:Weak[CC],QES;",
        InteractionNumber=7,
        Position=FourVector(10.0, 20.0, 30.0, 40.0),
    )

    result = neutrinos_from_event(SimpleNamespace(Primaries=[vertex]))

    assert len(result) == 1
    truth = result[0]
    assert truth.interaction_id == 0
    assert truth.position == (1.0, 2.0, 3.0, 40.0)
    assert truth.pdg_code == 14
    assert truth.target == 1000180400
    assert truth.nucleon == 2212
    assert truth.current_type == 0
    assert truth.interaction_mode == 1
    assert truth.lepton_pdg_code == 13
    assert truth.lepton_track_id == 4
    assert truth.energy_transfer == 140.0
    assert truth.momentum_transfer_mag == math.sqrt(130000.0)
    assert truth.momentum_transfer == 110400.0
    assert truth.inelasticity == 0.14


def test_ignores_primary_without_preserved_initial_state():
    event = SimpleNamespace(
        Primaries=[SimpleNamespace(Informational=[], Particles=[])]
    )
    assert neutrinos_from_event(event) == []


def test_tolerates_non_genie_reaction_metadata():
    incoming = particle(12, (0.0, 0.0, 200.0, 200.0))
    initial = SimpleNamespace(GeneratorName="initial-state", Particles=[incoming])
    vertex = SimpleNamespace(
        Informational=[initial],
        Particles=[],
        Reaction="generator-specific text",
        InteractionNumber=1,
        Position=FourVector(0.0, 0.0, 0.0, 0.0),
    )

    truth = neutrinos_from_event(SimpleNamespace(Primaries=[vertex]))[0]
    assert truth.pdg_code == 12
    assert truth.target == -1
    assert truth.current_type == -1
    assert truth.interaction_mode == -1


def test_converts_neutrino_truth_to_larcv():
    incoming = particle(-14, (0.0, 0.0, 500.0, 500.0))
    outgoing = particle(-14, (0.0, 100.0, 400.0, 420.0), track_id=2)
    initial = SimpleNamespace(GeneratorName="initial-state", Particles=[incoming])
    vertex = SimpleNamespace(
        Informational=[initial],
        Particles=[outgoing],
        Reaction="nu:-14;tgt:1000180400;N:2112;proc:Weak[NC],RES;",
        InteractionNumber=3,
        Position=FourVector(0.0, 0.0, 0.0, 1.0),
    )
    truth = neutrinos_from_event(SimpleNamespace(Primaries=[vertex]))[0]
    output = FakeNeutrino()

    with patch.object(utils.larcv, "Neutrino", return_value=output), patch.object(
        utils.larcv, "InstanceID_t", side_effect=lambda value: value
    ):
        utils.larcv_neutrino(truth)

    assert output.values["id"] == 0
    assert output.values["interaction_id"] == 0
    assert output.values["current_type"] == 1
    assert output.values["interaction_mode"] == 4
    assert output.values["pdg_code"] == -14
    assert output.values["lepton_pdg_code"] == -14
    assert output.values["position"] == (0.0, 0.0, 0.0, 1.0)


def test_interaction_ids_are_event_local_and_ordered():
    def vertex(source_id, pdg=14):
        incoming = particle(pdg, (0.0, 0.0, 1000.0, 1000.0))
        initial = SimpleNamespace(
            GeneratorName="initial-state", Particles=[incoming]
        )
        return SimpleNamespace(
            Informational=[initial],
            Particles=[],
            Reaction="nu:14;tgt:1000180400;proc:Weak[NC],QES;",
            InteractionNumber=source_id,
            Position=FourVector(0.0, 0.0, 0.0, 0.0),
        )

    first_event = neutrinos_from_event(
        SimpleNamespace(Primaries=[vertex(17), vertex(103)])
    )
    second_event = neutrinos_from_event(SimpleNamespace(Primaries=[vertex(104)]))

    assert [truth.interaction_id for truth in first_event] == [0, 1]
    assert [truth.interaction_id for truth in second_event] == [0]
