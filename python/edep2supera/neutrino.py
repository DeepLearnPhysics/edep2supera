"""Extract neutrino interaction truth from EDepSim primary vertices."""

from dataclasses import dataclass
import math
import re


NEUTRINO_PDGS = {12, 14, 16}
GENIE_MODES = {
    "QES": 1,
    "1Kaon": 2,
    "DIS": 3,
    "RES": 4,
    "COH": 5,
    "DFR": 6,
    "NuEEL": 7,
    "IMD": 8,
    "AMNuGamma": 9,
    "MEC": 10,
    "CEvNS": 11,
    "IBD": 12,
    "GLR": 13,
    "IMDAnh": 14,
    "PhotonCOH": 15,
    "PhotonRES": 16,
    "1Pion": 17,
}


@dataclass
class NeutrinoTruth:
    index: int
    interaction_id: int
    position: tuple
    pdg_code: int
    momentum: tuple
    energy_init: float
    target: int = -1
    nucleon: int = -1
    current_type: int = -1
    interaction_mode: int = -1
    interaction_type: int = -1
    lepton_pdg_code: int = 0
    lepton_track_id: int = -1
    lepton_p: float = 0.0
    theta: float = 0.0
    momentum_transfer: float = 0.0
    momentum_transfer_mag: float = 0.0
    energy_transfer: float = 0.0
    bjorken_x: float = 0.0
    inelasticity: float = 0.0
    hadronic_invariant_mass: float = 0.0
    reaction: str = ""


def _get(obj, accessor, field):
    method = getattr(obj, accessor, None)
    return method() if callable(method) else getattr(obj, field)


def _particle_pdg(particle):
    return int(_get(particle, "GetPDGCode", "PDGCode"))


def _particle_momentum(particle):
    return _get(particle, "GetMomentum", "Momentum")


def _components(momentum):
    return (
        float(momentum.X()),
        float(momentum.Y()),
        float(momentum.Z()),
        float(momentum.E()),
    )


def _parse_reaction(reaction):
    fields = {}
    for key, value in re.findall(r"(?:^|;)([^:;]+):([^;]+)", reaction):
        fields[key] = value

    process = fields.get("proc", "")
    current_type = 0 if "Weak[CC]" in process else 1 if "Weak[NC]" in process else -1
    mode_name = process.rsplit(",", 1)[-1] if "," in process else ""
    mode = GENIE_MODES.get(mode_name, -1)
    def integer(name):
        try:
            return int(fields.get(name, -1))
        except (TypeError, ValueError):
            return -1

    return {
        "target": integer("tgt"),
        "nucleon": integer("N"),
        "current_type": current_type,
        "interaction_mode": mode,
        # This is the best interaction code present in an EDepSim reaction
        # string. It matches flow2supera's fallback when no separate code exists.
        "interaction_type": mode,
    }


def _expected_lepton(neutrino_pdg, current_type):
    if current_type == 0:
        return (1 if neutrino_pdg > 0 else -1) * (abs(neutrino_pdg) - 1)
    if current_type == 1:
        return neutrino_pdg
    return None


def _kinematics(neutrino, lepton, nucleon):
    npx, npy, npz, ne = _components(_particle_momentum(neutrino))
    lpx, lpy, lpz, le = _components(_particle_momentum(lepton))
    npmag = math.sqrt(npx * npx + npy * npy + npz * npz)
    lpmag = math.sqrt(lpx * lpx + lpy * lpy + lpz * lpz)
    q0 = ne - le
    qx, qy, qz = npx - lpx, npy - lpy, npz - lpz
    q3 = math.sqrt(qx * qx + qy * qy + qz * qz)
    q2 = max(0.0, q3 * q3 - q0 * q0)
    theta = 0.0
    if npmag > 0.0 and lpmag > 0.0:
        cosine = (npx * lpx + npy * lpy + npz * lpz) / (npmag * lpmag)
        theta = math.acos(max(-1.0, min(1.0, cosine)))

    # Use the struck-nucleon mass for the conventional reconstructed x and W.
    mass = 938.272 if nucleon == 2212 else 939.565
    x = q2 / (2.0 * mass * q0) if q0 > 0.0 else 0.0
    y = q0 / ne if ne > 0.0 else 0.0
    w2 = mass * mass + 2.0 * mass * q0 - q2
    return {
        "lepton_p": lpmag,
        "theta": theta,
        "momentum_transfer": q2,
        "momentum_transfer_mag": q3,
        "energy_transfer": q0,
        "bjorken_x": x,
        "inelasticity": y,
        "hadronic_invariant_mass": math.sqrt(max(0.0, w2)),
    }


def neutrinos_from_event(event):
    """Return one :class:`NeutrinoTruth` for each recognizable interaction.

    RooTracker-backed EDepSim files preserve the status-zero incoming neutrino
    and target in an ``initial-state`` informational vertex. Other generators
    become supported automatically when they populate the same TG4 contract.
    """

    result = []
    for vertex in event.Primaries:
        initial = next(
            (
                item
                for item in vertex.Informational
                if str(_get(item, "GetGeneratorName", "GeneratorName"))
                == "initial-state"
            ),
            None,
        )
        if initial is None:
            continue
        neutrino = next(
            (p for p in initial.Particles if abs(_particle_pdg(p)) in NEUTRINO_PDGS),
            None,
        )
        if neutrino is None:
            continue

        reaction = str(_get(vertex, "GetReaction", "Reaction"))
        parsed = _parse_reaction(reaction)
        nu_pdg = _particle_pdg(neutrino)
        expected = _expected_lepton(nu_pdg, parsed["current_type"])
        lepton = next(
            (p for p in vertex.Particles if _particle_pdg(p) == expected), None
        )
        position = _get(vertex, "GetPosition", "Position")
        px, py, pz, energy = _components(_particle_momentum(neutrino))
        truth = NeutrinoTruth(
            index=len(result),
            interaction_id=int(
                _get(vertex, "GetInteractionNumber", "InteractionNumber")
            ),
            # EDepSim persists positions in mm; LArCV uses cm. Both use ns here.
            position=(
                float(position.X()) / 10.0,
                float(position.Y()) / 10.0,
                float(position.Z()) / 10.0,
                float(position.T()),
            ),
            pdg_code=nu_pdg,
            momentum=(px, py, pz),
            energy_init=energy,
            reaction=reaction,
            **parsed,
        )
        if lepton is not None:
            truth.lepton_pdg_code = _particle_pdg(lepton)
            truth.lepton_track_id = int(
                _get(lepton, "GetTrackId", "TrackId")
            )
            for key, value in _kinematics(neutrino, lepton, truth.nucleon).items():
                setattr(truth, key, value)
        result.append(truth)
    return result
