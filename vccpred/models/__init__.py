"""Model registry: add a new model by importing it and listing it in MODELS."""
from vccpred.models.base import PerturbationModel
from vccpred.models.control_resample import ControlResample

MODELS: dict[str, type[PerturbationModel]] = {m.name: m for m in [ControlResample]}


def get_model(name: str, **params) -> PerturbationModel:
    if name not in MODELS:
        raise KeyError(f"unknown model {name!r}; available: {sorted(MODELS)}")
    return MODELS[name](**params)
