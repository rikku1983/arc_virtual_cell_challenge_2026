"""Prediction pipeline for the Virtual Cell Challenge 2026.

Models live in vccpred.models and share one interface (PerturbationModel);
everything else — loading controls, writing the submission, packaging — is
model-agnostic. Run: .venv/bin/python -m vccpred --help
"""
