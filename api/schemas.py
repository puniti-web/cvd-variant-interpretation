"""Pydantic models for the read-only API. Field names follow api/API_CONTRACT.md."""

from __future__ import annotations

from typing import Any, Generic, Literal, TypeVar

from pydantic import BaseModel, Field

T = TypeVar("T")


class Maybe(BaseModel, Generic[T]):
    value: T | None = None
    reason: str | None = None


def maybe(value: T | None, reason: str | None = None) -> dict[str, Any]:
    """Null measurements carry a reason. A present value clears the reason."""
    if value is None:
        return {"value": None, "reason": reason or "Not computed."}
    return {"value": value, "reason": None}


def absent(reason: str) -> dict[str, Any]:
    return maybe(None, reason)


class ArtifactModel(BaseModel):
    key: str
    path: str
    present: bool
    loaded: bool
    reason: str | None = None


class HealthResponse(BaseModel):
    status: Literal["ok"] = "ok"
    artifacts: list[ArtifactModel]


class PredictionModel(BaseModel):
    label: Literal["pathogenic", "benign"] | None = None
    prob_pathogenic: float | None = None
    threshold: Literal[0.5] = 0.5
    reason: str | None = None


class ShapFeatureModel(BaseModel):
    feature: str
    raw_feature: str
    shap_value: float
    feature_value_transformed: float
    feature_value_raw: float | str | bool | None = None


class LocalShapModel(BaseModel):
    shap_space: Literal["probability"] = "probability"
    base_value: float | None = None
    output_value: float | None = None
    features: list[ShapFeatureModel] = Field(default_factory=list)
    reason: str | None = None


class ErrorDetail(BaseModel):
    value: None = None
    reason: str


class ErrorBody(BaseModel):
    detail: ErrorDetail
