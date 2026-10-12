"""Resource manifest loader (data/resource_manifest.yaml; schema 2)."""
from __future__ import annotations
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import yaml

DEFAULT_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "resource_manifest.yaml")


@dataclass
class Resource:
    id: str
    name: str
    roles: List[str]
    provides: List[str]
    default_enabled: bool
    access: str
    automated: bool
    get: str
    url: str
    format: Optional[str]
    cite: Optional[str]
    terms: Optional[str]
    warn: Optional[str]
    check: dict
    if_missing: dict
    fetch: dict
    flags: dict
    raw: dict = field(repr=False, default_factory=dict)


@dataclass
class Manifest:
    schema: int
    profiles: dict
    resources: Dict[str, Resource]
    metric_requirements: dict
    policy: dict
    use_cases: list

    def by_capability(self, cap: str) -> List[str]:
        return [r.id for r in self.resources.values() if cap in r.provides]


def load_manifest(path: Optional[str] = None) -> Manifest:
    d = yaml.safe_load(open(path or DEFAULT_PATH))
    res = {}
    for r in d["resources"]:
        role = r.get("role", [])
        res[r["id"]] = Resource(
            id=r["id"], name=r.get("name", r["id"]), roles=[role] if isinstance(role, str) else list(role),
            provides=list(r.get("provides", [])), default_enabled=bool(r.get("default_enabled", False)),
            access=r.get("access", "public_download"), automated=bool(r.get("automated", False)), get=r.get("get", ""),
            url=r.get("url", ""), format=r.get("format"), cite=r.get("cite"), terms=r.get("terms"), warn=r.get("warn"),
            check=r.get("check") or {}, if_missing=r.get("if_missing") or {}, fetch=r.get("fetch") or {}, flags=r.get("flags") or {}, raw=r)
    return Manifest(d.get("schema", 1), d.get("profiles", {}), res, d.get("metric_requirements", {}),
                    d.get("missing_resource_policy", {}), d.get("use_cases", []))
