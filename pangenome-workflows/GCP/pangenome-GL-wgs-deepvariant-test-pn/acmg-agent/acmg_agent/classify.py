"""ACMG/AMP 2015 Table 5 combining rules + Tavtigian 2020 Bayesian point scale."""
from __future__ import annotations

from acmg_agent.model import Classification, Criterion

_PATH_POINTS = {
    "very_strong": 8,
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
}

_BENIGN_POINTS = {
    "stand_alone": -8,
    "strong": -4,
    "moderate": -2,
    "supporting": -1,
}


def points(criteria: list[Criterion]) -> int:
    """Compute Tavtigian 2020 net Bayesian point score across all met criteria."""
    total = 0
    for c in criteria:
        if c.status != "met" or not c.strength:
            continue
        if c.code.startswith("P"):
            total += _PATH_POINTS.get(c.strength, 0)
        elif c.code.startswith("B"):
            total += _BENIGN_POINTS.get(c.strength, 0)
    return total


def point_class(pts: int, ba1_met: bool = False) -> str:
    """Map Tavtigian 2020 point score to 5-tier ACMG classification."""
    if ba1_met:
        return "Benign"
    if pts >= 10:
        return "Pathogenic"
    if pts >= 6:
        return "Likely pathogenic"
    if pts >= 0:
        return "VUS"
    if pts >= -5:
        return "Likely benign"
    return "Benign"


def rule_based_class(criteria: list[Criterion]) -> str:
    """Evaluate Richards et al. 2015 Table 5 + ClinGen SVI combining rules."""
    p_vs = 0
    p_s = 0
    p_m = 0
    p_p = 0
    b_sa = 0
    b_s = 0
    b_p = 0

    for c in criteria:
        if c.status != "met" or not c.strength:
            continue
        if c.code.startswith("P"):
            if c.strength == "very_strong":
                p_vs += 1
            elif c.strength == "strong":
                p_s += 1
            elif c.strength == "moderate":
                p_m += 1
            elif c.strength == "supporting":
                p_p += 1
        elif c.code.startswith("B"):
            if c.strength == "stand_alone" or c.code == "BA1":
                b_sa += 1
            elif c.strength == "strong":
                b_s += 1
            elif c.strength in ("supporting", "moderate"):
                b_p += 1

    if b_sa >= 1:
        return "Benign"

    # Pathogenic (Richards 2015 Table 5)
    is_path = (
        (p_vs >= 1 and (p_s >= 1 or p_m >= 2 or (p_m >= 1 and p_p >= 1) or p_p >= 2))
        or (p_s >= 2)
        or (p_s == 1 and (p_m >= 3 or (p_m == 2 and p_p >= 2) or (p_m == 1 and p_p >= 4)))
    )

    # Likely pathogenic (Richards 2015 Table 5 + ClinGen SVI PVS1 + PM2_Supporting)
    is_lp = (
        (p_vs == 1 and p_m == 1)
        or (p_vs == 1 and p_p >= 1)
        or (p_s == 1 and 1 <= p_m <= 2)
        or (p_s == 1 and p_p >= 2)
        or (p_m >= 3)
        or (p_m == 2 and p_p >= 2)
        or (p_m == 1 and p_p >= 4)
    )

    # Benign / Likely benign (Richards 2015 Table 5)
    is_ben = b_s >= 2
    is_lb = (b_s == 1 and b_p >= 1) or (b_p >= 2)

    if is_path:
        return "Pathogenic"
    if is_lp:
        return "Likely pathogenic"
    if is_ben:
        return "Benign"
    if is_lb:
        return "Likely benign"
    return "VUS"


def classify(criteria: list[Criterion]) -> Classification:
    """Combine met criteria using both Tavtigian 2020 points (primary) and Richards 2015 rules."""
    ba1_met = any(c.code == "BA1" and c.status == "met" for c in criteria)
    pts = points(criteria)
    p_cls = point_class(pts, ba1_met=ba1_met)
    rb_cls = rule_based_class(criteria)

    has_path = any(c.code.startswith("P") and c.status == "met" for c in criteria)
    has_benign = any(c.code.startswith("B") and c.status == "met" for c in criteria)
    conflict = (p_cls != rb_cls) or (has_path and has_benign)

    return Classification(
        class_=p_cls,
        points=pts,
        rule_based_class=rb_cls,
        conflict=conflict,
    )
