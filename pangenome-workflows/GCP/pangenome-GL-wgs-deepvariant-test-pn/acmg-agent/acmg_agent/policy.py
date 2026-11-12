"""Reviewer output policy: the LLM describes evidence and criteria only.

It must not recommend testing, counselling, referral, management, treatment,
surveillance or any other next clinical step. Those decisions belong to the
clinical team, and nothing in the evidence pack supports them.
"""
from __future__ import annotations

import re

RECOMMENDATION_RULE = (
    "Describe the evidence and the criteria only. Do NOT recommend or suggest any "
    "clinical action: no testing (predictive, confirmatory, carrier, cascade or "
    "family), no genetic counselling, referral, management, treatment, surveillance, "
    "screening or follow-up, and no phrases such as 'should be considered' or "
    "'is advised'. Those decisions belong to the clinical team."
)

REDACTION_NOTE = "Sentences containing clinical recommendations were removed (reviewer output policy)."

_PATTERNS = [
    r"\brecommend\w*",
    r"\badvis(?:e|ed|able)\b",
    r"\bsuggest(?:s|ed)?\b[^.]*\b(?:test\w*|counsel\w*|referr?\w*|evaluat\w*|screen\w*)",
    r"\bshould\b[^.]*\b(?:be\s+)?(?:offered|considered|referred|tested|counsel\w*|screened|evaluated|monitored|followed)\b",
    r"\b(?:consider|warrant\w*|merit\w*)\b[^.]*\b(?:test\w*|counsel\w*|referr?\w*|screen\w*|surveillance|evaluation)\b",
    r"\bcounsel(?:l)?ing\b",
    r"\bcounsel(?:l)?ed\b",
    r"\breferr(?:al|ed)\b",
    r"\b(?:predictive|presymptomatic|pre-symptomatic|cascade|carrier|confirmatory|prenatal|preimplantation)\s+(?:genetic\s+)?testing\b",
    r"\b(?:family|relatives?)\b[^.]*\btesting\b",
    r"\bsurveillance\b",
    r"\b(?:clinical|medical|disease|patient)\s+management\b",
    r"\btreatment\s+(?:options?|should|is\s+(?:advised|indicated))\b",
    r"\b(?:initiat|start|begin)\w*\s+(?:\w+\s+)?treatment\b",
    r"\bfollow[- ]up\b",
    r"\bamenable to\b",
]
_RE = re.compile("|".join(f"(?:{p})" for p in _PATTERNS), re.IGNORECASE)
_SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9\"'(])")


def split_sentences(text: str) -> list[str]:
    return [s for s in _SENT_SPLIT.split((text or "").strip()) if s.strip()]


def find_recommendations(text: str | None) -> list[str]:
    """Sentences in ``text`` that recommend a clinical action."""
    return [s for s in split_sentences(text or "") if _RE.search(s)]


def strip_recommendations(text: str | None) -> tuple[str | None, bool]:
    """Remove offending sentences. Returns (clean_text_or_None, removed_any)."""
    if not text:
        return text, False
    sents = split_sentences(text)
    keep = [s for s in sents if not _RE.search(s)]
    if len(keep) == len(sents):
        return text, False
    clean = " ".join(keep).strip()
    return (clean or None), True


def violation_message(sentences: list[str]) -> str:
    quoted = " | ".join(f'"{s}"' for s in sentences)
    return (
        "POLICY VIOLATION: these sentences recommend a clinical action, which is not allowed: "
        f"{quoted}. {RECOMMENDATION_RULE} Rewrite the complete text without them."
    )
