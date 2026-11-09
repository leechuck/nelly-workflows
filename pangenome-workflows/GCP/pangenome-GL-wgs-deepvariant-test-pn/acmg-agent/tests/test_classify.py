from acmg_agent.criteria import Criterion
from acmg_agent.classify import classify, points, rule_based_class, point_class


def test_points_hexa():
    cs = [
        Criterion("PVS1", "met", "very_strong", "null"),
        Criterion("PM2", "met", "supporting", "rare"),
        Criterion("PM3", "met", "supporting", "hom in AR (0.5 pt)"),
    ]
    assert points(cs) == 10
    out = classify(cs)
    assert out.class_ == "Pathogenic"
    assert out.rule_based_class == "Pathogenic"
    assert out.points == 10
    assert out.conflict is False


def test_rule_based_pathogenic_pvs1_plus_pm():
    cs = [
        Criterion("PVS1", "met", "very_strong", "null"),
        Criterion("PM3", "met", "moderate", "in trans with P"),
    ]
    out = classify(cs)
    assert out.class_ in ("Likely pathogenic", "Pathogenic")


def test_conflict_flag():
    cs = [
        Criterion("PVS1", "met", "very_strong", "null"),
        Criterion("PM2", "met", "supporting", "rare"),
        Criterion("BS2", "met", "strong", "seen in controls"),
    ]
    out = classify(cs)
    assert out.conflict is True


def test_ba1_benign():
    cs = [Criterion("BA1", "met", "stand_alone", "AF 0.12")]
    out = classify(cs)
    assert out.class_ == "Benign"


def test_point_class_boundaries():
    assert point_class(10) == "Pathogenic"
    assert point_class(6) == "Likely pathogenic"
    assert point_class(5) == "VUS"
    assert point_class(0) == "VUS"
    assert point_class(-1) == "Likely benign"
    assert point_class(-5) == "Likely benign"
    assert point_class(-6) == "Benign"
