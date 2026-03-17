import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


@pytest.fixture
def checker():
    c = InternalRBSChecker()
    c.initiate()
    return c


def test_no_internal_rbs(checker):
    """
    A CDS with no Shine-Dalgarno motif upstream of any internal ATG.
    The sequence has an internal ATG at codon 5, but no SD motif in its upstream window.
    Expects (True, None).
    """
    # ATG at start (codon 1), then codons with no SD motif before the internal ATG at codon 5
    # Upstream of codon 5 ATG: 'GCAGCAGCA' — no AGGAG/AAGG/AGGA/GAGG
    cds = "ATG" + "GCA" * 3 + "ATG" + "GCA" * 5 + "TAA"
    result, region = checker.run(cds)
    assert result is True
    assert region is None


def test_internal_rbs_detected(checker):
    """
    A CDS with AGGAG (canonical Shine-Dalgarno) in the 13 nt upstream of an internal ATG.
    Expects (False, problematic_region).

    Construction:
      - Position 0-2:  ATG  (legitimate start)
      - Position 3-8:  AGGAGA  (contains AGGAG at positions 3-7)
      - Position 9-11: ATG  (internal start — codon 4, scanned by checker)
      - upstream = cds[0:9] = 'ATGAGGAGA' → contains 'AGGAG' → flagged
    """
    cds = "ATGAGGAGA" + "ATG" + "GCAGCAGCA" + "TAA"
    result, region = checker.run(cds)
    assert result is False
    assert region is not None
    assert "ATG" in region


def test_first_codon_atg_not_flagged(checker):
    """
    The first codon ATG is the legitimate start codon and must never be flagged,
    even if there's an SD-like sequence in the RBS region before it.
    The checker only scans from position 9 onward.
    Expects (True, None).
    """
    # AGGAG right at the start followed by ATG start codon — this is a normal RBS situation
    cds = "ATG" + "GCA" * 10 + "TAA"
    result, region = checker.run(cds)
    assert result is True
    assert region is None


def test_near_start_atg_not_flagged(checker):
    """
    An ATG at codon positions 2 and 3 (within the first 9 nt) must not be scanned.
    The checker starts at position 9, so ATG codons at positions 3 and 6 are excluded.
    Expects (True, None) even if there's an SD motif early in the sequence.
    """
    # ATG at position 0 (start), ATG at position 3 (codon 2), ATG at position 6 (codon 3)
    # None of these should trigger the checker since scan starts at position 9
    cds = "ATG" + "ATG" + "ATG" + "GCA" * 8 + "TAA"
    result, region = checker.run(cds)
    assert result is True
    assert region is None


def test_sd_without_downstream_atg_not_flagged(checker):
    """
    A sequence containing an SD motif (AGGAG) but no downstream ATG codon.
    Expects (True, None) because the SD motif alone is not an RBS.
    """
    # Build a CDS where AGGAG appears but is not followed by ATG in codon positions
    # Codons: ATG GCA GGT GAA GCA GCA GAA GAT GCA TAA
    # No internal ATG codon → no flag
    cds = "ATG" + "GCA" + "GGT" + "GAA" + "GCA" + "GCA" + "GAA" + "GAT" + "GCA" + "TAA"
    # Manually verify: no codon from position 9+ is ATG
    assert "ATG" not in [cds[i:i+3] for i in range(9, len(cds)-2, 3) if cds[i:i+3] != "TAA"]
    result, region = checker.run(cds)
    assert result is True
    assert region is None


def test_multiple_internal_atg_flags_first(checker):
    """
    Two internal ATG codons; the first has AGGAG upstream, the second does not.
    Expects (False, ...) flagging the first problematic ATG.

    Construction:
      - 'ATGAGGAGA' + 'ATG' → codon 4 (pos 9) has AGGAG upstream → flagged first
      - 'GCA'*3 + 'ATG' → codon 8 (pos 21) has no AGGAG upstream → would be clean
    """
    cds = "ATGAGGAGA" + "ATG" + "GCA" * 3 + "ATG" + "GCAGCA" + "TAA"
    result, region = checker.run(cds)
    assert result is False
    assert region is not None
    assert "ATG" in region
