import random
import csv
from itertools import product

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    Uses guided random (Monte Carlo weighted by RSCU) with a sliding window fallback
    to satisfy hairpin, forbidden sequence, internal promoter, CAI, and internal RBS constraints.
    """

    MAX_ATTEMPTS = 1000

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.aaToSynonymousCodons = {}
        self.aaToWeights = {}
        self.codonToFreq = {}
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.codonChecker = None
        self.internalRbsChecker = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

        # Build RSCU-weighted synonymous codon tables from codon_usage.txt
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        aa_codon_freq = {}
        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 3:
                    continue
                codon = row[0].strip()
                aa = row[1].strip()
                freq = float(row[2].strip())
                if aa == '*':
                    continue
                aa_codon_freq.setdefault(aa, []).append((codon, freq))
                self.codonToFreq[codon] = freq

        for aa, codon_freqs in aa_codon_freq.items():
            self.aaToSynonymousCodons[aa] = [cf[0] for cf in codon_freqs]
            self.aaToWeights[aa] = [cf[1] for cf in codon_freqs]

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()
        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()
        self.internalRbsChecker = InternalRBSChecker()
        self.internalRbsChecker.initiate()

    def _passesAll(self, codons: list, utr: str) -> bool:
        """Returns True if the codon list satisfies all constraints."""
        cds = ''.join(codons)
        full_seq = utr.upper() + cds
        ok, _ = self.forbiddenChecker.run(full_seq)
        if not ok:
            return False
        ok, _ = self.promoterChecker.run(full_seq)
        if not ok:
            return False
        ok, _ = hairpin_checker(full_seq)
        if not ok:
            return False
        ok, _, _, _ = self.codonChecker.run(codons)
        if not ok:
            return False
        ok, _ = self.internalRbsChecker.run(cds)
        if not ok:
            return False
        return True

    def _scoreWindow(self, candidate_cds: str, utr: str) -> int:
        """Counts satisfied constraints for a candidate CDS (used in sliding window)."""
        score = 0
        full_seq = utr.upper() + candidate_cds
        ok, _ = self.forbiddenChecker.run(full_seq)
        if ok:
            score += 1
        ok, _ = self.promoterChecker.run(full_seq)
        if ok:
            score += 1
        ok, _ = hairpin_checker(full_seq)
        if ok:
            score += 1
        ok, _ = self.internalRbsChecker.run(candidate_cds)
        if ok:
            score += 1
        return score

    def _guidedRandom(self, peptide: str, utr: str):
        """
        Monte Carlo: sample codons weighted by RSCU, accept first fully valid sequence.
        Returns codon list (with stop codon) or None after MAX_ATTEMPTS.
        """
        for _ in range(self.MAX_ATTEMPTS):
            codons = [
                random.choices(self.aaToSynonymousCodons[aa], weights=self.aaToWeights[aa], k=1)[0]
                for aa in peptide
            ]
            codons.append('TAA')
            if self._passesAll(codons, utr):
                return codons
        return None

    def _slidingWindow(self, peptide: str, utr: str) -> list:
        """
        Greedy sliding window over non-overlapping 3-AA windows.
        For each window, enumerates all synonymous codon combos scored in context
        (preamble + window + downstream seed for next 6 AAs). Picks best by
        constraint score, then average RSCU as a tie-breaker.
        """
        n = len(peptide)
        codons = [self.aminoAcidToCodon[aa] for aa in peptide] + ['TAA']

        for i in range(0, n, 3):
            window_aas = peptide[i:i + 3]
            if not window_aas:
                break
            downstream = [self.aminoAcidToCodon[aa] for aa in peptide[i + 3:i + 9]]
            codon_options = [self.aaToSynonymousCodons[aa] for aa in window_aas]

            best_score = -1
            best_cai = -1.0
            best_combo = [self.aminoAcidToCodon[aa] for aa in window_aas]

            for combo in product(*codon_options):
                candidate_cds = ''.join(codons[:i] + list(combo) + downstream + ['TAA'])
                score = self._scoreWindow(candidate_cds, utr)
                cai = sum(self.codonToFreq.get(c, 0.01) for c in combo) / len(combo)
                if score > best_score or (score == best_score and cai > best_cai):
                    best_score = score
                    best_cai = cai
                    best_combo = list(combo)

            codons[i:i + len(window_aas)] = best_combo

        return codons

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # Get a tentative RBS for UTR context during codon optimisation
        dummyCds = ''.join(self.aminoAcidToCodon[aa] for aa in peptide) + 'TAA'
        tentativeRbs = self.rbsChooser.run(dummyCds, ignores)
        utr = tentativeRbs.utr

        # Phase 1: guided random
        codons = self._guidedRandom(peptide, utr)

        # Phase 2: sliding window fallback
        if codons is None:
            codons = self._slidingWindow(peptide, utr)

        # Select final RBS with optimised CDS
        cds = ''.join(codons)
        selectedRBS = self.rbsChooser.run(cds, ignores)

        return Transcript(selectedRBS, peptide, codons)


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    # Print out the transcript information
    print(transcript)
