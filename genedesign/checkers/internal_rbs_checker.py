class InternalRBSChecker:
    """
    Checks for internal ribosome binding sites (Shine-Dalgarno sequences) within a CDS.

    An internal RBS occurs when a Shine-Dalgarno-like motif appears 5–13 nt upstream
    of an internal ATG codon (i.e., any ATG after the first three codons). This can
    cause unintended translation initiation mid-sequence in E. coli.

    Input (run method):
        cds (str): The coding sequence string.

    Output:
        Tuple[bool, str | None]:
            - (True, None) if no internal RBS is detected.
            - (False, problematic_region) if a Shine-Dalgarno + ATG pair is found.
    """

    SD_PATTERNS = ['AGGAG']

    def initiate(self) -> None:
        """No state to initialize — patterns are class-level constants."""
        pass

    def run(self, cds: str) -> tuple[bool, str | None]:
        """
        Scans the CDS for internal Shine-Dalgarno + ATG combinations.

        Skips the first 9 nt (3 codons) to avoid flagging the legitimate start codon.
        For each subsequent ATG codon, checks 5–13 nt upstream for SD motifs.

        Parameters:
            cds (str): Coding sequence to scan.

        Returns:
            tuple: (True, None) if clean, (False, region) if internal RBS found.
        """
        cds = cds.upper()
        # Start scanning from codon 4 onward (position 9) to skip the real start codon
        for i in range(9, len(cds) - 2, 3):
            if cds[i:i + 3] == 'ATG':
                # Check 5–13 nt upstream for a Shine-Dalgarno motif
                upstream_start = max(0, i - 13)
                upstream = cds[upstream_start:i]
                for pattern in self.SD_PATTERNS:
                    if pattern in upstream:
                        return False, cds[upstream_start:i + 3]
        return True, None
