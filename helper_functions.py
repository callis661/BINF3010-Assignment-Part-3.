def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    valid = set("ACGTU")
    pairs = [(a, b) for a, b in zip(reference_sequence.upper(), distant_sequence.upper())
             if a in valid and b in valid]
    length = len(pairs)
    if length == 0:
        return 0.0
    mismatches = sum(1 for a, b in pairs if a != b)
    if mismatches == 0:
        return 0.0
    p = mismatches / length
    arg = 1 - 4 / 3 * p
    if arg <= 0:
        return float("inf")
    d = -0.75 * math.log(arg)       # substitutions per site
    return d * length               # same unit as uncorrected distance (per genome)

    raise NotImplementedError()


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
