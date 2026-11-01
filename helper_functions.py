import math 

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
    # Keep only positions where both sequences have a determined nucleotide
    # (skips indels "-" and ambiguous/unsequenced bases such as "N")
    valid = set("ACGTU")
    purines = set("AG")
    pairs = [(a.replace("U", "T"), b.replace("U", "T"))
             for a, b in zip(reference_sequence.upper(), distant_sequence.upper())
             if a in valid and b in valid]
    length = len(pairs)
    if length == 0:
        return 0.0
    transitions = sum(1 for a, b in pairs if a != b and ((a in purines) == (b in purines)))
    transversions = sum(1 for a, b in pairs if a != b and ((a in purines) != (b in purines)))
    if transitions == 0 and transversions == 0:
        return 0.0
    P = transitions / length       # proportion of transitions (A<->G, C<->T)
    Q = transversions / length     # proportion of transversions (purine <-> pyrimidine)
    arg1 = 1 - 2 * P - Q
    arg2 = 1 - 2 * Q
    if arg1 <= 0 or arg2 <= 0:
        return float("inf")
    d = -0.5 * math.log(arg1) - 0.25 * math.log(arg2)   # substitutions per site
    return d * length                                    # per genome, same unit as uncorrected
