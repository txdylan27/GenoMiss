"""
Scoring module for gene misannotation detection.

This module provides a composite scoring system (0-100 scale) for evaluating the
likelihood that a fused gene represents a true misannotation event. The score
combines multiple factors with adjustable weights.

Scoring Components:
1. Query Coverage: Percentage of the fused protein sequence covered by alignment
2. Bit Score Improvement: Relative improvement of fused vs control alignments
3. Overlap Equilibrium: Balance of alignment between the two gene parts
4. Organism Hit Count: Number of different organisms matching the fused gene

Weights can be easily adjusted via the coefficient constants below.
"""

# ============================================================================
# SCORING COEFFICIENTS - Adjust these to tune the scoring system
# ============================================================================
# All coefficients should sum to 1.0 for a proper weighted average

WEIGHT_QUERY_COVERAGE = 0.50      # Percentage of fused gene aligned
WEIGHT_BITSCORE_IMPROVEMENT = 0.25  # Bit score improvement over control
WEIGHT_OVERLAP_EQUILIBRIUM = 0.05   # Balance between gene part matches
WEIGHT_ORGANISM_COUNT = 0.15        # Number of organism hits
WEIGHT_EVALUE = 0.05                # E-value significance (lower is better)

# Sanity check for weights
assert abs(WEIGHT_QUERY_COVERAGE + WEIGHT_BITSCORE_IMPROVEMENT +
           WEIGHT_OVERLAP_EQUILIBRIUM + WEIGHT_ORGANISM_COUNT + WEIGHT_EVALUE - 1.0) < 0.001, \
       "Scoring weights must sum to 1.0"

# ============================================================================
# SCORING FUNCTIONS
# ============================================================================

def calculate_query_coverage_score(query_coverage):
    """
    Calculate score based on query coverage percentage.

    Args:
        query_coverage (float): Query coverage percentage (0-100)

    Returns:
        float: Score from 0-100

    Logic:
        - Direct mapping: query_coverage% maps to score 0-100
        - Higher coverage = higher confidence in the fusion
    """
    return min(100.0, max(0.0, query_coverage))


def calculate_bitscore_improvement_score(fused_bitscore, control_bitscore_1, control_bitscore_2):
    """
    Calculate score based on bit score improvement of fused gene vs controls.

    Args:
        fused_bitscore (float): Bit score of the fused protein alignment
        control_bitscore_1 (float): Bit score of gene part 1 alignment
        control_bitscore_2 (float): Bit score of gene part 2 alignment

    Returns:
        float: Score from 0-100

    Logic:
        - Compare fused score to the MAX of the two control scores
        - Calculate relative improvement ratio
        - Cap at 100% improvement for scoring purposes
        - Negative improvements get 0 score

    Examples:
        - Fused: 1000, Control max: 500 → 100% improvement → score 100
        - Fused: 750, Control max: 500 → 50% improvement → score 50
        - Fused: 400, Control max: 500 → negative improvement → score 0
    """
    if control_bitscore_1 is None or control_bitscore_2 is None:
        return 0.0

    max_control = max(control_bitscore_1, control_bitscore_2)

    if max_control == 0:
        return 0.0

    improvement_ratio = (fused_bitscore - max_control) / max_control

    # Map improvement to 0-100 scale
    # 0% improvement → 0 score
    # 100% improvement → 100 score
    # Cap at 100% for scoring purposes
    score = min(100.0, max(0.0, improvement_ratio * 100))

    return score


def calculate_overlap_equilibrium_score(start_query, end_query, gene_1_len):
    """
    Calculate score based on how evenly the alignment spans both gene parts.

    Args:
        start_query (int): Start position of alignment in query sequence
        end_query (int): End position of alignment in query sequence
        gene_1_len (int): Length of the first gene part (fusion junction point)

    Returns:
        float: Score from 0-100

    Logic:
        - Calculate how many AAs of the alignment fall in each gene part
        - Equilibrium ratio = min(part1, part2) / max(part1, part2)
        - 50-50 split → ratio 1.0 → score 100
        - 90-10 split → ratio 0.11 → score 11
        - Favors alignments that genuinely span both genes

    Examples:
        - Gene1: 500 AA, Gene2: 500 AA, Alignment: 250-750 → perfect 50-50 → score 100
        - Gene1: 500 AA, Gene2: 500 AA, Alignment: 450-550 → 50-50 → score 100
        - Gene1: 500 AA, Gene2: 500 AA, Alignment: 100-900 → 40-60 → score 66.7
        - Gene1: 500 AA, Gene2: 500 AA, Alignment: 10-550 → 49-1 → score 2
    """
    # Calculate AA in each gene part
    gene_1_coverage = max(0, min(end_query, gene_1_len) - start_query)
    gene_2_coverage = max(0, end_query - gene_1_len)

    if gene_1_coverage == 0 or gene_2_coverage == 0:
        # Alignment doesn't actually span both genes
        return 0.0

    # Calculate equilibrium ratio
    min_coverage = min(gene_1_coverage, gene_2_coverage)
    max_coverage = max(gene_1_coverage, gene_2_coverage)

    equilibrium_ratio = min_coverage / max_coverage

    # Convert to 0-100 scale
    return equilibrium_ratio * 100


def calculate_organism_count_score(organism_hit_count):
    """
    Calculate score based on number of different organisms with hits.

    Args:
        organism_hit_count (int): Number of different organisms matching this fused gene

    Returns:
        float: Score from 0-100

    Logic:
        - More organisms = higher confidence the fusion is real
        - Logarithmic scaling to avoid over-weighting this factor
        - 1 organism → score 50
        - 5 organisms → score 100
        - 10+ organisms → capped at 100

    Note: This score is calculated separately and passed in, as it requires
          grouping operations across the entire dataset.
    """
    if organism_hit_count <= 0:
        return 0.0

    # Logarithmic scaling with base case: 1 organism = 50 points
    import math
    score = 50 + (50 * math.log(organism_hit_count, 5))

    return min(100.0, max(0.0, score))


def calculate_evalue_score(evalue):
    """
    Calculate score based on e-value (expected value) of alignment.

    Args:
        evalue (float): E-value from DIAMOND alignment

    Returns:
        float: Score from 0-100

    Logic:
        - Lower e-values indicate more significant alignments
        - Logarithmic scaling based on e-value
        - e-value = 0 → score 100
        - e-value = 1e-50 → score 100
        - e-value = 1e-10 → score 70
        - e-value = 1e-5 → score 50 (DIAMOND cutoff)
        - e-value = 1 → score 0
    """
    if evalue <= 0:
        return 100.0

    import math

    # Logarithmic scoring: -log10(evalue) maps to score
    # E-value of 1e-50 gives log value of 50
    # E-value of 1e-10 gives log value of 10
    # E-value of 1e-5 gives log value of 5
    # E-value of 1 gives log value of 0

    log_evalue = -math.log10(evalue)

    # Scale to 0-100: map 0-10 range to 0-100
    # Values below 1e-10 get capped at 100
    score = min(100.0, max(0.0, log_evalue * 10))

    return score


def calculate_composite_score(query_coverage, fused_bitscore, control_bitscore_1,
                              control_bitscore_2, start_query, end_query,
                              gene_1_len, organism_hit_count, evalue):
    """
    Calculate the final composite score combining all factors.

    Args:
        query_coverage (float): Query coverage percentage (0-100)
        fused_bitscore (float): Bit score of fused protein alignment
        control_bitscore_1 (float): Bit score of gene part 1 alignment (or None)
        control_bitscore_2 (float): Bit score of gene part 2 alignment (or None)
        start_query (int): Start position of alignment in query
        end_query (int): End position of alignment in query
        gene_1_len (int): Length of first gene part
        organism_hit_count (int): Number of organisms with hits to this fused gene
        evalue (float): E-value of the alignment

    Returns:
        float: Composite score from 0-100

    Formula:
        score = (w1 × query_coverage_score) +
                (w2 × bitscore_improvement_score) +
                (w3 × overlap_equilibrium_score) +
                (w4 × organism_count_score) +
                (w5 × evalue_score)

    where w1, w2, w3, w4, w5 are the weight coefficients defined at module level.
    """
    # Calculate individual component scores
    coverage_score = calculate_query_coverage_score(query_coverage)
    bitscore_score = calculate_bitscore_improvement_score(
        fused_bitscore, control_bitscore_1, control_bitscore_2
    )
    equilibrium_score = calculate_overlap_equilibrium_score(
        start_query, end_query, gene_1_len
    )
    organism_score = calculate_organism_count_score(organism_hit_count)
    evalue_score = calculate_evalue_score(evalue)

    # Weighted combination
    composite = (
        WEIGHT_QUERY_COVERAGE * coverage_score +
        WEIGHT_BITSCORE_IMPROVEMENT * bitscore_score +
        WEIGHT_OVERLAP_EQUILIBRIUM * equilibrium_score +
        WEIGHT_ORGANISM_COUNT * organism_score +
        WEIGHT_EVALUE * evalue_score
    )

    return round(composite, 2)


def get_scoring_info():
    """
    Get a dictionary with current scoring configuration.

    Returns:
        dict: Scoring weights and methodology info
    """
    return {
        "weights": {
            "query_coverage": WEIGHT_QUERY_COVERAGE,
            "bitscore_improvement": WEIGHT_BITSCORE_IMPROVEMENT,
            "overlap_equilibrium": WEIGHT_OVERLAP_EQUILIBRIUM,
            "organism_count": WEIGHT_ORGANISM_COUNT,
            "evalue": WEIGHT_EVALUE
        },
        "descriptions": {
            "query_coverage": "Percentage of fused gene sequence covered by alignment",
            "bitscore_improvement": "Relative improvement of fused vs control bit scores",
            "overlap_equilibrium": "Balance of alignment between two gene parts (favors 50-50)",
            "organism_count": "Number of different organisms with hits to this fused gene",
            "evalue": "E-value significance of alignment (lower is better)"
        },
        "score_range": "0-100 (higher = more likely true misannotation)"
    }
