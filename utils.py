"""
utils.py
------------
GEM quality control utilities.
Used in obj1_gem_qc_lplantarum.ipynb and obj1_gem_qc_lmesenteroides.ipynb.
"""

import cobra
import cobra.flux_analysis
import memote
import memote.suite.api as memote_api


def gem_summary(model: cobra.Model, label: str) -> dict:
    """
    Standard QC checks matching the pipeline described in the article:
      1. Memote score
      2. Model stats (reactions, metabolites, genes)
      3. Growth rate in default medium (FBA)
      4. Biomass escapes (growth with all exchanges closed)
      5. Blocked reactions via FVA (fraction_of_optimum=0)

    Parameters
    ----------
    model : cobra.Model
        Loaded COBRApy model. Not modified in place.
    label : str
        Short name for display (e.g. 'iLP728', 'Koduru2022').

    Returns
    -------
    dict with all results. Prints a formatted summary as side effect.
    """
    results = {"label": label}


    # ── 2. Stats ──────────────────────────────────────────────────────────────
    results["n_reactions"]   = len(model.reactions)
    results["n_metabolites"] = len(model.metabolites)
    results["n_genes"]       = len(model.genes)

    # ── 3. Growth — default medium ────────────────────────────────────────────
    sol = model.optimize()
    results["growth_default"] = round(sol.objective_value, 4) if sol.status == "optimal" else None
    results["fba_status"]     = sol.status

    # ── 4. Biomass escapes — all exchanges closed ─────────────────────────────
    with model:
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)
        sol_closed = model.optimize()
        leak = sol_closed.objective_value if sol_closed.status == "optimal" else 0.0
    results["biomass_escape"] = round(leak, 6)

    # ── 5. Blocked reactions via FVA ──────────────────────────────────────────
    fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0)
    blocked = fva[(fva.maximum == 0) & (fva.minimum == 0)]
    results["n_blocked"]   = len(blocked)
    results["pct_blocked"] = round(len(blocked) / len(model.reactions) * 100, 1)

    _print_summary(results)
    return results


def _print_summary(r: dict) -> None:
    w = 30
    print(f"\n{'─' * 48}")
    print(f"  GEM QC — {r['label']}")
    print(f"{'─' * 48}")
    print(f"  {'Reactions':<{w}} {r['n_reactions']}")
    print(f"  {'Metabolites':<{w}} {r['n_metabolites']}")
    print(f"  {'Genes':<{w}} {r['n_genes']}")
    print(f"  {'Growth (default medium)':<{w}} {r['growth_default']} h⁻¹  [{r['fba_status']}]")
    print(f"  {'Biomass escape':<{w}} {r['biomass_escape']}")
    print(f"  {'Blocked reactions':<{w}} {r['n_blocked']} ({r['pct_blocked']}%)")
    print(f"{'─' * 48}\n")