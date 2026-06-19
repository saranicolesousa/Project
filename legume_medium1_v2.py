"""
legume_medium1.py  —  v2
------------------------
Shared module: nutritional composition of chickpea and fava bean,
unit-conversion functions, and exchange-bound mappings for iLP728
and Koduru2022.
"""

# ── System parameters ─────────────────────────────────────────────────────────
X                = 0.5   # gDW/L   (average cell density over fermentation)
T                = 24    # h       (fermentation duration)
RAZAO_HIDRATACAO = 5     # L medium / kg dry flour
HUMIDADE_CHICKPEA = 0.10  # 10% moisture in dry chickpea flour


# ── Molecular weights (g/mol) ─────────────────────────────────────────────────
MW = {
    # Free sugars
    "glicose":      180.16, "frutose":    180.16, "galactose":  180.16,
    "ribose":       150.13, "manose":     180.16,
    # Disaccharides
    "sacarose":     342.30, "maltose":    342.30,
    # Raffinose-family oligosaccharides (RFOs)
    "raffinose":    504.44, "stachyose":  666.58,
    # Amino acids (total protein hydrolysate; see note above)
    "glutamato":    147.13, "aspartato":  133.10, "arginina":   174.20,
    "leucina":      131.17, "lisina":     146.19, "fenilalanina":165.19,
    "serina":        87.08, "alanina":     89.09, "prolina":    115.13,
    "glicina":       75.03, "valina":     117.15, "isoleucina": 131.17,
    "treonina":     119.12, "tirosina":   181.19, "histidina":  155.15,
    "cisteina":     121.16, "metionina":  149.21, "triptofano": 204.23,
    # Free amino acids (reported separately from total protein AAs)
    # These are composition-derived and use efficiency=1.0
    "asparagina":   132.12, "glutamina":   146.15,
    # Vitamins (soluble fraction)
    "riboflavina":  376.36, "tiamina":    337.27, "niacina":    123.11,
    "piridoxina":   169.18, "biotina":    244.31, "folato":     441.40,
    "pantotenato":  219.23, "ascorbato":  176.12,
    # Minerals (ionic; always efficiency=1.0)
    "fosfato":       94.97, "sulfato":     96.06, "manganes":    54.94,
    "sodio":         22.99, "magnesio":    24.31,
    # Fatty acids
    "palmítico":    256.42, "esteárico":  284.48, "oleico":     282.46,
    # LOX-derived off-flavors (free volatile compounds)
    "hexanal":      100.16, "nonenal":    140.22,
}

# ── Efficiency fractions ───────────────────────────────────────────────────────
# fraction of the listed total composition that becomes available for uptake
# within the 24h fermentation window. Justification per category:
#   free_sugar   1.0  — already dissolved in slurry; no hydrolysis step required
#   disaccharide 0.90 — sucrose fully disappeared in L. plantarum-fermented chickpea flour within the fermentation period
#                       (Mansour & Khalil 1995, Menofiya J Agric Res 20:151-164). Same value applied to Lm in the absence of a
#                       species-specific number; both species are known sucrose/maltose fermenters.
#   rfo_lp       0.55 — raffinose reduced 52.9% by L. plantarum in fermented chickpea flour (Mansour & Khalil 1995). 
#                       Convergent literature across legume flours (chickpea, lentil, bean, pea) reports 53-81% raffinose 
#                       reduction by L. plantarum/L. brevis within ~24h (Coda et al. 2015, Eur Food Res Technol; 
#                       review in Sharma & Goyal 2021,
#                       IntechOpen ch. "Fermentation as Strategy..."). 0.55 is a conservative anchor within this range, 
#                       not the upper bound.
#   stachyose_lp 0.00 — structural limitation, not biological: iLP728 has no EX_styse_e / transport reaction for stachyose. 
#                       Note: real L. plantarum strains (incl. the chickpea study above) show complete stachyose elimination,
#                       so this is a documented gap in the GEM reconstruction, not a reflection of WCFS1 biology — flag in Limitations.
#   rfo_lm       0.60 — Leu. mesenteroides is the only species (vs L. lactis, S. thermophilus) carrying the complete RFO gene
#                       cluster (lacS-galA-galK-galT) with high alpha-galactosidase activity (Wang et al. 2021,
#                       LWT/ScienceDirect, "Capacity of soybean carbohydrate metabolism..."). No chickpea/fava-bean-specific
#                       percentage found; kept at parity with the revised Lp value given the qualitatively stronger enzymatic
#                       toolkit.
#   stachyose_lm 0.40 — kept below raffinose; stachyose requires an additional hydrolysis step and is consistently
#                       reported as harder to degrade than raffinose across LAB/legume fermentation studies. No Lm-specific
#                       number found; this is an estimate.
#   aa           0.20 — SHARED value for Lp and Lm (no organism-specific differentiation — see rationale below). Chickpea-
#                       specific: total essential AA increased 6.3-11.3% during L. plantarum fermentation (Mansour & Khalil
#                       1995). Fava bean: free AA content increased qualitatively after L. plantarum fermentation
#                       (Rizzello et al. 2019, Food Res Int, faba bean LAB fermentation), no extractable percentage. These
#                       metrics measure net compositional change, not the substrate-availability fraction this parameter
#                       represents, so 0.20 is an estimate anchored in the general magnitude reported, not a literal derivation.
#                       Tested as sensitivity scenario between 0.12-0.25 (see obj2, Section 6).
#   vitamin      1.0  — small soluble molecules
#   mineral      1.0  — ionic species
#   lipid        0.12 — LAB do not beta-oxidise fatty acids for energy, but readily and substantially incorporate exogenous
#                       unsaturated fatty acids (oleic, linoleic) into membrane phospholipids once liberated as free fatty
#                       acid (Enterococcus/Lactobacillus studies, e.g. Fabre et al. 2014 AEM/PMC4178632 — >2/3 of membrane
#                       content became linoleic acid upon supplementation).
#                       The real bottleneck is upstream: most legume crude fat is triglyceride-bound, and LAB generally lack
#                       extracellular lipase. Legume lipoxygenase activity requires free linoleic acid and is well documented
#                       as active upon tissue disruption (Shi et al. 2020,Legume Science), implying a non-trivial free-FA pool
#                       exists from slurry preparation onward. 0.12 raised from 0.05 to avoid artificially suppressing the
#                       linoleic acid consumption mechanism central to the off-flavor mitigation narrative — treat as an
#                       estimate, test in sensitivity given its importance to that conclusion.
#   off_flavor   1.0  — free volatile compounds already present

_AAs = [
    "glutamato", "aspartato", "arginina", "leucina", "lisina",
    "fenilalanina", "serina", "alanina", "prolina", "glicina",
    "valina", "isoleucina", "treonina", "tirosina", "histidina",
    "cisteina", "metionina", "triptofano",
]
_VITAMINS = [
    "riboflavina", "tiamina", "niacina", "piridoxina",
    "biotina", "folato", "pantotenato", "ascorbato",
]
_MINERALS = ["fosfato", "sulfato", "manganes", "sodio", "magnesio"]
_LIPIDS   = ["palmítico", "esteárico", "oleico"]
_FREE_SUG = ["glicose", "frutose", "galactose", "ribose", "manose"]
_DISAC    = ["sacarose", "maltose"]

EFICIENCIA_LP = {
    **{k: 1.00 for k in _FREE_SUG},
    **{k: 0.90 for k in _DISAC},
    "raffinose":  0.55,  # Mansour & Khalil 1995 — 52.9% reduction, chickpea-specific
    "stachyose":  0.00,  # no EX_styse_e in iLP728 — model cannot consume it
    **{k: 0.20 for k in _AAs},
    **{k: 1.00 for k in _VITAMINS},
    **{k: 1.00 for k in _MINERALS},
    **{k: 0.12 for k in _LIPIDS},
    "hexanal":    1.00,
    "nonenal":    1.00,
    # Free amino acids — already dissolved, efficiency=1.0
    "asparagina": 1.00,
    "glutamina":  1.00,
}

EFICIENCIA_LM = {
    **{k: 1.00 for k in _FREE_SUG},
    **{k: 0.90 for k in _DISAC},
    "raffinose":  0.60,  # Koduru has EX_raffin_e; Lm is saccharolytic
    "stachyose":  0.40,  # Koduru has EX_styse_e; lower efficiency than raffinose
    **{k: 0.20 for k in _AAs},
    **{k: 1.00 for k in _VITAMINS},
    **{k: 1.00 for k in _MINERALS},
    **{k: 0.12 for k in _LIPIDS},
    "hexanal":    1.00,
    "nonenal":    1.00,
    # Free amino acids — already dissolved, efficiency=1.0
    "asparagina": 1.00,
    "glutamina":  1.00,
}

# ── Conversion functions ──────────────────────────────────────────────────────
def calcular_bound(C_mgkg_dw, MW_gmol, efficiency=1.0,
                   X=X, t=T, razao=RAZAO_HIDRATACAO):
    """
    Convert mg/kg DW → uptake bound (mmol/gDW/h, negative).

    Implements the Sauer batch-culture specific consumption rate:
        qS = ∆S · (μ / X)

    Approximated here as: qS ≈ S0 / (X * t)  (assumes near-complete
    consumption or a conservative upper-bound on availability).

    Parameters
    ----------
    C_mgkg_dw  : float   nutrient concentration in mg/kg DW flour
    MW_gmol    : float   molecular weight in g/mol
    efficiency : float   fraction of C_mgkg_dw that becomes bioavailable
                         (0.0–1.0; see EFICIENCIA_LP / EFICIENCIA_LM)
    """
    if C_mgkg_dw <= 0 or efficiency <= 0:
        return 0.0
    C_eff  = C_mgkg_dw * efficiency      # mg/kg DW after accessibility correction
    C_meio = C_eff / razao               # mg/L (hydration dilution)
    C_gL   = C_meio / 1000              # g/L
    v      = C_gL / (MW_gmol * X * t)   # mol/L/h/gDW
    return -v * 1000                     # mmol/gDW/h (negative = uptake)


# ── Unit helpers for composition data ────────────────────────────────────────
def fw_to_dw_mgkg(C_g_per_100g_fw, humidade=HUMIDADE_CHICKPEA):
    """g/100g FW → mg/kg DW"""
    return C_g_per_100g_fw / (1 - humidade) * 10000

def mgper100g_to_mgperkg(C):
    """mg/100g → mg/kg (×10)"""
    return C * 10

def cp_fw(low, high):
    """Chickpea: g/100g FW → mg/kg DW  →  (min, med, max)"""
    mn = fw_to_dw_mgkg(low);  mx = fw_to_dw_mgkg(high)
    return mn, (mn + mx) / 2, mx

def cp_mg(low, high):
    """Chickpea: mg/100g FW → mg/kg DW  →  (min, med, max)"""
    f  = 1 / (1 - HUMIDADE_CHICKPEA)
    mn = mgper100g_to_mgperkg(low)  * f
    mx = mgper100g_to_mgperkg(high) * f
    return mn, (mn + mx) / 2, mx

def cp_ug(low, high):
    """Chickpea: µg/100g FW → mg/kg DW  →  (min, med, max)"""
    f  = 1 / (1 - HUMIDADE_CHICKPEA)
    mn = (low  * 1e-3 * 10) * f
    mx = (high * 1e-3 * 10) * f
    return mn, (mn + mx) / 2, mx

def fb_g(low, high):
    """Fava bean: g/100g DW → mg/kg DW  →  (min, med, max)"""
    mn = low * 10000; mx = high * 10000
    return mn, (mn + mx) / 2, mx

def fb_mg(low, high):
    """Fava bean: mg/100g DW → mg/kg DW  →  (min, med, max)"""
    mn = mgper100g_to_mgperkg(low); mx = mgper100g_to_mgperkg(high)
    return mn, (mn + mx) / 2, mx

def fb_ug(low, high):
    """Fava bean: µg/100g DW → mg/kg DW  →  (min, med, max)"""
    mn = low * 1e-3 * 10; mx = high * 1e-3 * 10
    return mn, (mn + mx) / 2, mx

def fb_mg_directo(low, high):
    """Fava bean: mg/100g DW (sugars) → mg/kg DW  →  (min, med, max)"""
    mn = low * 10; mx = high * 10
    return mn, (mn + mx) / 2, mx

def offlavor_dw(low, high):
    """
    LOX-derived off-flavors already in mg/kg DW (= µg/g DW).
    Consistent with volatile aldehyde profiles in raw legume flour
    (range ~0.1–5 mg/kg DW for hexanal; Verni 2022).
    Returns (min, med, max) in mg/kg DW — same units as all other compounds,
    so calcular_bound() handles them identically.
    """
    mn = low; mx = high
    return mn, (mn + mx) / 2, mx


# ── Nutritional composition ───────────────────────────────────────────────────
# Sources: FRIDA (frida.fooddata.dk), USDA FoodData, Tao 2022, Verni 2022,
#          Margier et al. (2018) Int J Mol Sci, Carbonaro et al. (2012)
#          Br J Nutr. Amino acid values are from total protein hydrolysate
#          tables (NOT free amino acids). See efficiency section above.

compostos_cp = {
    # Free sugars (g/100g FW → mg/kg DW via cp_fw)
    "glicose":       cp_fw(0.5,  1.0),
    "frutose":       cp_fw(0.0,  0.2),
    "galactose":     cp_fw(0.0,  1.0),
    "ribose":        cp_fw(0.0,  0.15),
    "manose":        cp_fw(0.0,  0.08),
    # Disaccharides
    "sacarose":      cp_fw(1.5,  3.0),
    "maltose":       cp_fw(0.0,  0.9),
    # RFOs — chickpea has trace raffinose; stachyose absent
    # (chickpea is not a major RFO source; fava bean is)
    # Amino acids — total protein hydrolysate (mg/100g FW → mg/kg DW)
    "glutamato":     cp_mg(3000, 4200),
    "aspartato":     cp_mg(2200, 3200),
    "arginina":      cp_mg(1700, 2500),
    "leucina":       cp_mg(1300, 2000),
    "lisina":        cp_mg(1200, 1800),
    "fenilalanina":  cp_mg(1000, 1600),
    "serina":        cp_mg(900,  1400),
    "alanina":       cp_mg(750,  1200),
    "prolina":       cp_mg(700,  1100),
    "glicina":       cp_mg(700,  1100),
    "valina":        cp_mg(800,  1200),
    "isoleucina":    cp_mg(700,  1100),
    "treonina":      cp_mg(650,  1000),
    "tirosina":      cp_mg(500,   900),
    "histidina":     cp_mg(450,   750),
    "cisteina":      cp_mg(250,   500),
    "metionina":     cp_mg(200,   450),
    "triptofano":    cp_mg(180,   350),
    # Vitamins (mg/100g FW → mg/kg DW)
    "riboflavina":   cp_mg(0.08,  0.21),
    "tiamina":       cp_mg(0.31,  0.48),
    "niacina":       cp_mg(1.5,   2.5),
    "piridoxina":    cp_mg(0.44,  0.65),
    "pantotenato":   cp_mg(0.9,   1.6),
    "biotina":       cp_ug(14,    25),
    "folato":        cp_ug(181,   557),
    "ascorbato":     cp_mg(0.0,   4.0),  # absent in both models; silently skipped
    # Minerals (mg/100g FW → mg/kg DW)
    "fosfato":       cp_mg(250,   366),
    "manganes":      cp_mg(0.45,  2.6),
    "sodio":         cp_mg(2.5,   10),
    "magnesio":      cp_mg(115,   164),  # no EX_mg2_e in iLP728 (implicit cofactor)
    # Fatty acids (g/100g FW → mg/kg DW)
    "palmítico":     cp_fw(0.5,   0.8),
    "esteárico":     cp_fw(0.1,   0.3),
    "oleico":        cp_fw(1.2,   2.1),
    # Free amino acids (separate from total protein pool above)
    # Source: total free AA in chickpea ~129 mg/100g DW (Celleno et al. 2021,
    # PMC8002183). Asn ~25-35% of free AA pool; Gln ~15-25%.
    # These are free (not protein-bound), so efficiency=1.0 in both strains.
    # Values in mg/100g FW (converted via cp_mg → mg/kg DW).
    "asparagina":    cp_mg(22, 45),    # ~25-50 mg/100g DW free Asn
    "glutamina":     cp_mg(16, 36),    # ~18-40 mg/100g DW free Gln
    # LOX-derived off-flavors (mg/kg DW — already final units)
    "hexanal":       offlavor_dw(0.5,  1.0),
    "nonenal":       offlavor_dw(0.025, 0.045),
}

compostos_fb = {
    # Free sugars (mg/100g DW → mg/kg DW via fb_mg_directo)
    "glicose":       fb_mg_directo(10,   80),
    "frutose":       fb_mg_directo(10,   60),
    "sacarose":      fb_g(0.8,  2.2),
    # RFOs — major in fava bean
    "raffinose":     fb_mg(110, 390),
    "stachyose":     fb_mg(440, 1370),
    # Amino acids — total protein hydrolysate (mg/100g DW → mg/kg DW)
    "glutamato":     fb_mg(3800, 5300),
    "aspartato":     fb_mg(2600, 3900),
    "arginina":      fb_mg(2200, 2700),
    "leucina":       fb_mg(1321, 2379),
    "lisina":        fb_mg(1165, 2468),
    "fenilalanina":  fb_mg(850,  1250),
    "serina":        fb_mg(700,  1400),
    "alanina":       fb_mg(370,  1410),
    "prolina":       fb_mg(600,  1200),
    "glicina":       fb_mg(500,  1200),
    "valina":        fb_mg(900,  1350),
    "isoleucina":    fb_mg(850,  1200),
    "treonina":      fb_mg(750,  1050),
    "tirosina":      fb_mg(500,  1000),
    "histidina":     fb_mg(430,   780),
    "cisteina":      fb_mg(100,   300),
    "metionina":     fb_mg(130,   330),
    "triptofano":    fb_mg(180,   310),
    # Vitamins
    "riboflavina":   fb_mg(0.2,   0.4),
    "tiamina":       fb_mg(0.4,   0.7),
    "niacina":       fb_mg(2.0,   3.5),
    "piridoxina":    fb_mg(0.2,   0.5),
    "pantotenato":   fb_mg(0.2,   0.6),
    "biotina":       fb_ug(1,     5),
    "folato":        fb_ug(320,   520),
    "ascorbato":     fb_mg(0.5,   2.0),
    # Minerals
    "fosfato":       fb_mg(421,   755),
    "manganes":      fb_mg(0.9,   2.7),
    "sodio":         fb_mg(10,    20),
    "magnesio":      fb_mg(130,   250),
    # Fatty acids
    "palmítico":     fb_mg(136,   374),
    "esteárico":     fb_mg(30,    140),
    "oleico":        fb_mg(150,   560),
    # Free amino acids (separate from total protein pool)
    # Fava bean accumulates more free Asn than chickpea (amide-exporter
    # from N2 fixation; Thavarajah et al. 2005). Range from literature:
    # Asn 200-400 mg/100g DW; Gln 50-150 mg/100g DW.
    "asparagina":    fb_mg(200, 400),
    "glutamina":     fb_mg(50,  150),
    # LOX-derived off-flavors (mg/kg DW)
    "hexanal":       offlavor_dw(1.0,  2.1),
    "nonenal":       offlavor_dw(0.05, 0.09),
}

# Scenario shortcuts
chickpea_min = {c: v[0] for c, v in compostos_cp.items()}
chickpea_med = {c: v[1] for c, v in compostos_cp.items()}
chickpea_max = {c: v[2] for c, v in compostos_cp.items()}
favabean_min = {c: v[0] for c, v in compostos_fb.items()}
favabean_med = {c: v[1] for c, v in compostos_fb.items()}
favabean_max = {c: v[2] for c, v in compostos_fb.items()}


# ── Exchange ID mappings ──────────────────────────────────────────────────────

# iLP728_curated_v2.xml — BiGG convention
# CORRECTIONS vs v1:
#   galactose  EX_galt_e → EX_gal_e  (EX_galt_e is galactitol, not galactose)
#   magnesio   removed   (EX_mg2_e does not exist; Mg implicit cofactor in iLP728)
#   raffinose  added     (EX_raffin_e + RAFFINt added to iLP728_curated_v2.xml)
#   stachyose  present but efficiency=0.0 (no exchange; silently returns 0 bound)
mapa_ilp728 = {
    "glicose":      "EX_glc__D_e",
    "frutose":      "EX_fru_e",
    "galactose":    "EX_gal_e",       # corrected (was EX_galt_e = galactitol)
    "ribose":       "EX_rib__D_e",
    "manose":       "EX_man_e",
    "sacarose":     "EX_sucr_e",
    "maltose":      "EX_malt_e",
    "raffinose":    "EX_raffin_e",    # added to model in iLP728_curated_v2.xml
    "stachyose":    None,  # no exchange in iLP728; eff=0.00 in EFICIENCIA_LP.
                       # Stachyose is bioavailable only to Lm (EX_styse_e in Koduru, eff=0.40). This is the primary RFO advantage of Lm in fava bean medium.
    "glutamato":    "EX_glu__L_e",
    "aspartato":    "EX_asp__L_e",
    "arginina":     "EX_arg__L_e",
    "leucina":      "EX_leu__L_e",
    "lisina":       "EX_lys__L_e",
    "fenilalanina": "EX_phe__L_e",
    "serina":       "EX_ser__L_e",
    "alanina":      "EX_ala__L_e",
    "prolina":      "EX_pro__L_e",
    "glicina":      "EX_gly_e",
    "valina":       "EX_val__L_e",
    "isoleucina":   "EX_ile__L_e",
    "treonina":     "EX_thr__L_e",
    "tirosina":     "EX_tyr__L_e",
    "histidina":    "EX_his__L_e",
    "cisteina":     "EX_cys__L_e",
    "metionina":    "EX_met__L_e",
    "triptofano":   "EX_trp__L_e",
    "riboflavina":  "EX_ribflv_e",
    "tiamina":      "EX_thm_e",
    "niacina":      "EX_nac_e",
    "piridoxina":   "EX_pydxn_e",
    "biotina":      "EX_btn_e",
    "folato":       "EX_fol_e",
    "pantotenato":  "EX_pnto__R_e",
    "ascorbato":    None,             # EX_ascb_e absent in iLP728
    "fosfato":      "EX_pi_e",
    "manganes":     "EX_mn2_e",
    "sodio":        "EX_na1_e",
    "magnesio":     None,             # EX_mg2_e absent; Mg modelled as implicit
    "palmítico":    "EX_hdca_e",
    "esteárico":    "EX_ocdca_e",
    "oleico":       "EX_ocdcea_e",
    "hexanal":      "EX_hxa_e",
    "nonenal":      "EX_nnl_e",
    # Free amino acids (composition-derived; overrides adicional bound)
    "asparagina":   "EX_asn__L_e",
    "glutamina":    "EX_gln__L_e",
}

# Koduru2022_curated_v2.xml
# Glucose ID harmonized to BiGG convention (EX_glc__D_e) in v2 of the Koduru model.
# glc-D_e → glc__D_e in metabolite IDs; EX_glc_e → EX_glc__D_e in exchange ID.
# Affected reactions: EX_glc__D_e, GLCpts, GLCt2, TREHe, MALTe.
# RFO exchanges (EX_raffin_e, EX_styse_e, EX_gal_e) are at (0,1000) by default
# and are correctly opened by aplicar_meio_leguminosa() via the composition bounds.
mapa_koduru = {
    "glicose":      "EX_glc__D_e",   # harmonized to BiGG convention in Koduru v2
    "frutose":      "EX_fru_e",
    "galactose":    "EX_gal_e",      # EX_gal_e exists in Koduru; was (0,1000) default
    "ribose":       "EX_rib__D_e",
    "manose":       "EX_man_e",
    "sacarose":     "EX_sucr_e",
    "maltose":      "EX_malt_e",
    "raffinose":    "EX_raffin_e",   # exists in Koduru; was (0,1000) default
    "stachyose":    "EX_styse_e",    # exists in Koduru; was (0,1000) default
    "glutamato":    "EX_glu__L_e",
    "aspartato":    "EX_asp__L_e",
    "arginina":     "EX_arg__L_e",
    "leucina":      "EX_leu__L_e",
    "lisina":       "EX_lys__L_e",
    "fenilalanina": "EX_phe__L_e",
    "serina":       "EX_ser__L_e",
    "alanina":      "EX_ala__L_e",
    "prolina":      "EX_pro__L_e",
    "glicina":      "EX_gly_e",
    "valina":       "EX_val__L_e",
    "isoleucina":   "EX_ile__L_e",
    "treonina":     "EX_thr__L_e",
    "tirosina":     "EX_tyr__L_e",
    "histidina":    "EX_his__L_e",
    "cisteina":     "EX_cys__L_e",
    "metionina":    "EX_met__L_e",
    "triptofano":   "EX_trp__L_e",
    "riboflavina":  "EX_ribflv_e",
    "tiamina":      "EX_thm_e",
    "niacina":      "EX_nac_e",
    "piridoxina":   "EX_pydxn_e",
    "biotina":      "EX_btn_e",
    "folato":       "EX_fol_e",
    "pantotenato":  "EX_pnto__R_e",
    "ascorbato":    None,             # EX_ascb_e absent in Koduru
    "fosfato":      "EX_pi_e",
    "manganes":     "EX_mn2_e",
    "sodio":        "EX_na1_e",
    "magnesio":     "EX_mg_e",       # exists in Koduru (EX_mg_e, not EX_mg2_e)
    "palmítico":    "EX_hdca_e",
    "esteárico":    "EX_ocdca_e",
    "oleico":       "EX_ocdcea_e",
    "hexanal":      "EX_hxa_e",
    "nonenal":      "EX_nnl_e",
    # Free amino acids (composition-derived; overrides adicional bound)
    "asparagina":   "EX_asn__L_e",
    "glutamina":    "EX_gln__L_e",
}
def aplicar_meio_leguminosa(model, composicao_mgkg_dw, mapa_exchanges,
                             nome_matriz="leguminosa", modelo_tipo="ilp728",
                             eficiencia=None, verbose=True):
    """
    Apply legume fermentation medium to a COBRApy model.

    Parameters
    ----------
    model             : cobra.Model (modified in-place; use `with model:` for context)
    composicao_mgkg_dw: dict  {compound_name: concentration_mg_per_kg_DW}
    mapa_exchanges    : dict  {compound_name: exchange_reaction_id | None}
    modelo_tipo       : "ilp728" or "koduru"
    eficiencia        : dict  {compound_name: float}; defaults to strain-specific map
    verbose           : bool
    """
    if eficiencia is None:
        eficiencia = EFICIENCIA_LP if modelo_tipo == "ilp728" else EFICIENCIA_LM

    # Basal inorganic supplements (not concentration-derived)
    # These are always available in excess and are not legume-specific.
    if modelo_tipo == "ilp728":
        universais = {
            "EX_h2o_e": (-1000, 1000),
            "EX_h_e":   (-1000, 1000),
            "EX_nh4_e": (-10,   1000),
            "EX_so4_e": (-10,   1000),
            "EX_pi_e":  (-10,   1000),
            "EX_mn2_e": (-10,   1000),
            "EX_na1_e": (-10,   1000),
            "EX_co2_e": (0,     1000),
            "EX_o2_e":  (-1,1000)
        }
        # Nucleobases and amino acid supplements not in legume composition
        # but required for growth (auxotrophic requirements of WCFS1).
        # Note: EX_asn__L_e and EX_gln__L_e are removed from adicionais —
        # they are now composition-derived (free AA values in compostos_cp/fb).
        adicionais = [
            "EX_ade_e", "EX_gua_e", "EX_ura_e", "EX_hxan_e",
            "EX_ins_e", "EX_thymd_e", "EX_pydam_e",
        ]
    elif modelo_tipo == "koduru":
        universais = {
            "EX_h2o_e": (-1000, 1000),
            "EX_h_e":   (-1000, 1000),
            "EX_nh4_e": (-10,   1000),
            "EX_so4_e": (-10,   1000),
            "EX_pi_e":  (-10,   1000),
            "EX_mn2_e": (-10,   1000),
            "EX_na1_e": (-10,   1000),
            "EX_co2_e": (0,     1000),
            "EX_mg_e":  (-10,   1000),
            "EX_fe2_e": (-10,   1000),
            "EX_k_e":   (-10,   1000),
            "EX_o2_e":  (-1,1000)
        }
        # Note: EX_asn__L_e and EX_gln__L_e removed from adicionais —
        # now composition-derived (free AA values in compostos_cp/fb).
        adicionais = [
            "EX_ade_e", "EX_gua_e", "EX_ura_e",
            "EX_thymd_e",
        ]
    else:
        raise ValueError(f"Unknown modelo_tipo: {modelo_tipo!r}. Use 'ilp728', 'koduru', or 'community'.")
    

    # Step 1: close all exchanges (allow secretion only)
    for rxn in model.exchanges:
        rxn.bounds = (0, 1000)

    # Step 2: open basal inorganics
    for ex_id, bounds in universais.items():
        try:
            model.reactions.get_by_id(ex_id).bounds = bounds
        except KeyError:
            pass

    # Step 3: open auxotrophic supplements
    for ex_id in adicionais:
        try:
            model.reactions.get_by_id(ex_id).bounds = (-1, 1000)
        except KeyError:
            pass

    # Step 4: apply composition-derived bounds with efficiency correction
    bounds_aplicados = {}
    nao_encontrados  = []
    skipped_none     = []

    for composto, C_mgkg in composicao_mgkg_dw.items():
        ex_id = mapa_exchanges.get(composto)
        if ex_id is None:
            skipped_none.append(composto)
            continue
        if composto not in MW:
            continue

        eff   = eficiencia.get(composto, 1.0)
        bound = calcular_bound(C_mgkg, MW[composto], efficiency=eff)

        if bound == 0.0:
            continue  # zero bound: leave exchange closed

        try:
            model.reactions.get_by_id(ex_id).bounds = (bound, 1000)
            bounds_aplicados[ex_id] = bound
        except KeyError:
            nao_encontrados.append(f"{composto} ({ex_id})")

    if verbose:
        print(f"\n=== Meio {nome_matriz} ({modelo_tipo}) ===")
        print(f"  Bounds aplicados: {len(bounds_aplicados)}")
        if nao_encontrados:
            print(f"  Não encontrados no modelo: {nao_encontrados}")
        if skipped_none:
            print(f"  Sem exchange mapeado (None): {skipped_none}")

    return bounds_aplicados

def build_smetana_environment(community, composicao, eficiencia, verbose=True):
    """
    Build a SMETANA Environment for a legume fermentation medium.
 
    Works directly with pool reactions in the merged model — no manual
    ID overrides needed. Pool reactions are discovered from the merged
    model itself, handling both R_EX_M_*_pool and R_EX_*_pool naming.
 
    Parameters
    ----------
    community   : smetana.interface.Community
    composicao  : dict  {compound_name: concentration_mg_per_kg_DW}
                  Use chickpea_med, fava_med, etc. from this module.
    eficiencia  : dict  {compound_name: float}
                  Use EFICIENCIA_LP or EFICIENCIA_LM.
    verbose     : bool
 
    Returns
    -------
    env         : smetana.interface.Environment
    log         : list of dicts  (one entry per bound applied)
 
    Notes
    -----
    Magnésio (Lm): EX_mg_e no Lm aponta para M_mg2_e no compartimento
    externo partilhado. O pool é R_EX_M_mg2_e_pool. O mapa_koduru mantém
    'EX_mg_e' (exchange real do organismo); aqui faz-se o override explícito
    para o pool correcto.
    """
    from smetana.interface import Environment
    merged = community.merged
 
    # ------------------------------------------------------------------
    # 1. Construir pool_map robusto
    #    Suporta R_EX_M_*_pool  (metabolitos com prefixo M_)
    #    e       R_EX_*_pool    (metabolitos sem prefixo M_)
    # ------------------------------------------------------------------
    pool_map = {}
    for r_id in merged.reactions:
        if not r_id.endswith('_pool'):
            continue
        if r_id.startswith('R_EX_M_'):
            met_id = r_id[7:-5]      # remove 'R_EX_M_' e '_pool'
        elif r_id.startswith('R_EX_'):
            met_id = r_id[5:-5]      # remove 'R_EX_' e '_pool'
        else:
            continue
        pool_map[met_id] = r_id
 
    # Override explícito: EX_mg_e (Lm) → pool M_mg2_e
    # O exchange do organismo usa mg_e mas o metabolito externo partilhado
    # é M_mg2_e — o mesmo ião, IDs diferentes entre os dois modelos.
    if 'mg2_e' in pool_map and 'mg_e' not in pool_map:
        pool_map['mg_e'] = pool_map['mg2_e']
 
    # ------------------------------------------------------------------
    # 2. Universais — sempre disponíveis, não derivados da composição
    # ------------------------------------------------------------------
    env = Environment()
 
    UNIVERSAIS = {
        'h2o_e' : -1000.0,
        'h_e'   : -1000.0,
        'nh4_e' :   -10.0,
        'pi_e'  :   -10.0,
        'so4_e' :   -10.0,
        'k_e'   :   -10.0,
        'na1_e' :   -10.0,
        'fe2_e' :    -5.0,
        'mn2_e' :    -5.0,
        'mg2_e' :    -5.0,   # cobre Lp (via mg2_e) e Lm (via override mg_e → mg2_e)
        'co2_e' :     0.0,
        'o2_e'  :    -1.0,
    }
    for met_id, lb in UNIVERSAIS.items():
        pool_id = pool_map.get(met_id)
        if pool_id:
            env[pool_id] = (lb, 1000.0)
 
    # ------------------------------------------------------------------
    # 3. Auxotrofias LAB — traço, não derivadas da composição
    #    Nucleobases + vitaminas que as LAB não sintetizam
    # ------------------------------------------------------------------
    AUXOTROFIAS = [
        'ade_e', 'gua_e', 'ura_e', 'thymd_e',
        'hxan_e', 'ins_e', 'pydam_e',
    ]
    for met_id in AUXOTROFIAS:
        pool_id = pool_map.get(met_id)
        if pool_id:
            env[pool_id] = (-1.0, 1000.0)
 
    # ------------------------------------------------------------------
    # 4. Nutrientes quantificados da composição leguminosa
    # ------------------------------------------------------------------
    # Construir mapa composto → met_id a partir dos dois mapas
    # (usa a entrada do mapa_ilp728 e mapa_koduru — os met_ids BiGG
    #  são os mesmos após harmonização; o Lp serve de referência)
    mapa_met = {}
    for mapa in [mapa_ilp728, mapa_koduru]:
        for composto, ex_id in mapa.items():
            if ex_id is None or composto in mapa_met:
                continue
            met_id = ex_id.replace('EX_', '').replace('R_EX_', '')
            mapa_met[composto] = met_id
 
    log = []
    sem_pool = []
    sem_mw   = []
 
    for composto, C_mgkg in composicao.items():
        if composto not in mapa_met:
            continue
        if composto not in MW:
            sem_mw.append(composto)
            continue
 
        met_id  = mapa_met[composto]
        pool_id = pool_map.get(met_id)
 
        if pool_id is None:
            sem_pool.append((composto, met_id))
            continue
 
        eff   = eficiencia.get(composto, 1.0)
        bound = calcular_bound(C_mgkg, MW[composto], efficiency=eff)
 
        if bound == 0.0:
            continue
 
        env[pool_id] = (bound, 1000.0)
        log.append({
            'composto' : composto,
            'met_id'   : met_id,
            'pool_id'  : pool_id,
            'eff'      : eff,
            'lb'       : round(bound, 6),
        })
 
    if verbose:
        print(f"  Bounds aplicados : {len(log)}")
        if sem_pool:
            print(f"  Sem pool reaction: {[c for c, _ in sem_pool]}")
        if sem_mw:
            print(f"  Sem MW definido  : {sem_mw}")
 
    return env, log
 

def build_smetana_environment_(community, composicao, eficiencia, verbose=True):
    """
    Build a SMETANA Environment for a legume fermentation medium.
 
    Works directly with pool reactions in the merged model. Handles two
    known naming issues in the Koduru model:
 
    1. __45__ encoding: Koduru uses hyphens in some metabolite IDs
       (e.g. gln-L_e), which the SMETANA merger encodes as __45__
       (e.g. gln__45__L_e). This creates separate pool reactions for
       the same metabolite. Both pools are opened with the same bound
       so that Lp (uses gln__L_e) and Lm (uses gln__45__L_e) can both
       access the nutrient.
 
    2. mg_e / mg2_e: Koduru exchange is EX_mg_e but the shared external
       metabolite is M_mg2_e. The pool is R_EX_M_mg2_e_pool. Handled
       via explicit override in pool_map.
 
    Parameters
    ----------
    community   : smetana.interface.Community
    composicao  : dict  {compound_name: concentration_mg_per_kg_DW}
                  Use chickpea_med, favabean_med, etc. from this module.
    eficiencia  : dict  {compound_name: float}
                  Use EFICIENCIA_LP or EFICIENCIA_LM.
    verbose     : bool
 
    Returns
    -------
    env         : smetana.interface.Environment
    log         : list of dicts  (one entry per bound applied)
    """
    from smetana.interface import Environment
 
    merged = community.merged
 
    # ------------------------------------------------------------------
    # 1. Construir pool_map robusto
    #    Suporta R_EX_M_*_pool  (metabolitos com prefixo M_)
    #    e       R_EX_*_pool    (metabolitos sem prefixo M_)
    # ------------------------------------------------------------------
    pool_map = {}
    for r_id in merged.reactions:
        if not r_id.endswith('_pool'):
            continue
        if r_id.startswith('R_EX_M_'):
            met_id = r_id[7:-5]
        elif r_id.startswith('R_EX_'):
            met_id = r_id[5:-5]
        else:
            continue
        pool_map[met_id] = r_id
 
    # Override: EX_mg_e (Lm) → pool mg2_e (metabolito externo partilhado)
    if 'mg2_e' in pool_map and 'mg_e' not in pool_map:
        pool_map['mg_e'] = pool_map['mg2_e']
 
    # ------------------------------------------------------------------
    # 2. Construir mapa __45__ → pool sem __45__
    #    Para cada pool com __45__, registar também a chave sem encoding.
    #    Assim gln__L_e → R_EX_M_gln__45__L_e_pool quando não existe
    #    R_EX_M_gln__L_e_pool (e vice-versa).
    # ------------------------------------------------------------------
    for met_id, pool_id in list(pool_map.items()):
        if '__45__' in met_id:
            # gln__45__L_e  →  gln__L_e
            plain = met_id.replace('__45__', '')
            if plain not in pool_map:
                pool_map[plain] = pool_id
        else:
            # gln__L_e  →  gln__45__L_e (se existir pool __45__)
            for sfx in ['__L_e', '__D_e', '__R_e', '__RR_e']:
                if met_id.endswith(sfx):
                    encoded = met_id[:-len(sfx)] + '__45__' + sfx[2:]
                    if encoded in pool_map and met_id not in pool_map:
                        pool_map[met_id] = pool_map[encoded]
                    break
 
    # ------------------------------------------------------------------
    # 3. Universais — sempre disponíveis, não derivados da composição
    # ------------------------------------------------------------------
    env = Environment()
 
    UNIVERSAIS = {
        'h2o_e' : -1000.0,
        'h_e'   : -1000.0,
        'nh4_e' :   -10.0,
        'pi_e'  :   -10.0,
        'so4_e' :   -10.0,
        'k_e'   :   -10.0,
        'na1_e' :   -10.0,
        'cl_e'  :   -10.0,
        'fe2_e' :    -5.0,
        'mn2_e' :    -5.0,
        'mg2_e' :    -5.0,
        'zn2_e' :    -1.0,
        'co2_e' :     0.0,
        'o2_e'  :    -1.0,
    }
    for met_id, lb in UNIVERSAIS.items():
        for key in [met_id, met_id.replace('__', '__45__')]:
            pool_id = pool_map.get(key)
            if pool_id:
                env[pool_id] = (lb, 1000.0)
                break
 
    # ------------------------------------------------------------------
    # 4. Auxotrofias LAB — traço, não derivadas da composição
    # ------------------------------------------------------------------
    AUXOTROFIAS = [
        'ade_e', 'gua_e', 'ura_e', 'thymd_e',
        'hxan_e', 'ins_e', 'pydam_e', 'nmn_e',
    ]
    for met_id in AUXOTROFIAS:
        pool_id = pool_map.get(met_id)
        if pool_id:
            env[pool_id] = (-1.0, 1000.0)
 
    # ------------------------------------------------------------------
    # 5. Nutrientes quantificados da composição leguminosa
    #    Para cada composto, abre a pool sem __45__ E a pool com __45__
    #    (se existir) com o mesmo bound — garante acesso a Lp e Lm.
    # ------------------------------------------------------------------
    mapa_met = {}
    for mapa in [mapa_ilp728, mapa_koduru]:
        for composto, ex_id in mapa.items():
            if ex_id is None or composto in mapa_met:
                continue
            met_id = ex_id.replace('EX_', '').replace('R_EX_', '')
            mapa_met[composto] = met_id
 
    log = []
    sem_pool = []
    sem_mw   = []
 
    for composto, C_mgkg in composicao.items():
        if composto not in mapa_met:
            continue
        if composto not in MW:
            sem_mw.append(composto)
            continue
 
        met_id = mapa_met[composto]
        eff    = eficiencia.get(composto, 1.0)
        bound  = calcular_bound(C_mgkg, MW[composto], efficiency=eff)
 
        if bound == 0.0:
            continue
 
        # Tentar pool sem __45__ e pool com __45__
        pools_abertas = []
        for sfx in ['__L_e', '__D_e', '__R_e']:
            if met_id.endswith(sfx):
                encoded = met_id[:-len(sfx)] + '__45__' + sfx[2:]
                for candidate in [met_id, encoded]:
                    pid = pool_map.get(candidate)
                    if pid and pid not in pools_abertas:
                        pools_abertas.append(pid)
                break
        else:
            pid = pool_map.get(met_id)
            if pid:
                pools_abertas.append(pid)
 
        if not pools_abertas:
            sem_pool.append((composto, met_id))
            continue
 
        for pid in pools_abertas:
            env[pid] = (bound, 1000.0)
 
        log.append({
            'composto'  : composto,
            'met_id'    : met_id,
            'pools'     : pools_abertas,
            'n_pools'   : len(pools_abertas),
            'eff'       : eff,
            'lb'        : round(bound, 6),
        })
 
    if verbose:
        duplos = sum(1 for e in log if e['n_pools'] > 1)
        print(f"  Bounds aplicados : {len(log)} compostos "
              f"({duplos} com pool dupla __45__)")
        if sem_pool:
            print(f"  Sem pool reaction: {[c for c, _ in sem_pool]}")
        if sem_mw:
            print(f"  Sem MW definido  : {sem_mw}")
 
    return env, log

def harmonize_koduru_ids(lm_model):
    """
    Rename __45__ encoded IDs in the Koduru model to BiGG standard before
    building the SMETANA community.
 
    The Koduru model uses hyphens in some metabolite IDs (e.g. gln-L_e),
    which reframed encodes as __45__ when loading with flavor='bigg'.
    This creates duplicate pool reactions in the merged model (one for
    gln__L_e from iLP728, one for gln__45__L_e from Koduru), inflating
    MIP, MRO and MU counts.
 
    This function renames all __45__ IDs to their BiGG equivalent
    (__45__ -> __) directly in the reframed model object, before passing
    it to Community(). The XML is not modified.
 
    Parameters
    ----------
    lm_model : reframed CBModel  (loaded with load_cbmodel)
 
    Returns
    -------
    lm_model : same object, modified in-place
    n_renamed : int  number of metabolites renamed
    """
    import re
 
    # Collect all metabolite IDs with __45__
    to_rename = {
        m_id: m_id.replace('__45__', '__')
        for m_id in list(lm_model.metabolites.keys())
        if '__45__' in m_id
    }
 
    if not to_rename:
        return lm_model, 0
 
    # 1. Rename metabolites
    for old_id, new_id in to_rename.items():
        if new_id in lm_model.metabolites:
            # Target ID already exists — skip to avoid collision
            # (should not happen after full harmonization but be safe)
            continue
        met = lm_model.metabolites[old_id]
        met.id = new_id
        lm_model.metabolites[new_id] = met
        del lm_model.metabolites[old_id]
 
    # 2. Update stoichiometry in all reactions that reference renamed mets
    for rxn in lm_model.reactions.values():
        new_stoich = {}
        changed = False
        for m_id, coeff in rxn.stoichiometry.items():
            if m_id in to_rename and to_rename[m_id] not in rxn.stoichiometry:
                new_stoich[to_rename[m_id]] = coeff
                changed = True
            else:
                new_stoich[m_id] = coeff
        if changed:
            rxn.stoichiometry = new_stoich
 
    # 3. Rename exchange reactions whose ID contains __45__
    ex_to_rename = {
        r_id: r_id.replace('__45__', '__')
        for r_id in list(lm_model.reactions.keys())
        if '__45__' in r_id
    }
    for old_id, new_id in ex_to_rename.items():
        if new_id in lm_model.reactions:
            continue
        rxn = lm_model.reactions[old_id]
        rxn.id = new_id
        lm_model.reactions[new_id] = rxn
        del lm_model.reactions[old_id]
 
    return lm_model, len(to_rename)