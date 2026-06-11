#!/usr/bin/env python3
"""
ct_aliases.py
=============
Hand-curated crosswalk from each annotation method's label vocabulary to a
*canonical Cell Ontology (CT) node name* as it appears in the
CT_Ontology_Spreadsheet ("Celltype" column).

WHY THIS FILE EXISTS
--------------------
Method labels (Monaco / CellTypist / Azimuth / AIFI / DICE / HPCA) do NOT use
the same strings as the HIPC CT tree. To map a per-cell label onto the tree we
need a label -> CT-node lookup. This file provides that lookup, plus a
normalisation routine and a fuzzy fallback for labels not listed here.

The CT node a label points to can be at ANY depth (a terminal type like
"Treg", or an internal node like "T Cell" / "Monocyte" / "DC"). The builder
then expands that node to its full treeLevel1..6 ancestor path.

HOW TO EXTEND
-------------
ALIASES is { CT_canonical_name : [ method-label strings ] }.
- Add the raw method label string (any case / punctuation) to the right list.
- Normalisation (_norm) lowercases, collapses spaces, and strips most
  punctuation, so you usually don't need exact casing.
- Anything still unmatched lands in unmapped_labels_report.csv with a fuzzy
  best-guess for you to confirm.
"""

import re

# ---------------------------------------------------------------------------
# Canonical CT node names — MUST match the "Celltype" column of the ontology.
# (Kept here only as a sanity reference; the builder validates against the xlsx.)
# ---------------------------------------------------------------------------
CT_NODES = [
    "Blood Cell", "Platelet", "RBC", "HSC", "Doublet", "Leukocyte",
    "Lymphoid Cell", "NK Cell", "T Cell", "CD4 T Cell (ab)",
    "CD4 Naive / T Central Memory", "CD4 T Effector Memory", "Treg",
    "CD8 T Cell (ab)", "CD8 Naive / T Central Memory",
    "CD8 Cytotoxic / T Effector Memory", "gdT Cell", "MAIT Cell", "NKT Cell",
    "B Cell", "Effector B", "Plasma Cell", "Plasmablast", "Naive B Cell",
    "Memory B Cell", "Myeloid Cell", "Monocyte", "Classical Monocyte",
    "Non-Classical Monocyte", "Intermediate Monocyte", "Granulocyte",
    "Neutrophil", "Eosinophil", "Basophil", "Mast Cell", "DC",
    "Plasmacytoid DC", "Conventional DC 1", "Conventional DC 2",
]


def _norm(s: str) -> str:
    """Lowercase, strip punctuation to spaces, collapse whitespace."""
    if s is None:
        return ""
    s = str(s).strip().lower()
    s = re.sub(r"[._/\\,\-+()]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


# ---------------------------------------------------------------------------
# ALIASES: CT canonical name -> list of raw method labels that map to it.
# Covers Monaco(SingleR fine), Azimuth PBMC l1/l2/l3, CellTypist
# Immune_All_High/Low, AIFI L2/L3, DICE, HPCA. Extend freely.
# ---------------------------------------------------------------------------
ALIASES = {
    # ---- top / structural ----
    "Blood Cell": ["blood cell"],
    "Leukocyte": ["leukocyte", "immune cell", "pbmc"],
    "Lymphoid Cell": ["lymphoid cell", "lymphocyte", "lymphocytes",
                      "innate lymphoid cells", "ilc", "ilc precursor",
                      "ilc3", "ilc1", "ilc2", "cycling t/nk cell",
                      "cycling t nk cell", "t/nk cell"],
    "Myeloid Cell": ["myeloid cell", "myeloid", "myeloid cells", "mono dc",
                     "mono/dc", "mnp", "cycling myeloid cell",
                     "g2 m phase myeloid cell", "g2/m phase myeloid cell",
                     # macrophages have no dedicated node -> conservative ancestor
                     "macrophage", "macrophages", "m1 macrophage",
                     "m2 macrophage", "alveolar macrophages",
                     "alveolar macrophage", "microglia", "kupffer cells",
                     "erythrophagocytic macrophages",
                     "intermediate macrophages", "lyve1 macrophage",
                     "lyve1+ macrophage"],

    # ---- non-leukocyte blood ----
    "Platelet": ["platelet", "platelets", "megakaryocyte", "megakaryocytes",
                 "mk", "megakaryocytes/platelets", "platelet/megakaryocyte",
                 "megakaryocyte precursor", "megakaryocyte progenitor",
                 "mkp", "early mk", "early megakaryocyte"],
    "RBC": ["rbc", "erythrocyte", "erythrocytes", "red blood cell",
            "erythroid", "erythroblast", "erythroid cells", "ery",
            "late erythroid", "early erythroid", "mid erythroid",
            "early sox4+ erythroblast",
            "fetal hbg+ erythrocyte", "proerythroblast", "normoblast",
            "intermediate epcam+ erythroblast",
            "late hemoglobin+ erythroblast"],
    "HSC": ["hsc", "hspc", "hsc/mpp", "hematopoietic stem cell", "stem cells",
            "progenitor cells", "progenitor", "progenitor cell",
            "cd34+ progenitor", "hpc", "early lymphoid/t lymphoid",
            "cmp", "gmp", "mep", "clp", "clp cell", "cmp cell", "mpp",
            "lmpp", "etp", "elp", "memp", "hematopoietic precursor cell",
            "hsc_-g-csf", "hsc -g-csf", "hsc_cd34+", "hsc cd34+",
            "bm & prog.", "bm & prog", "bm prog",
            "megakaryocyte-erythroid-mast cell progenitor",
            "monocyte precursor", "eobama progenitor",
            "neutrophil-myeloid progenitor"],
    "Doublet": ["doublet", "doublets"],

    # ---- NK / NKT ----
    "NK Cell": ["nk cell", "nk cells", "natural killer cells",
                "natural killer cell", "nk", "cd56bright nk", "cd56dim nk",
                "nk_cell", "nk_56hi", "cd56 bright nk", "cd56 dim nk",
                "adaptive nk", "adaptive nk cell", "nk proliferating",
                "nk_proliferating", "proliferating nk cell", "cycling nk cells",
                "transitional nk", "tissue-resident nk cell",
                "cd16 nk cell", "cd56 nk cell", "cd16+ nk cells",
                "cd16- nk cells", "cd56bright nk cell", "cd56dim nk cell",
                "gzmk+ cd56dim nk cell", "gzmk- cd56dim nk cell",
                "isg+ cd56dim nk cell"],
    "NKT Cell": ["nkt cell", "nkt cells", "nkt", "nk t cell"],

    # ---- T cell parent (generic / lineage-only / unconventional) ----
    "T Cell": ["t cell", "t cells", "t_cells", "tcell", "abt cell",
               "abt (entry) cell", "abt entry cell",
               "alpha beta t cell", "double negative t cell",
               "double positive t cell", "dn t cell", "dp t cell",
               "double-negative thymocytes", "double negative thymocytes",
               "double-positive thymocytes", "double positive thymocytes",
               "other t", "proliferating t", "t proliferating",
               "proliferating t cell", "cycling t cell", "cycling t cells",
               "naive t cell", "memory t cell", "inf-activated t cell",
               "ifn-activated t cell", "isg+ t cell",
               "t(agonist)", "t agonist",
               "cd8aa", "cd8a/a", "cd8a a",
               "t cells, cd4", "t cells, cd8"],

    # CD4 lineage
    "CD4 T Cell (ab)": ["cd4 t cell", "cd4 t cells", "cd4 t", "cd4+ t cells",
                        "helper t cells", "th cells", "cd4tcell",
                        "t helper cells", "memory cd4 t cell",
                        "memory cd4 t cells", "cd4 memory t cell",
                        "isg+ memory cd4 t cell",
                        # polarized helper subsets: effector/memory CD4, place at parent
                        "th1 cells", "th2 cells", "th17 cells",
                        "th1 th17 cells", "th1/th17 cells",
                        "type 1 helper t cells", "type 2 helper t cells",
                        "type 17 helper t cells", "follicular helper t cells",
                        "tfh cell", "tfh cells",
                        "t cells cd4 tfh", "t cells cd4 th1", "t cells cd4 th2",
                        "t cells cd4 th17", "t cells cd4 th1 17",
                        "klrf1- gzmb+ cd27- memory cd4 t cell"],
    "CD4 Naive / T Central Memory": [
        "cd4 naive", "cd4 tcm", "cd4 naive t central memory",
        "naive cd4 t cells", "naive cd4 t cell", "central memory cd4 t cells",
        "tcm naive helper t cells", "tcm/naive helper t cells",
        "naive cd4 cells", "cd4 t naive", "core naive cd4 t cell",
        "cm cd4 t cell", "central memory cd4 t cell",
        "t cells cd4 naive", "t cells, cd4, naive",
        "t cells, cd4, naive, stimulated",
        "isg+ naive cd4 t cell", "sox4+ naive cd4 t cell"],
    "CD4 T Effector Memory": [
        "cd4 tem", "cd4 t effector memory", "effector memory cd4 t cells",
        "tem effector helper t cells", "tem/effector helper t cells",
        "terminal effector cd4 t cells", "cd4 ctl", "cd4 tem cell",
        "em cd4 t cell", "effector memory cd4 t cell",
        "cd4 cytotoxic t cell", "cd4 tcm/tem",
        "gzmb- cd27+ em cd4 t cell", "gzmb- cd27- em cd4 t cell",
        "memory cd4+ cytotoxic t cells", "memory cd4 cytotoxic t cells",
        "klrb1 cytotoxic cd4 t", "klrb1+ cytotoxic cd4 t",
        "tem/effector helper t cells pd1+", "tem effector helper t cells pd1"],
    # NOTE: the HIPC tree nests Treg under CD4 T Cell (ab) (CL:0000815 here sits
    # on the CD4 branch). CD4 Tregs map here. CD8 FOXP3+ Tregs are bona fide but
    # genuinely CD8 lineage (Churlaud 2015; Yu 2025 Trends Immunol) -> placing
    # them here would mis-call CD4 at treeLevel5, so they go under CD8 T Cell (ab).
    "Treg": ["treg", "tregs", "regulatory t cells", "regulatory t cell",
             "t regulatory cells", "treg cell", "treg(diff)", "treg diff",
             "t cells cd4 memory treg", "t cells cd4 naive treg",
             "t cells, cd4, memory treg", "t cells, cd4, naive treg",
             "memory treg", "naive treg", "memory cd4 treg", "naive cd4 treg",
             "gzmk+ memory cd4 treg", "klrb1+ memory cd4 treg"],

    # CD8 lineage
    "CD8 T Cell (ab)": ["cd8 t cell", "cd8 t cells", "cd8 t", "cd8+ t cells",
                        "cytotoxic t cells", "cd8tcell",
                        "memory cd8 t cell", "memory cd8 t cells",
                        "isg+ memory cd8 t cell",
                        "klrb1 cd8 t cell", "klrb1+ cd8 t cell",
                        # CD8 FOXP3+ Tregs: bona fide regulatory but CD8 lineage
                        "memory cd8 treg", "klrb1+ memory cd8 treg",
                        "cd8a/b(entry)", "cd8a/b entry", "cd8ab entry"],
    "CD8 Naive / T Central Memory": [
        "cd8 naive", "cd8 tcm", "cd8 naive t central memory",
        "naive cd8 t cells", "naive cd8 t cell", "central memory cd8 t cells",
        "tcm naive cytotoxic t cells", "tcm/naive cytotoxic t cells",
        "core naive cd8 t cell", "cm cd8 t cell",
        "central memory cd8 t cell", "t cells cd8 naive",
        "t cells, cd8, naive", "t cells, cd8, naive, stimulated",
        "cd8 t naive", "isg+ naive cd8 t cell", "sox4+ naive cd8 t cell"],
    "CD8 Cytotoxic / T Effector Memory": [
        "cd8 tem", "cd8 cytotoxic t effector memory",
        "effector memory cd8 t cells", "terminal effector cd8 t cells",
        "tem effector cytotoxic t cells", "tem/effector cytotoxic t cells",
        "tem temra cytotoxic t cells", "tem/temra cytotoxic t cells",
        "tem trm cytotoxic t cells", "tem/trm cytotoxic t cells",
        "trm cytotoxic t cells",
        "cd8 temra", "cd8 teff", "gzmk cd8 tem", "gzmb cd8 tem",
        "gzmb cd8 t cell", "gzmk cd8 t cell",
        "em cd8 t cell", "effector memory cd8 t cell", "cd8 ctl",
        "cytotoxic cd8 t cell",
        "gzmk+ cd27+ em cd8 t cell", "gzmk- cd27+ em cd8 t cell",
        "klrf1+ gzmb+ cd27- em cd8 t cell",
        "klrf1- gzmb+ cd27- em cd8 t cell"],

    # innate-like T
    "gdT Cell": ["gdt cell", "gdt cells", "gdt", "gd t cells", "gd t cell",
                 "gamma delta t cells", "gamma-delta t cells",
                 "crtam+ gamma-delta t cells", "crtam+ gamma delta t cells",
                 "vd2 gd t cells", "non vd2 gd t cells", "non-vd2 gd t cells",
                 "vd1 gd t cell", "vd2 gd t cell",
                 "gzmb+ vd2 gdt", "gzmk+ vd2 gdt", "naive vd1 gdt",
                 "sox4+ vd1 gdt", "klrf1+ effector vd1 gdt",
                 "klrf1- effector vd1 gdt"],
    "MAIT Cell": ["mait cell", "mait cells", "mait", "mait_cell",
                  "cd4 mait", "cd8 mait", "isg+ mait"],

    # ---- B cells ----
    "B Cell": ["b cell", "b cells", "b_cells", "bcell", "b lineage",
               "b-cell lineage", "b cell lineage",
               "proliferating b", "b proliferating", "b cells, naive",
               "follicular b cells", "germinal center b cells",
               "germinal center b cell", "proliferative germinal center b cells",
               "preb cell", "pre-b cell", "pre-b_cell_cd34-",
               "pre b cell cd34-", "large pre-b cells", "small pre-b cells",
               "prob cell", "pro-b cell", "pro-b cells", "pro-b_cell_cd34+",
               "pre-pro-b cells", "pre pro b cells",
               "double-negative b cell"],
    "Naive B Cell": ["naive b cell", "naive b cells", "b naive",
                     "core naive b cell", "transitional b cells",
                     "transitional b cell", "b cells naive",
                     "isg+ naive b cell"],
    "Memory B Cell": ["memory b cell", "memory b cells", "b memory",
                      "b intermediate", "switched memory b cells",
                      "non switched memory b cells",
                      "non-switched memory b cells",
                      "unswitched memory b cells", "exhausted b cells",
                      "core memory b cell", "age associated b cells",
                      "abc", "atypical memory b cells",
                      "activated memory b cell", "early memory b cell",
                      "cd95 memory b cell", "type 2 polarized memory b cell"],
    "Effector B": ["effector b", "effector b cell",
                   "cd27+ effector b cell", "cd27- effector b cell",
                   "antibody secreting cells", "asc"],
    "Plasma Cell": ["plasma cell", "plasma cells", "plasma", "long lived plasma cell"],
    "Plasmablast": ["plasmablast", "plasmablasts", "plasmablast cell"],

    # ---- Monocytes ----
    "Monocyte": ["monocyte", "monocytes", "mono", "cd14 mono", "cd16 mono",
                 "mono-mac", "mono mac", "cycling monocytes", "cycling monocyte"],
    "Classical Monocyte": ["classical monocyte", "classical monocytes",
                           "cd14 monocyte", "cd14+ monocyte", "cd14 mono",
                           "cd14+ mono", "monocytes cd14+", "monocytes, cd14+",
                           "core cd14 monocyte", "cm cd14 monocyte",
                           "il1b+ cd14 monocyte", "isg+ cd14 monocyte"],
    "Non-Classical Monocyte": ["non classical monocyte",
                               "non-classical monocyte",
                               "non classical monocytes",
                               "non-classical monocytes", "cd16 monocyte",
                               "cd16+ monocyte", "cd16 mono", "cd16+ mono",
                               "monocytes cd16+", "monocytes, cd16+",
                               "core cd16 monocyte", "ncm",
                               "c1q+ cd16 monocyte", "isg+ cd16 monocyte"],
    "Intermediate Monocyte": ["intermediate monocyte",
                              "intermediate monocytes", "cd14 cd16 monocyte",
                              "im"],

    # ---- Granulocytes ----
    "Granulocyte": ["granulocyte", "granulocytes", "myelocytes", "myelocyte",
                    "promyelocyte", "promyelocytes",
                    "baeomap cell", "eobama progenitor cell"],
    "Neutrophil": ["neutrophil", "neutrophils", "low density neutrophils",
                   "low-density neutrophils", "neut", "pmn",
                   "neutrophil precursor"],
    "Eosinophil": ["eosinophil", "eosinophils", "eos"],
    "Basophil": ["basophil", "basophils", "low density basophils",
                 "low-density basophils", "baso"],
    "Mast Cell": ["mast cell", "mast cells", "mast", "mastcell"],

    # ---- Dendritic cells ----
    "DC": ["dc", "dendritic cell", "dendritic cells", "myeloid dendritic cells",
           "mdc", "conventional dc", "cdc", "dc precursor", "asdc",
           "ax dendritic cell", "mregdc", "mreg dc", "cycling dcs", "cycling dc",
           "migratory dcs", "migratory dc", "transitional dc"],
    "Plasmacytoid DC": ["plasmacytoid dc", "plasmacytoid dcs",
                        "plasmacytoid dendritic cells",
                        "plasmacytoid dendritic cell", "pdc", "pdcs",
                        "pdc precursor"],
    "Conventional DC 1": ["conventional dc 1", "cdc1", "cdc 1", "dc1",
                          "conventional dendritic cell 1", "clec9a dc",
                          "bdca3 dc"],
    "Conventional DC 2": ["conventional dc 2", "cdc2", "cdc 2", "dc2",
                          "conventional dendritic cell 2", "cd1c dc",
                          "bdca1 dc", "dc3", "inflammatory dc",
                          "cd14+ cdc2", "hla-drhi cdc2", "hla drhi cdc2",
                          "isg+ cdc2"],
}

# Build normalised lookup: normalised-label -> CT canonical name
_LOOKUP = {}
for ct_name, labels in ALIASES.items():
    _LOOKUP[_norm(ct_name)] = ct_name          # the CT name itself maps to itself
    for lab in labels:
        _LOOKUP[_norm(lab)] = ct_name


def exact_match(label: str):
    """Return CT canonical name for a method label via the alias table, else None."""
    return _LOOKUP.get(_norm(label))


def fuzzy_match(label: str, ct_names, cutoff: float = 0.84):
    """
    Best-guess CT node for an unaliased label using token-set + difflib ratio.
    Returns (ct_name, score) or (None, best_score). Conservative by design;
    anything below `cutoff` should be human-reviewed.
    """
    import difflib
    nl = _norm(label)
    if not nl:
        return None, 0.0
    nl_tokens = set(nl.split())
    best, best_score = None, 0.0
    # candidate pool = CT names + every alias string, each pointing to a CT name
    candidates = {_norm(c): c for c in ct_names}
    for k, ct in _LOOKUP.items():
        candidates[k] = ct
    for cand_norm, ct in candidates.items():
        seq = difflib.SequenceMatcher(None, nl, cand_norm).ratio()
        ct_tokens = set(cand_norm.split())
        jac = (len(nl_tokens & ct_tokens) / len(nl_tokens | ct_tokens)
               if (nl_tokens | ct_tokens) else 0.0)
        score = 0.6 * seq + 0.4 * jac
        if score > best_score:
            best, best_score = ct, score
    if best_score >= cutoff:
        return best, round(best_score, 3)
    return None, round(best_score, 3)
