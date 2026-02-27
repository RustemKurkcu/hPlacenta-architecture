# ==============================================================================
# config/misi_v2_gene_sets.R
# MISI v2.0 — Metabolic Immune Spatial Susceptibility Index
# Gene Set Definitions for Dual-Axis Scoring Framework
#
# Architecture:
#   Axis 1: S-MISI (Susceptibility) — 5 dimensions
#   Axis 2: V-MISI (Severity)       — 5 dimensions
#   Composite: MISI_v2 = S-MISI × (1 + V-MISI)
#
# Author: Shan Kurkcu, UConn Health (Dr. Das Lab)
# Version: 2.0
# Date: 2025-07
#
# References:
#   Franklin et al. 2024 (mBio) — Toxic switch / MegL
#   Parhi et al. 2022 (Cell Reports) — Fap2 / Gal-GalNAc
#   Hoo et al. 2024 (Cell Systems) — Nutritional immunity / PLD1
#   Greenbaum et al. 2023 (Nature) — Spatial atlas / EVT zone
#   Britton et al. 2022 (mBio) — Rnf / FadA
#   Suwatthee et al. 2025 (Nat Commun) — Fap2 structure
# ==============================================================================

# ── AXIS 1: SUSCEPTIBILITY (S-MISI) ──────────────────────────────────────────
# Measures conditions that make a cell/niche permissive to Fn colonization

# Dimension 1: Tropism Potential (weight = 0.25)
# Rationale: Fap2 is the primary tropism determinant (Parhi 2022);
#            FadA provides secondary adhesion and barrier breach (Rubinstein 2013).
#            The 65:35 split reflects primacy of Fap2 in placental colonization.
SMISI_TROPISM_FAP2 <- c(
  "GALNT1",     # ★ Primary 'address label' — Gal-GalNAc display
  "GALNT2",     # Ubiquitous GalNAc-T
  "GALNT3",     # Expressed in trophoblasts
  "GALNT7",     # Mucin-type O-glycosylation
  "GALNT10",    # Tissue-specific expression
  "GALNT14",    # Placental expression
  "C1GALT1",    # Core 1 synthase (T-synthase)
  "C1GALT1C1",  # Core 1 synthase chaperone (Cosmc)
  "GCNT1",      # Core 2 GlcNAc-T
  "GCNT3",      # Core 2/4 GlcNAc-T
  "ST6GALNAC1", # Sialyltransferase
  "ST6GALNAC2", # Sialyltransferase
  "B3GNT6"      # Core 3 synthase
)

SMISI_TROPISM_FADA <- c(
  "CDH1",    # ★ E-cadherin — primary FadA receptor
  "CTNNB1",  # Beta-catenin — signaling mediator
  "CTNNA1",  # Alpha-catenin — actin linker
  "CTNND1",  # Delta-catenin (p120) — stability
  "JUP",     # Junction plakoglobin (gamma-catenin)
  "CDH2",    # N-cadherin
  "CDH3"     # P-cadherin (placental)
)

# Dimension 2: Nutrient Availability — "Dinner Bell" (weight = 0.25)
# Rationale: PLD1 is the primary EA producer and most direct measure of the
#            'dinner bell' signal. GDPD1 recycles membrane scraps into EA.
#            PLD1 downregulation is the most consistent marker of nutritional
#            immunity (Hoo 2024; Bhalla 2024).
SMISI_NUTRIENT_PLD1  <- c("PLD1")   # ★ KEY BAIT MARKER — 50% sub-weight
SMISI_NUTRIENT_GDPD1 <- c("GDPD1") # ★ Scavenger enzyme — 20% sub-weight
SMISI_NUTRIENT_EA_PATHWAY <- c(
  "ETNK1",   # Ethanolamine kinase 1
  "ETNK2",   # Ethanolamine kinase 2
  "PCYT2",   # CTP:phosphoethanolamine cytidylyltransferase
  "SELENOI", # Ethanolamine phosphotransferase (EPT)
  "EPT1",    # Ethanolamine phosphotransferase 1
  "CEPT1",   # Choline/ethanolamine phosphotransferase
  "PEMT"     # PE N-methyltransferase
)

# Dimension 3: Remodeling Highway (weight = 0.20) [NEW in v2.0]
# Rationale: MMP-mediated ECM degradation creates physical pathways that Fn
#            exploits for deeper tissue penetration. MMP-2 and MMP-9 are the
#            primary gelatinases active in placental invasion (Greenbaum 2023).
SMISI_REMODELING_MMP_CORE <- c(
  "MMP2",   # ★ Gelatinase A — constitutive
  "MMP9",   # ★ Gelatinase B — inducible
  "MMP14"   # ★ MT1-MMP — membrane-anchored
)

SMISI_REMODELING_ECM_DEGRADATION <- c(
  "ADAMTS4",  # Aggrecanase-1
  "ADAMTS5",  # Aggrecanase-2
  "CTSK",     # Cathepsin K
  "CTSL",     # Cathepsin L
  "CTSB",     # Cathepsin B
  "PLAUR",    # uPAR — plasminogen activator receptor
  "SERPINE1"  # PAI-1 — plasminogen activator inhibitor
)

# Dimension 4: Immune Tolerance — "Permissive Shield" (weight = 0.20) [NEW in v2.0]
# Rationale: Immune tolerance programs create a permissive environment for Fn
#            colonization. The TIGIT axis is specifically relevant because Fap2
#            directly activates TIGIT to suppress NK and T cell killing (Parhi 2022).
#            NK cytotoxicity is subtracted because active surveillance eliminates Fn.
SMISI_TOLERANCE_CHECKPOINT <- c(
  "HLA-G",    # ★ EVT immune tolerance — primary checkpoint
  "CD274",    # PD-L1
  "PDCD1LG2", # PD-L2
  "IDO1",     # Indoleamine 2,3-dioxygenase 1
  "TGFB1",    # TGF-beta 1
  "IL10",     # Interleukin-10
  "LGALS9",   # Galectin-9 (TIM-3 ligand)
  "HAVCR2",   # TIM-3
  "VSIR",     # VISTA
  "LILRB1",   # ILT2
  "LILRB2"    # ILT4
)

SMISI_TOLERANCE_TIGIT_AXIS <- c(
  "TIGIT",  # ★ Fap2-activated killing-suppressor on NK/T cells
  "PVR",    # CD155 — TIGIT ligand
  "PVRL2"   # CD112 — TIGIT ligand
)

SMISI_NK_CYTOTOXIC_NEGATIVE <- c(
  "NKG7",  # NK granule protein 7
  "GNLY",  # Granulysin
  "PRF1",  # Perforin
  "GZMB",  # Granzyme B
  "GZMA"   # Granzyme A
)

# Dimension 5: Barrier Integrity (weight = -0.10, NEGATIVE — protective)
# Rationale: An intact syncytiotrophoblast barrier and tight junctions physically
#            prevent Fn access to underlying tissue.
SMISI_BARRIER_STB <- c(
  "CSH1",    # Chorionic somatomammotropin hormone 1
  "CSH2",    # Chorionic somatomammotropin hormone 2
  "CYP19A1", # Aromatase
  "PSG1",    # Pregnancy-specific glycoprotein 1
  "PSG3"     # Pregnancy-specific glycoprotein 3
)

SMISI_BARRIER_TIGHT_JUNCTIONS <- c(
  "TJP1",  # ZO-1
  "TJP2",  # ZO-2
  "OCLN",  # Occludin
  "CLDN1", # Claudin-1
  "CLDN3", # Claudin-3
  "CLDN4"  # Claudin-4
)


# ── AXIS 2: SEVERITY (V-MISI) ────────────────────────────────────────────────
# Measures potential for pathological consequences once Fn infection is established

# Dimension 1: Toxic Switch Potential (weight = 0.30) [↑ from 0.05 in v1.0]
# Rationale: The toxic switch is the critical mechanistic link to preeclampsia.
#            Tissues with high baseline H₂S production and low detoxification
#            capacity are most vulnerable to additional H₂S from bacterial MegL.
#            (Franklin 2024; H₂S review PMC12486883)
VMISI_H2S_PRODUCTION <- c(
  "CBS",   # ★ Cystathionine beta-synthase
  "CTH",   # ★ Cystathionine gamma-lyase (dysregulated in PE)
  "MPST"   # Mercaptopyruvate sulfurtransferase
)

VMISI_H2S_DETOX <- c(
  "SQOR",  # Sulfide quinone oxidoreductase
  "ETHE1", # Persulfide dioxygenase
  "TST"    # Thiosulfate sulfurtransferase
)

VMISI_AMMONIA_CLEARANCE <- c(
  "GLUL",  # Glutamine synthetase
  "GLUD1", # Glutamate dehydrogenase 1
  "CPS1",  # Carbamoyl phosphate synthetase 1
  "OTC",   # Ornithine transcarbamylase
  "ASS1",  # Argininosuccinate synthase
  "ASL",   # Argininosuccinate lyase
  "ARG1"   # Arginase 1
)

# Dimension 2: Inflammatory Amplification (weight = 0.25) [NEW in v2.0]
# Rationale: Fn triggers placental inflammation through TLR4-mediated NF-κB
#            activation (Garcia-So 2019). The inflammatory cascade amplifies
#            tissue damage and recruits additional immune cells.
VMISI_NFKB_TLR4 <- c(
  "TLR4",  # ★ Fn LPS receptor
  "NFKB1", # NF-κB p50
  "NFKB2", # NF-κB p52
  "RELA",  # NF-κB p65
  "MYD88", # MyD88 adaptor
  "IRAK4", # IL-1 receptor-associated kinase 4
  "TRAF6"  # TNF receptor-associated factor 6
)

VMISI_MYELOID_INFLAMMATION <- c(
  "IL1B",   # Interleukin-1 beta
  "TNF",    # Tumor necrosis factor alpha
  "NFKBIA", # IκBα
  "CXCL8",  # IL-8
  "S100A8", # Calgranulin A
  "S100A9", # Calgranulin B
  "CCL2",   # MCP-1
  "CCL3",   # MIP-1α
  "CCL4"    # MIP-1β
)

VMISI_IFN_RESPONSE <- c(
  "ISG15", # ISG15 ubiquitin-like modifier
  "IFIT1", # Interferon-induced protein with tetratricopeptide repeats 1
  "IFIT2", # IFIT2
  "IFIT3", # IFIT3
  "MX1",   # MX dynamin-like GTPase 1
  "OAS1",  # 2'-5'-oligoadenylate synthetase 1
  "OAS2",  # OAS2
  "IFI6",  # Interferon alpha-inducible protein 6
  "IFI27"  # Interferon alpha-inducible protein 27
)

# Dimension 3: Vascular Damage Potential (weight = 0.25) [NEW in v2.0]
# Rationale: Endothelial activation and angiogenic imbalance are the proximate
#            causes of preeclampsia symptoms. H₂S-induced endothelial damage
#            manifests as increased VCAM1/ICAM1 and elevated sFlt-1.
VMISI_ENDOTHELIAL_ACTIVATION <- c(
  "VCAM1", # ★ Vascular cell adhesion molecule 1
  "ICAM1", # Intercellular adhesion molecule 1
  "SELE",  # E-selectin
  "ENG",   # Endoglin (CD105)
  "EDN1"   # Endothelin-1
)

VMISI_ANGIOGENIC_IMBALANCE <- c(
  "FLT1",  # ★ VEGFR1 / sFlt-1 (anti-angiogenic in PE)
  "VEGFA", # VEGF-A (pro-angiogenic — high = protective)
  "PGF",   # Placental growth factor (pro-angiogenic — high = protective)
  "KDR",   # VEGFR2
  "NOS3"   # eNOS
)

# Dimension 4: Immune Dysregulation (weight = 0.15) [NEW in v2.0]
# Rationale: While immune tolerance enables Fn colonization (susceptibility),
#            subsequent immune activation in response to infection and toxins
#            drives tissue damage (severity).
VMISI_NK_EXCESS <- c(
  "GZMB",  # Granzyme B
  "PRF1",  # Perforin
  "GNLY",  # Granulysin
  "IFNG",  # Interferon gamma
  "KLRK1" # NKG2D
)

VMISI_M1_MACROPHAGE <- c(
  "IL1B",  # IL-1β
  "TNF",   # TNF-α
  "NOS2",  # iNOS
  "CD80",  # B7-1
  "CD86"   # B7-2
)

# Dimension 5: Hypoxia Response (weight = 0.05) [NEW in v2.0]
# Rationale: Hypoxia is both a consequence of vascular damage and an amplifier
#            of pathology. HIF1A activation is a hallmark of preeclampsia and
#            may create a metabolic environment favoring Fn survival (anaerobe).
VMISI_HYPOXIA <- c(
  "HIF1A",  # Hypoxia-inducible factor 1-alpha
  "VEGFA",  # VEGF-A
  "SLC2A1", # GLUT1
  "LDHA",   # Lactate dehydrogenase A
  "CA9",    # Carbonic anhydrase IX
  "BNIP3",  # BCL2/adenovirus E1B 19 kDa protein-interacting protein 3
  "EGLN1"   # Prolyl hydroxylase domain protein 2
)


# ── COMBINED MISI v2.0 GENE SET LIST ─────────────────────────────────────────

MISI_V2_GENESETS <- list(

  # ── S-MISI components ──
  S_Tropism_Fap2          = SMISI_TROPISM_FAP2,
  S_Tropism_FadA          = SMISI_TROPISM_FADA,
  S_Nutrient_PLD1         = SMISI_NUTRIENT_PLD1,
  S_Nutrient_GDPD1        = SMISI_NUTRIENT_GDPD1,
  S_Nutrient_EA_Pathway   = SMISI_NUTRIENT_EA_PATHWAY,
  S_Remodeling_MMP_Core   = SMISI_REMODELING_MMP_CORE,
  S_Remodeling_ECM        = SMISI_REMODELING_ECM_DEGRADATION,
  S_Tolerance_Checkpoint  = SMISI_TOLERANCE_CHECKPOINT,
  S_Tolerance_TIGIT       = SMISI_TOLERANCE_TIGIT_AXIS,
  S_NK_Cytotoxic_Neg      = SMISI_NK_CYTOTOXIC_NEGATIVE,
  S_Barrier_STB           = SMISI_BARRIER_STB,
  S_Barrier_TJ            = SMISI_BARRIER_TIGHT_JUNCTIONS,

  # ── V-MISI components ──
  V_H2S_Production        = VMISI_H2S_PRODUCTION,
  V_H2S_Detox             = VMISI_H2S_DETOX,
  V_Ammonia_Clearance     = VMISI_AMMONIA_CLEARANCE,
  V_NFkB_TLR4             = VMISI_NFKB_TLR4,
  V_Myeloid_Inflammation  = VMISI_MYELOID_INFLAMMATION,
  V_IFN_Response          = VMISI_IFN_RESPONSE,
  V_Endothelial_Act       = VMISI_ENDOTHELIAL_ACTIVATION,
  V_Angiogenic_Imbalance  = VMISI_ANGIOGENIC_IMBALANCE,
  V_NK_Excess             = VMISI_NK_EXCESS,
  V_M1_Macrophage         = VMISI_M1_MACROPHAGE,
  V_Hypoxia               = VMISI_HYPOXIA
)

# ── MISI v2.0 WEIGHTS ────────────────────────────────────────────────────────

# S-MISI dimension weights (must sum to 1.0 including negative barrier)
SMISI_DIM_WEIGHTS <- c(
  Tropism    = 0.25,   # Fap2 (65%) + FadA (35%)
  Nutrient   = 0.25,   # PLD1 (50%) + EA_pathway (30%) + GDPD1 (20%)
  Remodeling = 0.20,   # MMP_core (60%) + ECM_degradation (40%)
  Tolerance  = 0.20,   # Checkpoint (50%) + TIGIT (30%) - NK_cytotoxic (20%)
  Barrier    = -0.10   # STB + TJ (negative — protective)
)

# Sub-weights within Tropism dimension
SMISI_TROPISM_SUBWEIGHTS <- c(Fap2 = 0.65, FadA = 0.35)

# Sub-weights within Nutrient dimension
SMISI_NUTRIENT_SUBWEIGHTS <- c(PLD1 = 0.50, EA_Pathway = 0.30, GDPD1 = 0.20)

# Sub-weights within Remodeling dimension
SMISI_REMODELING_SUBWEIGHTS <- c(MMP_Core = 0.60, ECM_Degradation = 0.40)

# Sub-weights within Tolerance dimension (NK is negative)
SMISI_TOLERANCE_SUBWEIGHTS <- c(Checkpoint = 0.50, TIGIT = 0.30, NK_Neg = -0.20)

# V-MISI dimension weights (must sum to 1.0)
VMISI_DIM_WEIGHTS <- c(
  ToxicSwitch  = 0.30,  # H2S production + detox deficit + ammonia vulnerability
  Inflammation = 0.25,  # NF-κB/TLR4 + myeloid + IFN
  Vascular     = 0.25,  # Endothelial activation + angiogenic imbalance
  ImmuneDysreg = 0.15,  # NK excess + M1 macrophage
  Hypoxia      = 0.05   # HIF1A pathway
)

# Sub-weights within Toxic Switch dimension
VMISI_TOXICSWITCH_SUBWEIGHTS <- c(
  H2S_Production    = 0.40,
  H2S_Detox_Deficit = 0.30,  # NOTE: scored as INVERSE (low detox = high vulnerability)
  Ammonia_Vuln      = 0.30   # NOTE: scored as INVERSE (low clearance = high vulnerability)
)

# Sub-weights within Inflammation dimension
VMISI_INFLAMMATION_SUBWEIGHTS <- c(NFkB_TLR4 = 0.50, Myeloid = 0.30, IFN = 0.20)

# Sub-weights within Vascular dimension
VMISI_VASCULAR_SUBWEIGHTS <- c(Endothelial = 0.50, Angiogenic = 0.50)

# Sub-weights within Immune Dysregulation dimension
VMISI_IMMUNEDYSREG_SUBWEIGHTS <- c(NK_Excess = 0.50, M1_Macro = 0.50)

# ── TEMPORAL GATING WEIGHTS ──────────────────────────────────────────────────
# Gestational age-specific multipliers for S-MISI
# Based on PLD1 expression dynamics and EVT invasion activity
# (Greenbaum 2023; pipeline timecourse analysis)
MISI_V2_TEMPORAL_WEIGHTS <- c(
  "6"  = 0.90,
  "7"  = 1.10,
  "8"  = 1.20,   # Peak vulnerability window
  "9"  = 1.10,
  "10" = 0.90,
  "11" = 0.75,
  "12" = 0.60,
  "13" = 0.50,
  "14" = 0.45,
  "15" = 0.40,
  "16" = 0.35,
  "17" = 0.30,
  "18" = 0.25,
  "19" = 0.20,
  "20" = 0.20
)

# ── HELPER: get all genes for a MISI v2.0 axis ───────────────────────────────

#' Get all unique genes in S-MISI
get_smisi_genes <- function() {
  unique(unlist(list(
    SMISI_TROPISM_FAP2, SMISI_TROPISM_FADA,
    SMISI_NUTRIENT_PLD1, SMISI_NUTRIENT_GDPD1, SMISI_NUTRIENT_EA_PATHWAY,
    SMISI_REMODELING_MMP_CORE, SMISI_REMODELING_ECM_DEGRADATION,
    SMISI_TOLERANCE_CHECKPOINT, SMISI_TOLERANCE_TIGIT_AXIS,
    SMISI_NK_CYTOTOXIC_NEGATIVE,
    SMISI_BARRIER_STB, SMISI_BARRIER_TIGHT_JUNCTIONS
  )))
}

#' Get all unique genes in V-MISI
get_vmisi_genes <- function() {
  unique(unlist(list(
    VMISI_H2S_PRODUCTION, VMISI_H2S_DETOX, VMISI_AMMONIA_CLEARANCE,
    VMISI_NFKB_TLR4, VMISI_MYELOID_INFLAMMATION, VMISI_IFN_RESPONSE,
    VMISI_ENDOTHELIAL_ACTIVATION, VMISI_ANGIOGENIC_IMBALANCE,
    VMISI_NK_EXCESS, VMISI_M1_MACROPHAGE,
    VMISI_HYPOXIA
  )))
}

#' Get all unique genes in MISI v2.0
get_misi_v2_genes <- function() {
  unique(c(get_smisi_genes(), get_vmisi_genes()))
}

#' Print MISI v2.0 summary
print_misi_v2_summary <- function() {
  cat("\n", strrep("=", 65), "\n")
  cat("MISI v2.0 — Gene Set Summary\n")
  cat(strrep("=", 65), "\n\n")
  cat("S-MISI (Susceptibility Axis):\n")
  cat(sprintf("  Dim 1 Tropism    (0.25): Fap2=%d genes, FadA=%d genes\n",
              length(SMISI_TROPISM_FAP2), length(SMISI_TROPISM_FADA)))
  cat(sprintf("  Dim 2 Nutrient   (0.25): PLD1=%d, EA=%d, GDPD1=%d genes\n",
              length(SMISI_NUTRIENT_PLD1), length(SMISI_NUTRIENT_EA_PATHWAY),
              length(SMISI_NUTRIENT_GDPD1)))
  cat(sprintf("  Dim 3 Remodeling (0.20): MMP=%d, ECM=%d genes\n",
              length(SMISI_REMODELING_MMP_CORE), length(SMISI_REMODELING_ECM_DEGRADATION)))
  cat(sprintf("  Dim 4 Tolerance  (0.20): Checkpoint=%d, TIGIT=%d, NK_neg=%d genes\n",
              length(SMISI_TOLERANCE_CHECKPOINT), length(SMISI_TOLERANCE_TIGIT_AXIS),
              length(SMISI_NK_CYTOTOXIC_NEGATIVE)))
  cat(sprintf("  Dim 5 Barrier   (-0.10): STB=%d, TJ=%d genes\n",
              length(SMISI_BARRIER_STB), length(SMISI_BARRIER_TIGHT_JUNCTIONS)))
  cat(sprintf("  Total S-MISI unique genes: %d\n\n", length(get_smisi_genes())))

  cat("V-MISI (Severity Axis):\n")
  cat(sprintf("  Dim 1 ToxicSwitch  (0.30): H2S_prod=%d, H2S_detox=%d, NH3=%d genes\n",
              length(VMISI_H2S_PRODUCTION), length(VMISI_H2S_DETOX),
              length(VMISI_AMMONIA_CLEARANCE)))
  cat(sprintf("  Dim 2 Inflammation (0.25): NFkB=%d, Myeloid=%d, IFN=%d genes\n",
              length(VMISI_NFKB_TLR4), length(VMISI_MYELOID_INFLAMMATION),
              length(VMISI_IFN_RESPONSE)))
  cat(sprintf("  Dim 3 Vascular     (0.25): Endothelial=%d, Angiogenic=%d genes\n",
              length(VMISI_ENDOTHELIAL_ACTIVATION), length(VMISI_ANGIOGENIC_IMBALANCE)))
  cat(sprintf("  Dim 4 ImmuneDysreg (0.15): NK_excess=%d, M1=%d genes\n",
              length(VMISI_NK_EXCESS), length(VMISI_M1_MACROPHAGE)))
  cat(sprintf("  Dim 5 Hypoxia      (0.05): %d genes\n", length(VMISI_HYPOXIA)))
  cat(sprintf("  Total V-MISI unique genes: %d\n\n", length(get_vmisi_genes())))

  cat(sprintf("Total MISI v2.0 unique genes: %d\n", length(get_misi_v2_genes())))
  cat(strrep("=", 65), "\n\n")
}

if (interactive()) print_misi_v2_summary()

cat("✓ MISI v2.0 gene sets loaded\n")
cat("  S-MISI genes: ", length(get_smisi_genes()), "\n")
cat("  V-MISI genes: ", length(get_vmisi_genes()), "\n")
cat("  Total unique: ", length(get_misi_v2_genes()), "\n\n")