"""
Themed GO GSEA NES panels for OSD-564 (read from GO_results/GSEA_{BP,CC,MF}_results.csv).

GO MF panels (matched by GO ID):
  * TNT / F-actin signature
  * Extracellular vesicle signature
  * Amyloid-beta / tau signature (legacy framing; kept only so the existing file is
    regenerated with the harmonized style and correct OSD-564 label)

Supplementary harmonized panels (shared with OSD-562):
  * GO BP "protrusion and exocytosis" (matched by term name):
      Actin Filament Organization; Axonogenesis; Axon Guidance;
      Neuron Projection Morphogenesis; Regulated Exocytosis; Synaptic Vesicle Recycling
  * GO CC compartments (matched by term name):
      Recycling Endosome Membrane; Endocytic Vesicle; Clathrin-Coated Vesicle;
      Exocytic Vesicle; Mitochondrial Inner Membrane; Ribosome
"""
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()

from srep_nes_panel import nes_panel, XLIM_DEFAULT

warnings.filterwarnings('ignore')

# every panel in this script shares the (-3.2, 3.2) axis of the mitochondrial /
# metabolic panels (the CC panel holds Ribosome at NES = -2.27)
TITLE_XLIM = XLIM_DEFAULT

OUT_DIR = 'GO_results'
CSV = {
    'BP': 'GO_results/GSEA_BP_results.csv',
    'CC': 'GO_results/GSEA_CC_results.csv',
    'MF': 'GO_results/GSEA_MF_results.csv',
}

# ──────────────────────────────────────────────────────────────────
# Load GSEA results
# ──────────────────────────────────────────────────────────────────
print("Loading GSEA GO results ...")
GSEA = {}
for ont, path in CSV.items():
    d = pd.read_csv(path, index_col=0)
    d['Term'] = d['Term'].str.strip()
    GSEA[ont] = d
    print(f"  GO {ont}: {len(d)} terms")

# ──────────────────────────────────────────────────────────────────
# Themed pathway sets
#   'go_ids' -> matched by GO ID inside the Term string
#   'names'  -> matched by term name (case-insensitive; Term = "<name> (GO:xxxxxxx)")
# ──────────────────────────────────────────────────────────────────
THEMES = [
    {
        'key':   'TNT_Factin',
        'ont':   'MF',
        'title': 'GO MF TNT / F-actin terms, OSD-564 (n = 5 vs 5)',
        'go_ids': [
            'GO:0003779',   # Actin Binding
            'GO:0003785',   # Actin Monomer Binding
            'GO:0031005',   # Filamin Binding
            'GO:0051393',   # Alpha-actinin Binding
            'GO:0017022',   # Myosin Binding
            'GO:0031489',   # Myosin V Binding
            'GO:0017166',   # Vinculin Binding
            'GO:0031267',   # Small GTPase Binding
            'GO:0005085',   # Guanyl-nucleotide Exchange Factor Activity
            'GO:0005546',   # Phosphatidylinositol-4,5-bisphosphate Binding
        ],
        'fname': f'{OUT_DIR}/GSEA_MF_TNT_Factin_NES_barplot.png',
    },
    {
        'key':   'Extracellular_Vesicles',
        'ont':   'MF',
        'title': 'GO MF extracellular vesicle terms, OSD-564 (n = 5 vs 5)',
        'go_ids': [
            'GO:0035615',   # Clathrin-cargo Adaptor Activity
            'GO:0019905',   # Syntaxin Binding
            'GO:0031267',   # Small GTPase Binding
            'GO:0005085',   # Guanyl-nucleotide Exchange Factor Activity
            'GO:0097001',   # Ceramide Binding
            'GO:0005544',   # Calcium-dependent Phospholipid Binding
            'GO:0032266',   # Phosphatidylinositol-3-phosphate Binding
            'GO:0080025',   # Phosphatidylinositol-3,5-bisphosphate Binding
            'GO:0043130',   # Ubiquitin Binding
            'GO:0043395',   # Heparan Sulfate Proteoglycan Binding
        ],
        'fname': f'{OUT_DIR}/GSEA_MF_Extracellular_Vesicles_NES_barplot.png',
    },
    {
        'key':   'Amyloid_Tau',
        'ont':   'MF',
        'title': 'GO MF amyloid-β / tau terms, OSD-564 (n = 5 vs 5)',
        'go_ids': [
            'GO:0048156',   # Tau Protein Binding
            'GO:0050321',   # Tau-protein Kinase Activity
            'GO:0008017',   # Microtubule Binding
            'GO:0051721',   # Protein Phosphatase 2A Binding
            'GO:0004674',   # Protein Serine/Threonine Kinase Activity
            'GO:0004190',   # Aspartic-type Endopeptidase Activity
            'GO:0004222',   # Metalloendopeptidase Activity
            'GO:0043395',   # Heparan Sulfate Proteoglycan Binding
            'GO:0015485',   # Cholesterol Binding
            'GO:0070325',   # Lipoprotein Particle Receptor Binding
        ],
        'fname': f'{OUT_DIR}/GSEA_MF_Amyloid_Tau_NES_barplot.png',
    },
    {
        'key':   'BP_protrusion_exocytosis',
        'ont':   'BP',
        'title': 'GO BP protrusion and exocytosis terms, OSD-564 (n = 5 vs 5)',
        'names': [
            'Actin Filament Organization',
            'Axonogenesis',
            'Axon Guidance',
            'Neuron Projection Morphogenesis',
            'Regulated Exocytosis',
            'Synaptic Vesicle Recycling',
        ],
        'fname': f'{OUT_DIR}/GSEA_BP_protrusion_exocytosis_NES_barplot.png',
    },
    {
        'key':   'CC_compartments',
        'ont':   'CC',
        'title': 'GO CC compartment terms, OSD-564 (n = 5 vs 5)',
        'names': [
            'Recycling Endosome Membrane',
            'Endocytic Vesicle',
            'Clathrin-Coated Vesicle',
            'Exocytic Vesicle',
            'Mitochondrial Inner Membrane',
            'Ribosome',
        ],
        'fname': f'{OUT_DIR}/GSEA_CC_compartments_NES_barplot.png',
    },
]


def find_by_id(go_id, df):
    hits = df[df['Term'].str.contains(f'({go_id})', regex=False)]
    return hits.iloc[0] if len(hits) else None


def find_by_name(name, df):
    """Exact term-name match (case-insensitive) on 'Name (GO:xxxxxxx)'."""
    stem = df['Term'].str.replace(r'\s*\(GO:\d+\)\s*$', '', regex=True).str.lower()
    hits = df[stem == name.strip().lower()]
    return hits.iloc[0] if len(hits) else None


# ──────────────────────────────────────────────────────────────────
# Generate one figure per theme
# ──────────────────────────────────────────────────────────────────
for theme in THEMES:
    print(f"\nGenerating: {theme['key']}")
    df = GSEA[theme['ont']]
    records = []
    keys = theme.get('go_ids') or theme.get('names')
    for key in keys:
        row = find_by_id(key, df) if key.startswith('GO:') else find_by_name(key, df)
        if row is None:
            print(f"  [MISSING] {key}")
            continue
        records.append({
            'label': row['Term'],
            'nes':   float(row['NES']),
            'fdr':   float(row['FDR q-val']),
        })
        print(f"  {row['Term']:<62s} NES={row['NES']:+.3f}  FDR={row['FDR q-val']:.4f}")

    # listed order, first term at the top; 10-term MF panels are drawn taller
    # (90 x 78 mm) by nes_panel, 6/7-term panels at 90 x 60 mm
    nes_panel(records, theme['title'], os.path.splitext(theme['fname'])[0], xlim=TITLE_XLIM)

print("\nDone.")
