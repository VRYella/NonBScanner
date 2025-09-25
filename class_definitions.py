"""
Class and Subclass Definitions for NBDFinder
=============================================

Centralized definition of all motif classes, subclasses, and their visual properties
for use in the selection interface.
"""

# Class definitions with comprehensive subclass information
CLASS_DEFINITIONS = {
    "Curved_DNA": {
        "id": "Curved_DNA",
        "name": "Curved DNA",
        "description": "A-tract mediated DNA curvature",
        "color": "#FF9AA2",
        "subclasses": [
            {"id": "global_curvature", "name": "Global Curvature", "description": "Widespread curvature patterns"},
            {"id": "local_curvature", "name": "Local Curvature", "description": "Localized bent regions"},
            {"id": "phased_a_tracts", "name": "Phased A-tracts", "description": "Periodic A-tract sequences"}
        ]
    },
    "Slipped_DNA": {
        "id": "Slipped_DNA", 
        "name": "Slipped DNA",
        "description": "Direct/tandem repeats forming slipped structures",
        "color": "#FFDAC1",
        "subclasses": [
            {"id": "direct_repeat", "name": "Direct Repeat", "description": "Direct tandem repeats"},
            {"id": "str", "name": "STR (Short Tandem Repeat)", "description": "Short tandem repeat sequences"},
            {"id": "slipped_str", "name": "Slipped DNA [STR]", "description": "STR-based slipped structures"}
        ]
    },
    "Cruciform": {
        "id": "Cruciform",
        "name": "Cruciform DNA",
        "description": "Inverted repeats forming four-way junctions",
        "color": "#E2F0CB",
        "subclasses": [
            {"id": "palindromic_ir", "name": "Palindromic Inverted Repeat", "description": "Perfect palindromic sequences"},
            {"id": "cruciform_hairpin", "name": "Cruciform Hairpin", "description": "Hairpin-forming inverted repeats"}
        ]
    },
    "R-Loop": {
        "id": "R-Loop",
        "name": "R-Loop",
        "description": "RNA-DNA hybrids with displaced ssDNA",
        "color": "#FFD3B6",
        "subclasses": [
            {"id": "r_loop_formation", "name": "R-Loop Formation Site", "description": "Sites prone to R-loop formation"},
            {"id": "stable_r_loop", "name": "Stable R-Loop", "description": "Thermodynamically stable R-loops"}
        ]
    },
    "Triplex": {
        "id": "Triplex",
        "name": "Triplex DNA",
        "description": "Three-stranded DNA structures",
        "color": "#B5EAD7",
        "subclasses": [
            {"id": "mirror_repeat", "name": "Mirror Repeat Triplex", "description": "Mirror repeat-based triplex"},
            {"id": "g_triplex", "name": "G-Triplex intermediate", "description": "G-rich triplex structures"},
            {"id": "sticky_dna", "name": "Sticky DNA", "description": "GA-rich triplex-forming sequences"}
        ]
    },
    "G-Quadruplex": {
        "id": "G-Quadruplex",
        "name": "G-Quadruplex Family", 
        "description": "G-quadruplex structures and variants",
        "color": "#A2D7D8",
        "subclasses": [
            {"id": "canonical_g4", "name": "Canonical G4", "description": "Standard G-quadruplex structures"},
            {"id": "bulged_g4", "name": "Bulged G4", "description": "G4 with bulge loops"},
            {"id": "relaxed_g4", "name": "Relaxed G4", "description": "Less stringent G4 patterns"},
            {"id": "bipartite_g4", "name": "Bipartite G4", "description": "Two-block G4 structures"},
            {"id": "multimeric_g4", "name": "Multimeric G4", "description": "Multiple G4 units"},
            {"id": "imperfect_g4", "name": "Imperfect G4", "description": "G4 with mismatches"},
            {"id": "g_triplex_int", "name": "G-Triplex", "description": "G-rich triplex intermediate"}
        ]
    },
    "i-Motif": {
        "id": "i-Motif",
        "name": "i-Motif",
        "description": "C-rich structures complementary to G4",
        "color": "#B0C4DE",
        "subclasses": [
            {"id": "canonical_i_motif", "name": "Canonical i-Motif", "description": "Standard i-motif structures"},
            {"id": "extended_i_motif", "name": "Extended i-Motif", "description": "Longer C-rich regions"},
            {"id": "ac_motif", "name": "AC-Motif", "description": "Alternating A-rich/C-rich regions"}
        ]
    },
    "Z-DNA": {
        "id": "Z-DNA",
        "name": "Z-DNA",
        "description": "Left-handed double helix",
        "color": "#FFB7B2",
        "subclasses": [
            {"id": "z_dna_classic", "name": "Z-DNA", "description": "Classic alternating purine-pyrimidine Z-DNA"},
            {"id": "egz", "name": "eGZ (Extruded-G Z-DNA)", "description": "Extruded-G Z-DNA variant"},
            {"id": "gc_rich_z", "name": "GC-rich Z-DNA", "description": "GC-rich Z-DNA forming regions"}
        ]
    },
    "A-philic_DNA": {
        "id": "A-philic_DNA",
        "name": "A-philic DNA",
        "description": "A-rich DNA with unique structural properties",
        "color": "#E6B8F7",
        "subclasses": [
            {"id": "a_philic", "name": "A-philic DNA", "description": "A-rich structural motifs"}
        ]
    },
    "Hybrid": {
        "id": "Hybrid",
        "name": "Hybrid Motifs",
        "description": "Overlapping/composite motifs",
        "color": "#C1A192",
        "subclasses": [
            {"id": "g4_triplex_overlap", "name": "G-Quadruplex_Triplex_DNA_Overlap", "description": "G4 and Triplex overlaps"},
            {"id": "slipped_z_overlap", "name": "Slipped_DNA_Z-DNA_Overlap", "description": "Slipped DNA and Z-DNA overlaps"},
            {"id": "multi_class_overlap", "name": "Multi-class Overlap", "description": "Multiple class overlaps"}
        ]
    },
    "Cluster": {
        "id": "Cluster",
        "name": "Non-B DNA Clusters",
        "description": "High-density motif regions",
        "color": "#A2C8CC",
        "subclasses": [
            {"id": "motif_hotspot", "name": "Motif Hotspot", "description": "High motif density regions"},
            {"id": "dense_cluster", "name": "Dense Cluster", "description": "Densely packed motifs"},
            {"id": "mixed_cluster", "name": "Mixed Cluster", "description": "Multiple motif types clustered"}
        ]
    }
}

# Helper function to get all class IDs
def get_all_class_ids():
    """Get list of all class IDs"""
    return list(CLASS_DEFINITIONS.keys())

# Helper function to get all subclass IDs for a class
def get_subclass_ids(class_id):
    """Get list of subclass IDs for a given class"""
    if class_id in CLASS_DEFINITIONS:
        return [sub["id"] for sub in CLASS_DEFINITIONS[class_id]["subclasses"]]
    return []

# Helper function to map selection to analysis format
def map_selection_to_analysis_format(selected_classes, selected_subclasses):
    """Convert UI selection to format expected by analysis functions"""
    # This will be used to filter which detectors to run
    analysis_classes = []
    
    for class_id in selected_classes:
        if class_id in CLASS_DEFINITIONS:
            analysis_classes.append(class_id)
    
    return analysis_classes

# Default selection (all classes selected)
DEFAULT_SELECTED_CLASSES = list(CLASS_DEFINITIONS.keys())
DEFAULT_SELECTED_SUBCLASSES = {}
for class_id, class_info in CLASS_DEFINITIONS.items():
    DEFAULT_SELECTED_SUBCLASSES[class_id] = [sub["id"] for sub in class_info["subclasses"]]