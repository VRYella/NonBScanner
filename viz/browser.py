"""
Genome Browser Track Generation
==============================

Generate tracks for IGV, JBrowse, and UCSC Genome Browser from NBDFinder results.
"""

import json
from typing import List, Dict, Any, Optional, Union
from pathlib import Path
import tempfile
import subprocess
from ..nbdio.writers import export_to_bed, export_to_gff3

class TrackGenerator:
    """
    Generate genome browser tracks from motif detection results.
    """
    
    def __init__(self):
        """Initialize track generator."""
        self.color_scheme = {
            'Curved_DNA': '#FF9AA2',      # Pink
            'Slipped_DNA': '#FFDABC',     # Light orange
            'Cruciform': '#E2F0CB',       # Light green
            'R-Loop': '#FFD3B6',          # Peach
            'Triplex': '#B5EAD7',         # Mint
            'G-Quadruplex': '#A2D7D8',    # Light blue
            'i-Motif': '#C7CEDB',         # Light purple
            'Z-DNA': '#FFAAA5',           # Coral
            'Hybrid': '#FF8B94',          # Rose
            'Cluster': '#B4A7D6'          # Lavender
        }
    
    def create_igv_session(self, motifs: List[Dict[str, Any]], 
                          genome_build: str = "hg38",
                          output_path: Optional[Path] = None) -> str:
        """
        Create IGV session XML for motif visualization.
        
        Args:
            motifs: List of motif detection results
            genome_build: Genome build identifier
            output_path: Output path for session file
            
        Returns:
            IGV session XML content
        """
        if output_path is None:
            output_path = Path("nbdfinder_session.xml")
        
        # Create BED track for each motif class
        class_motifs = {}
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            if motif_class not in class_motifs:
                class_motifs[motif_class] = []
            class_motifs[motif_class].append(motif)
        
        # Generate IGV session
        session_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<Session genome="{genome_build}" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
'''
        
        track_files = []
        for motif_class, class_motifs_list in class_motifs.items():
            # Create BED file for this class
            bed_content = export_to_bed(
                class_motifs_list, 
                score_type="normalized",
                include_subclass=True
            )
            
            bed_file = output_path.parent / f"{motif_class.lower()}_motifs.bed"
            with open(bed_file, 'w') as f:
                f.write(bed_content)
            
            track_files.append((bed_file, motif_class))
            
            session_xml += f'''        <Resource path="{bed_file.name}"/>
'''
        
        session_xml += '''    </Resources>
    <Panel height="800" name="DataPanel" width="1200">
'''
        
        # Add tracks
        for bed_file, motif_class in track_files:
            color = self.color_scheme.get(motif_class, '#000000')
            session_xml += f'''        <Track altColor="0,0,178" autoScale="false" color="{color}" 
               displayMode="COLLAPSED" featureVisibilityWindow="-1" 
               fontSize="10" id="{bed_file.name}" name="{motif_class} Motifs" 
               renderer="BASIC_FEATURE" sortable="false" visible="true" 
               windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" 
                      maximum="1.0" minimum="0.0" type="LINEAR"/>
        </Track>
'''
        
        session_xml += '''    </Panel>
</Session>'''
        
        # Write session file
        with open(output_path, 'w') as f:
            f.write(session_xml)
        
        return session_xml
    
    def create_jbrowse_config(self, motifs: List[Dict[str, Any]], 
                             output_dir: Path) -> Dict[str, Any]:
        """
        Create JBrowse 2 configuration for motif tracks.
        
        Args:
            motifs: List of motif detection results
            output_dir: Output directory for track files
            
        Returns:
            JBrowse configuration dictionary
        """
        output_dir.mkdir(exist_ok=True)
        
        # Group motifs by class
        class_motifs = {}
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            if motif_class not in class_motifs:
                class_motifs[motif_class] = []
            class_motifs[motif_class].append(motif)
        
        config = {
            "assemblies": [
                {
                    "name": "NBDFinder_Analysis",
                    "sequence": {
                        "type": "ReferenceSequenceTrack",
                        "trackId": "reference_sequence",
                        "adapter": {
                            "type": "FromConfigSequenceAdapter",
                            "features": []
                        }
                    }
                }
            ],
            "tracks": []
        }
        
        # Create tracks for each motif class
        for motif_class, class_motifs_list in class_motifs.items():
            # Export to GFF3
            gff_content = export_to_gff3(class_motifs_list)
            gff_file = output_dir / f"{motif_class.lower()}_motifs.gff3"
            
            with open(gff_file, 'w') as f:
                f.write(gff_content)
            
            # Create track configuration
            track_config = {
                "type": "FeatureTrack",
                "trackId": f"{motif_class.lower()}_track",
                "name": f"{motif_class} Motifs",
                "assemblyNames": ["NBDFinder_Analysis"],
                "category": ["NBDFinder"],
                "adapter": {
                    "type": "Gff3Adapter",
                    "gffLocation": {
                        "uri": str(gff_file.name)
                    }
                },
                "displays": [
                    {
                        "type": "LinearBasicDisplay",
                        "displayId": f"{motif_class.lower()}_display",
                        "renderer": {
                            "type": "SvgFeatureRenderer",
                            "color1": self.color_scheme.get(motif_class, '#000000')
                        }
                    }
                ]
            }
            
            config["tracks"].append(track_config)
        
        # Write configuration
        config_file = output_dir / "jbrowse_config.json"
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        return config
    
    def create_ucsc_track_hub(self, motifs: List[Dict[str, Any]], 
                             hub_name: str = "NBDFinder_Motifs",
                             output_dir: Path = None) -> Dict[str, str]:
        """
        Create UCSC Track Hub for motif visualization.
        
        Args:
            motifs: List of motif detection results
            hub_name: Name for the track hub
            output_dir: Output directory for hub files
            
        Returns:
            Dictionary of generated file paths
        """
        if output_dir is None:
            output_dir = Path("nbdfinder_track_hub")
        
        output_dir.mkdir(exist_ok=True)
        
        # Create hub.txt
        hub_content = f"""hub {hub_name}
shortLabel NBDFinder Non-B DNA Motifs
longLabel Non-B DNA Motif Detection Results from NBDFinder
genomesFile genomes.txt
email support@nbdfinder.org
descriptionUrl http://nbdfinder.org/hub_description.html
"""
        
        hub_file = output_dir / "hub.txt"
        with open(hub_file, 'w') as f:
            f.write(hub_content)
        
        # Create genomes.txt
        genomes_content = """genome hg38
trackDb hg38/trackDb.txt
"""
        
        genomes_file = output_dir / "genomes.txt"
        with open(genomes_file, 'w') as f:
            f.write(genomes_content)
        
        # Create track directory
        track_dir = output_dir / "hg38"
        track_dir.mkdir(exist_ok=True)
        
        # Group motifs by class and create BigBed files
        class_motifs = {}
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            if motif_class not in class_motifs:
                class_motifs[motif_class] = []
            class_motifs[motif_class].append(motif)
        
        trackdb_content = ""
        
        for motif_class, class_motifs_list in class_motifs.items():
            track_name = motif_class.lower().replace('-', '_')
            
            # Create BED file
            bed_content = export_to_bed(class_motifs_list, score_type="normalized")
            bed_file = track_dir / f"{track_name}.bed"
            
            with open(bed_file, 'w') as f:
                f.write(bed_content)
            
            # Add to trackDb
            color = self.color_scheme.get(motif_class, '#000000').replace('#', '')
            color_rgb = ','.join(str(int(color[i:i+2], 16)) for i in (0, 2, 4))
            
            trackdb_content += f"""
track {track_name}
shortLabel {motif_class} Motifs
longLabel {motif_class} Non-B DNA Motifs detected by NBDFinder
type bed 9
itemRgb on
color {color_rgb}
priority {len(trackdb_content.split('track')) * 10}
visibility dense

"""
        
        # Write trackDb.txt
        trackdb_file = track_dir / "trackDb.txt"
        with open(trackdb_file, 'w') as f:
            f.write(trackdb_content)
        
        return {
            'hub_file': str(hub_file),
            'genomes_file': str(genomes_file),
            'trackdb_file': str(trackdb_file),
            'track_directory': str(track_dir)
        }
    
    def create_circos_config(self, motifs: List[Dict[str, Any]], 
                           output_dir: Path) -> str:
        """
        Create Circos configuration for circular genome visualization.
        
        Args:
            motifs: List of motif detection results
            output_dir: Output directory for Circos files
            
        Returns:
            Path to main Circos configuration file
        """
        output_dir.mkdir(exist_ok=True)
        
        # Create data files for each motif class
        class_motifs = {}
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            if motif_class not in class_motifs:
                class_motifs[motif_class] = []
            class_motifs[motif_class].append(motif)
        
        # Main configuration
        circos_config = """
<colors>
<<include etc/colors.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<ideogram>
<spacing>
default = 0.005r
</spacing>

radius    = 0.90r
thickness = 20p
fill      = yes
</ideogram>

show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6

<tick>
spacing        = 5u
size           = 10p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>

<plots>
"""
        
        plot_radius = 0.85
        for motif_class, class_motifs_list in class_motifs.items():
            # Create data file
            data_content = ""
            for motif in class_motifs_list:
                chrom = motif.get('Chromosome', 'chr1')
                start = motif.get('Start', 0)
                end = motif.get('End', 0)
                score = motif.get('Normalized_Score', 0.5)
                
                data_content += f"{chrom} {start} {end} {score}\n"
            
            data_file = output_dir / f"{motif_class.lower()}_data.txt"
            with open(data_file, 'w') as f:
                f.write(data_content)
            
            # Add plot configuration
            color = self.color_scheme.get(motif_class, '#000000')
            circos_config += f"""
<plot>
type = histogram
file = {data_file.name}
r1   = {plot_radius:.2f}r
r0   = {plot_radius - 0.08:.2f}r
fill_color = {color}
</plot>
"""
            plot_radius -= 0.10
        
        circos_config += """
</plots>

karyotype = data/karyotype/karyotype.human.txt

<image>
<<include etc/image.conf>>
</image>

<<include etc/housekeeping.conf>>
"""
        
        # Write main config
        config_file = output_dir / "circos.conf"
        with open(config_file, 'w') as f:
            f.write(circos_config)
        
        return str(config_file)

__all__ = [
    'TrackGenerator'
]