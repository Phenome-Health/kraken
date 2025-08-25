"""
Metagraph generation utilities for Biolink knowledge graphs
Analyzes node categories, edge predicates, and connectivity patterns
"""

from pathlib import Path
from collections import defaultdict, Counter
import sys
from typing import Dict, Set, Tuple, Any, Iterator
import json
import logging

from .kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl


class MetagraphStats:
    """Container for metagraph statistics"""
    
    def __init__(self, source_name: str = "unknown"):
        self.source_name = source_name
        
        # Node statistics
        self.node_categories = Counter()  # category -> count
        self.total_nodes = 0
        
        # Edge statistics  
        self.edge_predicates = Counter()  # predicate -> count
        self.category_pairs = Counter()  # (subject_category, object_category) -> count
        self.predicate_category_pairs = Counter()  # (predicate, subject_cat, object_cat) -> count
        self.total_edges = 0
        
        # Connectivity statistics
        self.node_degrees = defaultdict(int)  # node_id -> degree
        self.category_connectivity = defaultdict(lambda: defaultdict(int))  # subj_cat -> obj_cat -> count
        
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary for JSON serialization"""
        return {
            'source': self.source_name,
            'summary': {
                'total_nodes': self.total_nodes,
                'total_edges': self.total_edges,
                'unique_node_categories': len(self.node_categories),
                'unique_edge_predicates': len(self.edge_predicates),
                'unique_category_pairs': len(self.category_pairs)
            },
            'node_categories': dict(self.node_categories.most_common()),
            'edge_predicates': dict(self.edge_predicates.most_common()),
            'category_pairs': {
                f"{subj}--{obj}": count 
                for (subj, obj), count in self.category_pairs.most_common()
            },
            'predicate_category_distribution': {
                predicate: {
                    f"{subj}--{obj}": count
                    for (pred, subj, obj), count in self.predicate_category_pairs.items()
                    if pred == predicate
                }
                for predicate in self.edge_predicates.keys()
            }
        }


def generate_metagraph_streaming(nodes_file: Path, edges_file: Path, source_name: str = None) -> MetagraphStats:
    """Generate metagraph statistics from JSONL files using streaming"""
    
    if source_name is None:
        source_name = nodes_file.parent.name
    
    logging.info(f"Generating metagraph for {source_name}")
    
    stats = MetagraphStats(source_name)
    
    # Phase 1: Analyze nodes and build category mapping
    categories_map = {}  # node_id -> categories
    
    logging.info("Analyzing nodes...")
    for node in stream_nodes_from_jsonl(nodes_file):
        node_id = node['id']
        categories = node['categories']
        for category in categories:
            categories_map[node_id] = category
            stats.node_categories[category] += 1
        stats.total_nodes += 1
        
        if stats.total_nodes % 1000000 == 0:
            logging.info(f"Processed {stats.total_nodes} nodes for metagraph")
    
    logging.info(f"Found {len(stats.node_categories)} unique node categories")

    if not categories_map:
        logging.error(f"Categories map is empty.")
        sys.exit(1)
    
    # Phase 2: Analyze edges
    logging.info("Analyzing edges...")
    for edge in stream_edges_from_jsonl(edges_file):
        subject_id = edge['subject']
        object_id = edge['object']
        predicate = edge['predicate']
        if subject_id in categories_map and object_id in categories_map:
            # Get categories (default to NamedThing if not found)
            subject_category = categories_map[subject_id]
            object_category = categories_map[object_id]
            
            # Update statistics
            stats.edge_predicates[predicate] += 1
            stats.category_pairs[(subject_category, object_category)] += 1
            stats.predicate_category_pairs[(predicate, subject_category, object_category)] += 1
            
            # Update node degrees
            stats.node_degrees[subject_id] += 1
            stats.node_degrees[object_id] += 1
            
            # Update category connectivity
            stats.category_connectivity[subject_category][object_category] += 1

            stats.total_edges += 1
        else:
            logging.warning(f"Orphan edge: Edge between {subject_id} and {object_id} is missing from categories map")
        
        if stats.total_edges % 2000000 == 0:
            logging.info(f"Processed {stats.total_edges} edges for metagraph")
    
    logging.info(f"Metagraph analysis complete: {stats.total_nodes} nodes, {stats.total_edges} edges")
    return stats


def save_metagraph(stats: MetagraphStats, output_file: Path):
    """Save metagraph statistics to JSON file"""
    logging.info(f"Saving metagraph to {output_file}")
    
    with open(output_file, 'w') as f:
        json.dump(stats.to_dict(), f, indent=2)
    
    logging.info(f"Metagraph saved: {stats.total_nodes} nodes, {stats.total_edges} edges")


def generate_metagraph_for_source(nodes_file: Path, edges_file: Path, output_dir: Path, source_name: str = None):
    """Generate and save complete metagraph suite for a single source"""
    if source_name is None:
        source_name = nodes_file.parent.name
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate core statistics
    stats = generate_metagraph_streaming(nodes_file, edges_file, source_name)
    
    generated_files = []
    
    # 1. Save main JSON statistics
    json_file = output_dir / f"{source_name}_metagraph.json"
    save_metagraph(stats, json_file)
    generated_files.append(json_file)
    
    # 2. Generate human-readable summary
    summary = generate_metagraph_summary(stats)
    summary_file = output_dir / f"{source_name}_metagraph_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(summary)
    generated_files.append(summary_file)
    logging.info(f"Summary saved: {summary_file}")
    
    # 3. Generate Cytoscape files with different thresholds
    cytoscape_files = []
    thresholds = [1, 5, 10]
    
    for threshold in thresholds:
        if threshold == 1:
            cyto_file = output_dir / f"{source_name}_cytoscape.json"
        else:
            cyto_file = output_dir / f"{source_name}_cytoscape_min{threshold}.json"
        
        create_cytoscape_metagraph(stats, cyto_file, min_edge_count=threshold)
        cytoscape_files.append(cyto_file)
        generated_files.append(cyto_file)
    
    # 4. Generate HTML viewer
    html_file = create_html_viewer(output_dir, cytoscape_files, source_name)
    generated_files.append(html_file)
    
    logging.info(f"Metagraph suite generated for {source_name}: {len(generated_files)} files")
    return generated_files


def compare_metagraphs(metagraph_files: list, output_file: Path):
    """Compare multiple metagraphs and generate comparison report"""
    logging.info(f"Comparing {len(metagraph_files)} metagraphs")
    
    metagraphs = []
    for file_path in metagraph_files:
        with open(file_path, 'r') as f:
            metagraphs.append(json.load(f))
    
    comparison = {
        'sources_compared': [mg['source'] for mg in metagraphs],
        'summary_comparison': {},
        'category_overlap': {},
        'predicate_overlap': {},
        'unique_to_source': {}
    }
    
    # Summary comparison
    for mg in metagraphs:
        source = mg['source']
        comparison['summary_comparison'][source] = mg['summary']
    
    # Find overlaps and unique elements
    all_categories = set()
    all_predicates = set()
    
    for mg in metagraphs:
        categories = set(mg['node_categories'].keys())
        predicates = set(mg['edge_predicates'].keys())
        
        all_categories.update(categories)
        all_predicates.update(predicates)
        
        source = mg['source']
        comparison['unique_to_source'][source] = {
            'categories': list(categories),
            'predicates': list(predicates)
        }
    
    # Calculate overlaps
    comparison['category_overlap'] = {
        'total_unique_categories': len(all_categories),
        'categories': list(all_categories)
    }
    
    comparison['predicate_overlap'] = {
        'total_unique_predicates': len(all_predicates), 
        'predicates': list(all_predicates)
    }
    
    # Save comparison
    with open(output_file, 'w') as f:
        json.dump(comparison, f, indent=2)
    
    logging.info(f"Metagraph comparison saved to {output_file}")


def generate_metagraph_summary(stats: MetagraphStats) -> str:
    """Generate human-readable summary of metagraph"""
    summary_lines = [
        f"=== Metagraph Summary: {stats.source_name} ===",
        f"Total Nodes: {stats.total_nodes:,}",
        f"Total Edges: {stats.total_edges:,}",
        f"Unique Node Categories: {len(stats.node_categories)}",
        f"Unique Edge Predicates: {len(stats.edge_predicates)}",
        "",
        "Top Node Categories:",
    ]
    
    for category, count in stats.node_categories.most_common(10):
        percentage = (count / stats.total_nodes) * 100
        summary_lines.append(f"  {category}: {count:,} ({percentage:.1f}%)")
    
    summary_lines.extend([
        "",
        "Top Edge Predicates:",
    ])
    
    for predicate, count in stats.edge_predicates.most_common(10):
        percentage = (count / stats.total_edges) * 100
        summary_lines.append(f"  {predicate}: {count:,} ({percentage:.1f}%)")
    
    summary_lines.extend([
        "",
        "Top Category Pairs:",
    ])
    
    for (subj_cat, obj_cat), count in stats.category_pairs.most_common(10):
        percentage = (count / stats.total_edges) * 100
        summary_lines.append(f"  {subj_cat} -> {obj_cat}: {count:,} ({percentage:.1f}%)")
    
    return "\n".join(summary_lines)


def create_cytoscape_metagraph(stats: MetagraphStats, output_file: Path, min_edge_count: int = 1):
    """Create Cytoscape-compatible metagraph visualization file"""
    logging.info(f"Creating Cytoscape metagraph for {stats.source_name}")
    
    # Create nodes (categories)
    nodes = []
    for category, count in stats.node_categories.items():
        nodes.append({
            'data': {
                'id': category,
                'label': category.replace('biolink:', ''),
                'node_count': count,
                'size': min(100, max(10, count // 1000))  # Scale node size
            }
        })
    
    # Create edges (category relationships)
    edges = []
    edge_id = 0
    for (source_cat, target_cat), count in stats.category_pairs.items():
        if count >= min_edge_count:
            edges.append({
                'data': {
                    'id': f'edge_{edge_id}',
                    'source': source_cat,
                    'target': target_cat,
                    'edge_count': count,
                    'weight': min(10, max(1, count // 1000))  # Scale edge weight
                }
            })
            edge_id += 1
    
    cytoscape_data = {
        'elements': {
            'nodes': nodes,
            'edges': edges
        },
        'metadata': {
            'source': stats.source_name,
            'total_nodes': stats.total_nodes,
            'total_edges': stats.total_edges,
            'min_edge_count_filter': min_edge_count
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(cytoscape_data, f, indent=2)
    
    logging.info(f"Cytoscape metagraph saved: {len(nodes)} category nodes, {len(edges)} category edges")


def create_html_viewer(output_dir: Path, metagraph_files: list, source_name: str = "metagraph"):
    """Create HTML viewer for interactive metagraph visualization"""
    logging.info(f"Creating HTML viewer for {source_name}")
    
    html_file = output_dir / f"{source_name}_viewer.html"
    
    # Find all Cytoscape files across the entire artifacts structure
    artifacts_root = Path("artifacts/metagraphs")
    all_cytoscape_files = []
    
    if artifacts_root.exists():
        # Find all cytoscape files in the entire metagraphs directory
        for cyto_file in artifacts_root.rglob("*_cytoscape*.json"):
            try:
                # Get relative path from current viewer location to the artifacts root
                # then navigate to the specific file
                rel_to_artifacts = cyto_file.relative_to(artifacts_root)
                
                # Calculate how many levels up to get to artifacts from output_dir
                try:
                    output_rel_to_artifacts = output_dir.relative_to(artifacts_root)
                    levels_up = len(output_rel_to_artifacts.parts)
                    up_path = "/".join([".."] * levels_up)
                    rel_path = f"{up_path}/{rel_to_artifacts}" if up_path else str(rel_to_artifacts)
                except ValueError:
                    # output_dir is not under artifacts_root, use absolute reference
                    rel_path = f"../../{rel_to_artifacts}"
                
                # Create descriptive name from path
                path_parts = cyto_file.relative_to(artifacts_root).parts
                if len(path_parts) >= 2:
                    section = path_parts[0]  # harmonized, unified, etc.
                    source = path_parts[1] if len(path_parts) > 2 else path_parts[-1].split('_')[0]
                    
                    # Handle different threshold files
                    filename = cyto_file.stem
                    if 'min' in filename:
                        threshold = filename.split('min')[-1]
                        display_name = f"{section.title()} - {source.title()} (Min {threshold} edges)"
                    else:
                        display_name = f"{section.title()} - {source.title()}"
                else:
                    display_name = cyto_file.stem.replace('_cytoscape', '').replace('_', ' ').title()
                
                all_cytoscape_files.append({
                    'path': str(rel_path),
                    'name': display_name,
                    'section': section if len(path_parts) >= 2 else 'other',
                    'source': source if len(path_parts) >= 2 else 'unknown'
                })
            except ValueError:
                # Skip files that can't be made relative to output_dir
                continue
    
    # Also include local files (fallback)
    for file_path in metagraph_files:
        if file_path.name.endswith('_cytoscape.json'):
            try:
                rel_path = file_path.relative_to(output_dir)
                display_name = f"Local - {file_path.stem.replace('_cytoscape', '').replace('_', ' ').title()}"
                
                # Only add if not already in the global list
                if not any(opt['path'] == str(rel_path) for opt in all_cytoscape_files):
                    all_cytoscape_files.append({
                        'path': str(rel_path),
                        'name': display_name,
                        'section': 'local',
                        'source': source_name
                    })
            except ValueError:
                continue
    
    # Sort by section and then by name
    all_cytoscape_files.sort(key=lambda x: (x['section'], x['name']))
    
    file_options = all_cytoscape_files
    
    def build_dropdown_options(options):
        if not options:
            return ""
        
        # Group by section
        sections = {}
        for opt in options:
            section = opt['section']
            if section not in sections:
                sections[section] = []
            sections[section].append(opt)
        
        # Build HTML with optgroups
        html_parts = []
        section_labels = {
            'harmonized': 'Source Graphs',
            'unified': 'Integrated Graph',
            'local': 'Current Graph',
            'other': 'Other'
        }
        
        for section in ['harmonized', 'unified', 'local', 'other']:
            if section in sections and sections[section]:
                label = section_labels.get(section, section.title())
                html_parts.append(f'<optgroup label="{label}">')
                for opt in sections[section]:
                    html_parts.append(f'  <option value="{opt["path"]}">{opt["name"]}</option>')
                html_parts.append('</optgroup>')
        
        return chr(10).join(html_parts)
    
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Phenome-KG Metagraph Viewer - {source_name.title()}</title>
    <script src="https://unpkg.com/cytoscape@3.21.0/dist/cytoscape.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f8f9fa;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            text-align: center;
        }}
        
        .controls {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            align-items: center;
        }}
        
        #cy {{
            width: 100%;
            height: 70vh;
            border: 1px solid #ddd;
            background-color: white;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        button {{
            padding: 10px 16px;
            background: #4CAF50;
            color: white;
            border: none;
            cursor: pointer;
            border-radius: 6px;
            font-size: 14px;
            transition: background-color 0.3s;
        }}
        
        button:hover {{
            background: #45a049;
        }}
        
        .secondary-btn {{
            background: #2196F3;
        }}
        
        .secondary-btn:hover {{
            background: #1976D2;
        }}
        
        select, input[type="file"] {{
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 6px;
            font-size: 14px;
        }}
        
        .info {{
            background: white;
            margin-top: 20px;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 15px 0;
        }}
        
        .stat-card {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            text-align: center;
        }}
        
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: #667eea;
        }}
        
        .file-selector {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        
        .instructions {{
            background: #e3f2fd;
            border-left: 4px solid #2196F3;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 6px 6px 0;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Phenome-KG Metagraph Visualization</h1>
        <p>Interactive network view of knowledge graph structure - {source_name.title()}</p>
    </div>
    
    <div class="controls">
        <div class="file-selector">
            <label for="file-select">Metagraph:</label>
            <select id="file-select" onchange="loadSelectedFile()">
                <option value="">Select a metagraph...</option>
                {build_dropdown_options(file_options)}
            </select>
        </div>
        
        <div>
            <label for="layout-select">Layout:</label>
            <select id="layout-select">
                <option value="cose">Force-directed (COSE)</option>
                <option value="circle">Circle</option>
                <option value="grid">Grid</option>
                <option value="breadthfirst">Hierarchical</option>
                <option value="concentric">Concentric</option>
                <option value="cola">CoLa</option>
            </select>
            <button onclick="applyLayout()">Apply</button>
        </div>
        
        <div>
            <button onclick="fitGraph()" class="secondary-btn">Fit to Screen</button>
            <button onclick="resetZoom()" class="secondary-btn">Reset Zoom</button>
            <button onclick="exportImage()" class="secondary-btn">Export PNG</button>
        </div>
        
        <div>
            <input type="file" id="file-input" accept=".json" onchange="loadCustomFile(event)" style="display: none;">
            <button onclick="document.getElementById('file-input').click()" class="secondary-btn">Load Custom File</button>
        </div>
    </div>
    
    <div id="cy"></div>
    
    <div class="info">
        <div class="instructions">
            <h3>How to Use:</h3>
            <ul>
                <li><strong>Pan:</strong> Click and drag on empty space</li>
                <li><strong>Zoom:</strong> Mouse wheel or pinch gesture</li>
                <li><strong>Node Details:</strong> Click on nodes to see category information</li>
                <li><strong>Edge Details:</strong> Click on edges to see relationship counts</li>
                <li><strong>Node Sizes:</strong> Proportional to number of entities in each category</li>
                <li><strong>Edge Thickness:</strong> Proportional to number of connections between categories</li>
            </ul>
        </div>
        
        <div class="stats" id="stats-container">
            <div class="stat-card">
                <div class="stat-value" id="node-count">-</div>
                <div>Categories</div>
            </div>
            <div class="stat-card">
                <div class="stat-value" id="edge-count">-</div>
                <div>Connections</div>
            </div>
            <div class="stat-card">
                <div class="stat-value" id="total-nodes">-</div>
                <div>Total Entities</div>
            </div>
            <div class="stat-card">
                <div class="stat-value" id="total-edges">-</div>
                <div>Total Relations</div>
            </div>
        </div>
        
        <p id="status-message">Select a metagraph file to begin visualization...</p>
    </div>

    <script>
        let cy;
        let currentData = null;
        
        function initializeCytoscape(data) {{
            if (cy) {{
                cy.destroy();
            }}
            
            currentData = data;
            
            cy = cytoscape({{
                container: document.getElementById('cy'),
                elements: data.elements,
                style: [
                    {{
                        selector: 'node',
                        style: {{
                            'background-color': function(ele) {{
                                const category = ele.data('id') || '';
                                const count = ele.data('node_count') || 1;
                                
                                // Determine theme based on category
                                let baseColor = [200, 200, 200]; // Default gray
                                
                                // Genomic/Genetic entities (blues)
                                if (category.includes('Gene') || category.includes('Protein') || 
                                    category.includes('Transcript') || category.includes('MicroRNA') ||
                                    category.includes('RNA') || category.includes('Genomic') ||
                                    category.includes('Polypeptide') || category.includes('NucleicAcidEntity') ||
                                    category.includes('SequenceVariant') || category.includes('Haplotype') ||
                                    category.includes('MacromolecularComplex'))  {{
                                    baseColor = [173, 216, 230]; // Light blue
                                }}
                                // Chemical/Drug entities (greens) 
                                else if (category.includes('Chemical') || category.includes('Drug') ||
                                         category.includes('SmallMolecule') || category.includes('Compound') ||
                                         category.includes('MolecularMixture') || category.includes('Metabolite') ||
                                         category.includes('Food') || category.includes('MolecularEntity')) {{
                                    baseColor = [180, 215, 180]; // Softer sage green
                                }}
                                // Disease/Phenotype entities (reds/pinks)
                                else if (category.includes('Disease') || category.includes('Phenotypic') ||
                                         category.includes('Symptom') || category.includes('ClinicalFinding') ||
                                         category.includes('BehavioralFeature')) {{
                                    baseColor = [255, 182, 193]; // Light pink
                                }}
                                // Anatomy/Biology entities (purples)
                                else if (category.includes('Anatomical') || category.includes('Cell') ||
                                         category.includes('Tissue') || category.includes('Organ') ||
                                         category.includes('OrganismTaxon') || category.includes('Cellular')) {{
                                    baseColor = [221, 160, 221]; // Plum
                                }}
                                // Pathway/Process entities (oranges)
                                else if (category.includes('Pathway') || category.includes('Process') ||
                                         category.includes('Activity') || category.includes('Function') ||
                                         category.includes('Event') || category.includes('BiologicalEntity')) {{
                                    baseColor = [255, 218, 185]; // Peach
                                }}
                                
                                // Adjust intensity based on count (darker for more entities)
                                const intensity = Math.min(0.4, Math.max(0.1, count / 10000));
                                const r = Math.floor(baseColor[0] - (baseColor[0] - 255) * intensity);
                                const g = Math.floor(baseColor[1] - (baseColor[1] - 255) * intensity);
                                const b = Math.floor(baseColor[2] - (baseColor[2] - 255) * intensity);
                                
                                return `rgb(${{r}}, ${{g}}, ${{b}})`;
                            }},
                            'label': 'data(label)',
                            'color': '#000',
                            'text-valign': 'center',
                            'text-halign': 'center',
                            'font-size': function(ele) {{
                                const count = ele.data('node_count') || 1;
                                return Math.max(14, Math.min(20, 12 + count / 50)) + 'px';
                            }},
                            'font-weight': 'bold',
                            'text-outline-width': 1,
                            'text-outline-color': '#fff',
                            'width': function(ele) {{
                                const count = ele.data('node_count') || 1;
                                return Math.max(40, Math.min(120, 30 + count / 3));
                            }},
                            'height': function(ele) {{
                                const count = ele.data('node_count') || 1;
                                return Math.max(40, Math.min(120, 30 + count / 3));
                            }}
                        }}
                    }},
                    {{
                        selector: 'edge',
                        style: {{
                            'width': function(ele) {{
                                const count = ele.data('edge_count') || 1;
                                return Math.max(1, Math.min(8, count / 8));
                            }},
                            'line-color': '#888',
                            'target-arrow-color': '#888',
                            'target-arrow-shape': 'triangle',
                            'target-arrow-size': function(ele) {{
                                const count = ele.data('edge_count') || 1;
                                return Math.max(6, Math.min(12, count / 10));
                            }},
                            'curve-style': 'bezier',
                            'opacity': 0.8
                        }}
                    }},
                    {{
                        selector: 'node:selected',
                        style: {{
                            'border-width': 3,
                            'border-color': '#FF4444'
                        }}
                    }},
                    {{
                        selector: 'edge:selected',
                        style: {{
                            'line-color': '#FF4444',
                            'target-arrow-color': '#FF4444',
                            'opacity': 1
                        }}
                    }}
                ],
                layout: {{
                    name: 'cose',
                    idealEdgeLength: 120,
                    nodeOverlap: 30,
                    refresh: 20,
                    fit: true,
                    padding: 50,
                    randomize: false,
                    componentSpacing: 150,
                    nodeRepulsion: 800000,
                    edgeElasticity: 200,
                    nestingFactor: 5,
                    gravity: 100,
                    numIter: 2000,
                    initialTemp: 300,
                    coolingFactor: 0.95,
                    minTemp: 1.0
                }}
            }});
            
            // Add interaction events
            cy.on('tap', 'node', function(evt) {{
                const node = evt.target;
                const info = {{
                    category: node.data('label'),
                    entity_count: node.data('node_count'),
                    id: node.data('id')
                }};
                
                alert(`Category: ${{info.category}}\\nEntities: ${{info.entity_count.toLocaleString()}}\\nID: ${{info.id}}`);
                console.log('Node details:', info);
            }});
            
            cy.on('tap', 'edge', function(evt) {{
                const edge = evt.target;
                const sourceLabel = cy.getElementById(edge.data('source')).data('label');
                const targetLabel = cy.getElementById(edge.data('target')).data('label');
                
                const info = {{
                    connection: `${{sourceLabel}} → ${{targetLabel}}`,
                    edge_count: edge.data('edge_count')
                }};
                
                alert(`Connection: ${{info.connection}}\\nRelations: ${{info.edge_count.toLocaleString()}}`);
                console.log('Edge details:', info);
            }});
            
            // Update statistics
            updateStats(data);
        }}
        
        function updateStats(data) {{
            const nodeCount = data.elements.nodes.length;
            const edgeCount = data.elements.edges.length;
            const totalNodes = data.elements.nodes.reduce((sum, n) => sum + (n.data.node_count || 0), 0);
            const totalEdges = data.elements.edges.reduce((sum, e) => sum + (e.data.edge_count || 0), 0);
            
            document.getElementById('node-count').textContent = nodeCount.toLocaleString();
            document.getElementById('edge-count').textContent = edgeCount.toLocaleString();
            document.getElementById('total-nodes').textContent = totalNodes.toLocaleString();
            document.getElementById('total-edges').textContent = totalEdges.toLocaleString();
        }}
        
        function applyLayout() {{
            if (!cy) return;
            const layoutName = document.getElementById('layout-select').value;
            
            let layoutOptions = {{
                name: layoutName,
                fit: true,
                padding: 50,
                animate: true,
                animationDuration: 1000
            }};
            
            // Use the same detailed configuration for COSE as the initial layout
            if (layoutName === 'cose') {{
                layoutOptions = {{
                    name: 'cose',
                    idealEdgeLength: 120,
                    nodeOverlap: 30,
                    refresh: 20,
                    fit: true,
                    padding: 50,
                    randomize: false,
                    componentSpacing: 150,
                    nodeRepulsion: 800000,
                    edgeElasticity: 200,
                    nestingFactor: 5,
                    gravity: 100,
                    numIter: 2000,
                    initialTemp: 300,
                    coolingFactor: 0.95,
                    minTemp: 1.0,
                    animate: true,
                    animationDuration: 1000
                }};
            }}
            
            const layout = cy.layout(layoutOptions);
            layout.run();
        }}
        
        function fitGraph() {{
            if (cy) cy.fit(null, 50);
        }}
        
        function resetZoom() {{
            if (cy) {{
                cy.zoom(1);
                cy.center();
            }}
        }}
        
        function exportImage() {{
            if (cy) {{
                const png64 = cy.png({{scale: 2, full: true}});
                const link = document.createElement('a');
                link.href = png64;
                link.download = 'metagraph.png';
                link.click();
            }}
        }}
        
        function loadSelectedFile() {{
            const select = document.getElementById('file-select');
            const filePath = select.value;
            
            if (!filePath) return;
            
            document.getElementById('status-message').textContent = `Loading ${{select.options[select.selectedIndex].text}}...`;
            
            fetch(filePath)
                .then(response => {{
                    if (!response.ok) throw new Error(`HTTP ${{response.status}}`);
                    return response.json();
                }})
                .then(data => {{
                    initializeCytoscape(data);
                    document.getElementById('status-message').textContent = 
                        `Loaded: ${{select.options[select.selectedIndex].text}}`;
                }})
                .catch(error => {{
                    console.error('Error loading file:', error);
                    document.getElementById('status-message').textContent = 
                        `Error loading file: ${{error.message}}`;
                }});
        }}
        
        function loadCustomFile(event) {{
            const file = event.target.files[0];
            if (!file) return;
            
            const reader = new FileReader();
            reader.onload = function(e) {{
                try {{
                    const data = JSON.parse(e.target.result);
                    initializeCytoscape(data);
                    document.getElementById('status-message').textContent = `Loaded custom file: ${{file.name}}`;
                }} catch (error) {{
                    alert('Error loading file: ' + error.message);
                }}
            }};
            reader.readAsText(file);
        }}
        
        // Auto-load the current graph's default file
        window.addEventListener('DOMContentLoaded', function() {{
            const select = document.getElementById('file-select');
            
            // First try to find the "Local" option (current graph)
            let defaultOption = null;
            for (let i = 0; i < select.options.length; i++) {{
                const option = select.options[i];
                if (option.text.includes('Local - {source_name.title()}')) {{
                    defaultOption = option;
                    break;
                }}
            }}
            
            // If no local option found, try to find current source in the list
            if (!defaultOption) {{
                for (let i = 0; i < select.options.length; i++) {{
                    const option = select.options[i];
                    if (option.text.toLowerCase().includes('{source_name.lower()}') && !option.text.includes('Min')) {{
                        defaultOption = option;
                        break;
                    }}
                }}
            }}
            
            // Fallback to first available option
            if (!defaultOption && select.options.length > 1) {{
                defaultOption = select.options[1];
            }}
            
            if (defaultOption) {{
                select.value = defaultOption.value;
                loadSelectedFile();
            }}
        }});
    </script>
</body>
</html>'''
    
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logging.info(f"HTML viewer created: {html_file}")
    return html_file