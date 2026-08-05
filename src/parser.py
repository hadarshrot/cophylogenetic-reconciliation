from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
import os
import json

def convert_tree_to_nwk(filepath, host_name, parasite_name):
    """
    Convert files in the format tree to the format nwk
    """
    host_edges = {}
    parasite_edges = {}
    host_names = {}
    parasite_names = {}
    phi_mapping = {}

    current_section = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line in ["HOSTTREE", "HOSTNAMES", "PARASITETREE", "PARASITENAMES", 
                        "PHI", "HOSTRANKS", "PARASITERANKS", "HOSTREGIONS", "REGIONCOSTS"]:
                current_section = line
                continue

            if current_section == "HOSTTREE":
                parts = line.split()
                host_edges[parts[0]] = (None if parts[1] == 'null' else parts[1], 
                                        None if parts[2] == 'null' else parts[2])
            elif current_section == "PARASITETREE":
                parts = line.split()
                parasite_edges[parts[0]] = (None if parts[1] == 'null' else parts[1], 
                                            None if parts[2] == 'null' else parts[2])
            elif current_section == "HOSTNAMES":
                parts = line.split('\t')
                if len(parts) >= 2: host_names[parts[0].strip()] = parts[1].strip()
            elif current_section == "PARASITENAMES":
                parts = line.split('\t')
                if len(parts) >= 2: parasite_names[parts[0].strip()] = parts[1].strip()
            elif current_section == "PHI":
                parts = line.split('\t')
                if len(parts) >= 2:
                    h_id = parts[0].strip()
                    p_ids = parts[1].strip().split()
                    h_name = host_names.get(h_id, h_id)
                    for p_id in p_ids:
                        p_name = parasite_names.get(p_id, p_id)
                        phi_mapping[p_name] = h_name

    def build_clade(node_id, edges, names):
        c1_id, c2_id = edges[node_id]
        node_name = names.get(node_id)
        
        if c1_id is None and c2_id is None:
            return Clade(name=node_name)
        
        clades = []
        if c1_id: clades.append(build_clade(c1_id, edges, names))
        if c2_id: clades.append(build_clade(c2_id, edges, names))
        return Clade(clades=clades, name=node_name if node_name != node_id else None)

    def find_root(edges):
        children = set()
        for c1, c2 in edges.values():
            if c1: children.add(c1)
            if c2: children.add(c2)
        return list(set(edges.keys()) - children)[0]

    host_root = build_clade(find_root(host_edges), host_edges, host_names)
    parasite_root = build_clade(find_root(parasite_edges), parasite_edges, parasite_names)
    
    host_tree = Tree(root=host_root, rooted=True)
    parasite_tree = Tree(root=parasite_root, rooted=True)

    output_dir = f"data/{host_name}_{parasite_name}"
    os.makedirs(output_dir, exist_ok=True)

    host_file_path = os.path.join(output_dir, f"{host_name}.nwk")
    parasite_file_path = os.path.join(output_dir, f"{parasite_name}.nwk")
    mapping_file_path = os.path.join(output_dir, f"mapping.json")

    Phylo.write(host_tree, host_file_path, "newick")
    Phylo.write(parasite_tree, parasite_file_path, "newick")

    with open(mapping_file_path, "w") as f:
        json.dump(phi_mapping, f, indent=4)

    print(f"Files have been saved")
    return host_tree, parasite_tree, phi_mapping
