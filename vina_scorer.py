#!/usr/bin/env python3
"""
Standalone Vina-like Scoring Functions for Protein-DNA Complexes

This version only depends on BioPython and NumPy, avoiding MDAnalysis
and SciPy compatibility issues.

Author: Manus AI Assistant and Xiaopeng Xu
Date: 2025-07-09
"""
import numpy as np
import sys
import os
from pathlib import Path
from tqdm import tqdm
from typing import Dict, List, Tuple, Optional, Any

try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


class AtomProperties:
    """Atom properties and force field parameters"""
    
    # Van der Waals radii (Angstroms)
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80,
        'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98, 'Mg': 1.73, 'Ca': 2.31,
        'Mn': 1.73, 'Fe': 1.72, 'Zn': 1.39, 'Se': 1.90
    }
    
    # AMBER-based partial charges for DNA nucleotides
    DNA_CHARGES = {
        'DA': {  # Adenine
            'P': 1.17, 'O1P': -0.78, 'O2P': -0.78, 'O5\'': -0.50,
            'C5\'': 0.06, 'C4\'': 0.11, 'O4\'': -0.35, 'C1\'': 0.04,
            'N9': -0.03, 'C8': 0.16, 'N7': -0.62, 'C5': 0.07,
            'C6': 0.69, 'N6': -0.91, 'N1': -0.76, 'C2': 0.57,
            'N3': -0.74, 'C4': 0.38, 'C2\'': -0.09, 'C3\'': 0.07, 'O3\'': -0.65
        },
        'DT': {  # Thymine
            'P': 1.17, 'O1P': -0.78, 'O2P': -0.78, 'O5\'': -0.50,
            'C5\'': 0.06, 'C4\'': 0.11, 'O4\'': -0.35, 'C1\'': 0.07,
            'N1': -0.02, 'C6': -0.22, 'C2': 0.52, 'O2': -0.43,
            'N3': -0.43, 'C4': 0.60, 'O4': -0.58, 'C5': -0.16,
            'C7': -0.23, 'C2\'': -0.09, 'C3\'': 0.07, 'O3\'': -0.65
        },
        'DG': {  # Guanine
            'P': 1.17, 'O1P': -0.78, 'O2P': -0.78, 'O5\'': -0.50,
            'C5\'': 0.06, 'C4\'': 0.11, 'O4\'': -0.35, 'C1\'': 0.02,
            'N9': 0.05, 'C8': 0.14, 'N7': -0.57, 'C5': 0.20,
            'C6': 0.49, 'O6': -0.57, 'N1': -0.51, 'C2': 0.74,
            'N2': -0.92, 'N3': -0.66, 'C4': 0.18, 'C2\'': -0.09,
            'C3\'': 0.07, 'O3\'': -0.65
        },
        'DC': {  # Cytosine
            'P': 1.17, 'O1P': -0.78, 'O2P': -0.78, 'O5\'': -0.50,
            'C5\'': 0.06, 'C4\'': 0.11, 'O4\'': -0.35, 'C1\'': 0.01,
            'N1': -0.05, 'C6': -0.02, 'C2': 0.60, 'O2': -0.65,
            'N3': -0.76, 'C4': 0.82, 'N4': -0.95, 'C5': -0.52,
            'C2\'': -0.09, 'C3\'': 0.07, 'O3\'': -0.65
        }
    }
    
    # Protein residue charges (key atoms)
    PROTEIN_CHARGES = {
        'ARG': {'NH1': -0.86, 'NH2': -0.86, 'NE': -0.70, 'CZ': 0.64},
        'LYS': {'NZ': -0.30, 'CE': 0.25},
        'ASP': {'OD1': -0.80, 'OD2': -0.80, 'CG': 0.75},
        'GLU': {'OE1': -0.80, 'OE2': -0.80, 'CD': 0.75},
        'HIS': {'ND1': -0.40, 'NE2': -0.40, 'CE1': 0.25, 'CD2': 0.13},
        'SER': {'OG': -0.65, 'CB': 0.21},
        'THR': {'OG1': -0.65, 'CB': 0.30},
        'TYR': {'OH': -0.54, 'CZ': 0.32},
        'CYS': {'SG': -0.23, 'CB': 0.19},
        'ASN': {'OD1': -0.61, 'ND2': -0.92, 'CG': 0.71},
        'GLN': {'OE1': -0.61, 'NE2': -0.92, 'CD': 0.71}
    }
    
    @classmethod
    def get_vdw_radius(cls, element: str) -> float:
        """Get van der Waals radius for element"""
        return cls.VDW_RADII.get(element.upper(), 1.70)
    
    @classmethod
    def get_partial_charge(cls, atom: Atom, residue_name: str) -> float:
        """Get partial charge for atom"""
        atom_name = atom.get_name()
        
        # DNA charges
        if residue_name in cls.DNA_CHARGES:
            return cls.DNA_CHARGES[residue_name].get(atom_name, 0.0)
        
        # Protein charges
        if residue_name in cls.PROTEIN_CHARGES:
            return cls.PROTEIN_CHARGES[residue_name].get(atom_name, 0.0)
        
        # Default charges by element
        element = atom.element.upper()
        if element == 'N':
            return -0.30
        elif element == 'O':
            return -0.40
        elif element == 'P':
            return 1.17
        else:
            return 0.0
    
    @classmethod
    def is_hydrophobic(cls, atom: Atom, residue_name: str) -> bool:
        """Check if atom is hydrophobic"""
        if atom.element.upper() == 'C':
            # More sophisticated logic could be added
            hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET'}
            return residue_name in hydrophobic_residues
        return False
    
    @classmethod
    def can_hydrogen_bond(cls, atom: Atom) -> bool:
        """Check if atom can participate in hydrogen bonding"""
        element = atom.element.upper()
        return element in {'N', 'O', 'F'}


class StandaloneVinaScorer:
    """
    Standalone Vina-like scoring function for protein-DNA complexes
    
    This implementation avoids MDAnalysis and SciPy dependencies,
    using only BioPython and NumPy for maximum compatibility.
    """
    
    def __init__(self, 
                 cutoff_distance: float = 8.0,
                 include_electrostatics: bool = True,
                 dielectric_constant: float = 4.0,
                 include_dna_specific: bool = True):
        """
        Initialize the scorer
        
        Args:
            cutoff_distance: Maximum interaction distance (Angstroms)
            include_electrostatics: Include electrostatic terms
            dielectric_constant: Dielectric constant for electrostatics
            include_dna_specific: Include DNA-specific interaction terms
        """
        self.cutoff_distance = cutoff_distance
        self.include_electrostatics = include_electrostatics
        self.dielectric_constant = dielectric_constant
        self.include_dna_specific = include_dna_specific
        
        # Vina scoring function weights
        self.weights = {
            'gauss1': -0.0356,
            'gauss2': -0.00516,
            'repulsion': 0.840,
            'hydrophobic': -0.0351,
            'hydrogen_bonding': -0.587,
            'electrostatic': -0.2,
            'major_groove': -0.1,
            'minor_groove': -0.05,
            'base_stacking': -0.08,
            'phosphate_contact': -0.15,
            'n_rot': 0.0585
        }
        
        # Scoring function parameters
        self.gauss1_width = 0.5
        self.gauss2_offset = 3.0
        self.gauss2_width = 2.0
        self.hydrophobic_cutoff = (0.5, 1.5)
        self.hb_cutoff = (-0.7, 0.0)
        self.electrostatic_constant = 332.0636
    
    def load_structure(self, structure_file: str) -> Optional[Structure]:
        """Load structure from PDB or CIF file"""
        structure_path = Path(structure_file)
        if not structure_path.exists():
            raise FileNotFoundError(f"Structure file not found: {structure_file}")
        
        try:
            if structure_path.suffix.lower() == '.pdb':
                parser = PDBParser(QUIET=True)
            elif structure_path.suffix.lower() in ['.cif', '.mmcif']:
                parser = MMCIFParser(QUIET=True)
            else:
                raise ValueError(f"Unsupported file format: {structure_path.suffix}")
            
            structure = parser.get_structure('complex', str(structure_file))
            return structure
            
        except Exception as e:
            print(f"Error loading structure: {e}")
            return None
    
    def identify_chains(self, structure: Structure) -> Dict[str, List[str]]:
        """Identify protein and DNA chains"""
        protein_chains = []
        dna_chains = []
        
        for model in structure:
            for chain in model:
                residue_names = [residue.get_resname() for residue in chain 
                               if residue.id[0] == ' ']  # Standard residues only
                
                if not residue_names:
                    continue
                
                # DNA nucleotides
                dna_residues = {'DA', 'DT', 'DG', 'DC', 'A', 'T', 'G', 'C', 'DU', 'U'}
                
                dna_count = sum(1 for name in residue_names if name in dna_residues)
                total_residues = len(residue_names)
                
                if total_residues > 0:
                    dna_fraction = dna_count / total_residues
                    
                    if dna_fraction > 0.5:
                        dna_chains.append(chain.id)
                    else:
                        protein_chains.append(chain.id)
        
        return {'protein': protein_chains, 'dna': dna_chains}
    
    def extract_atoms(self, structure: Structure, chain_ids: List[str]) -> List[Tuple[Atom, str, str]]:
        """Extract atoms from specified chains"""
        atoms = []
        
        for model in structure:
            for chain in model:
                if chain.id in chain_ids:
                    for residue in chain:
                        if residue.id[0] == ' ':  # Standard residues only
                            for atom in residue:
                                if atom.element != 'H':  # Skip hydrogens
                                    atoms.append((atom, residue.get_resname(), chain.id))
        
        return atoms
    
    def calculate_surface_distance(self, atom1: Atom, atom2: Atom) -> float:
        """Calculate surface distance between atoms"""
        coord1 = atom1.get_coord()
        coord2 = atom2.get_coord()
        
        center_distance = np.linalg.norm(coord1 - coord2)
        
        vdw1 = AtomProperties.get_vdw_radius(atom1.element)
        vdw2 = AtomProperties.get_vdw_radius(atom2.element)
        
        return center_distance - vdw1 - vdw2
    
    def calculate_vina_terms(self, surface_distance: float) -> Dict[str, float]:
        """Calculate standard Vina scoring terms"""
        return {
            'gauss1': np.exp(-(surface_distance / self.gauss1_width) ** 2),
            'gauss2': np.exp(-((surface_distance - self.gauss2_offset) / self.gauss2_width) ** 2),
            'repulsion': surface_distance ** 2 if surface_distance < 0 else 0.0
        }
    
    def calculate_hydrophobic_term(self, surface_distance: float, 
                                  atom1: Atom, atom2: Atom,
                                  res1: str, res2: str) -> float:
        """Calculate hydrophobic interaction term"""
        if (AtomProperties.is_hydrophobic(atom1, res1) and 
            AtomProperties.is_hydrophobic(atom2, res2)):
            
            if surface_distance < self.hydrophobic_cutoff[0]:
                return 1.0
            elif surface_distance > self.hydrophobic_cutoff[1]:
                return 0.0
            else:
                return (self.hydrophobic_cutoff[1] - surface_distance) / (
                    self.hydrophobic_cutoff[1] - self.hydrophobic_cutoff[0])
        
        return 0.0
    
    def calculate_hydrogen_bonding(self, surface_distance: float,
                                  atom1: Atom, atom2: Atom) -> float:
        """Calculate hydrogen bonding term"""
        if (AtomProperties.can_hydrogen_bond(atom1) and 
            AtomProperties.can_hydrogen_bond(atom2)):
            
            if surface_distance < self.hb_cutoff[0]:
                return 1.0
            elif surface_distance > self.hb_cutoff[1]:
                return 0.0
            else:
                return (self.hb_cutoff[1] - surface_distance) / (
                    self.hb_cutoff[1] - self.hb_cutoff[0])
        
        return 0.0
    
    def calculate_electrostatic_term(self, atom1: Atom, atom2: Atom,
                                   res1_name: str, res2_name: str,
                                   center_distance: float) -> float:
        """Calculate electrostatic interaction"""
        if not self.include_electrostatics:
            return 0.0
        
        charge1 = AtomProperties.get_partial_charge(atom1, res1_name)
        charge2 = AtomProperties.get_partial_charge(atom2, res2_name)
        
        if abs(charge1) < 0.01 or abs(charge2) < 0.01:
            return 0.0
        
        # Distance-dependent dielectric
        effective_dielectric = max(1.0, self.dielectric_constant * center_distance / 5.0)
        
        return (self.electrostatic_constant * charge1 * charge2) / (
            effective_dielectric * center_distance)
    
    def is_in_major_groove(self, protein_atom: Atom, dna_residue: Residue) -> bool:
        """Simple heuristic for major groove detection"""
        protein_coord = protein_atom.get_coord()
        
        # Get base atoms
        base_atoms = []
        for atom in dna_residue:
            if atom.get_name() in ['N1', 'N3', 'N7', 'N9', 'C2', 'C4', 'C6', 'C8']:
                base_atoms.append(atom)
        
        if len(base_atoms) < 2:
            return False
        
        base_center = np.mean([atom.get_coord() for atom in base_atoms], axis=0)
        distance = np.linalg.norm(protein_coord - base_center)
        
        # Major groove: intermediate distance from base center
        return 4.0 <= distance <= 8.0
    
    def is_in_minor_groove(self, protein_atom: Atom, dna_residue: Residue) -> bool:
        """Simple heuristic for minor groove detection"""
        protein_coord = protein_atom.get_coord()
        
        # Get sugar atoms
        sugar_atoms = []
        for atom in dna_residue:
            if atom.get_name() in ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'O4\'']:
                sugar_atoms.append(atom)
        
        if len(sugar_atoms) < 2:
            return False
        
        sugar_center = np.mean([atom.get_coord() for atom in sugar_atoms], axis=0)
        distance = np.linalg.norm(protein_coord - sugar_center)
        
        # Minor groove: closer to sugar backbone
        return 3.0 <= distance <= 6.0
    
    def calculate_dna_specific_terms(self, protein_atom: Atom, dna_atom: Atom,
                                   dna_residue: Residue, surface_distance: float) -> Dict[str, float]:
        """Calculate DNA-specific interaction terms"""
        if not self.include_dna_specific:
            return {'major_groove': 0.0, 'minor_groove': 0.0, 
                   'base_stacking': 0.0, 'phosphate_contact': 0.0}
        
        terms = {}
        
        # Major groove interactions
        if self.is_in_major_groove(protein_atom, dna_residue):
            terms['major_groove'] = 1.5 * np.exp(-surface_distance**2 / 4.0)
        else:
            terms['major_groove'] = 0.0
        
        # Minor groove interactions
        if self.is_in_minor_groove(protein_atom, dna_residue):
            terms['minor_groove'] = 1.2 * np.exp(-surface_distance**2 / 4.0)
        else:
            terms['minor_groove'] = 0.0
        
        # Phosphate backbone interactions
        if dna_atom.get_name() in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'']:
            distance = np.linalg.norm(protein_atom.get_coord() - dna_atom.get_coord())
            if distance < 5.0:
                terms['phosphate_contact'] = np.exp(-distance**2 / 8.0)
            else:
                terms['phosphate_contact'] = 0.0
        else:
            terms['phosphate_contact'] = 0.0
        
        # Base stacking interactions
        if (protein_atom.element == 'C' and 
            dna_atom.get_name() in ['C2', 'C4', 'C5', 'C6', 'C8']):
            if 3.0 <= surface_distance <= 4.5:
                terms['base_stacking'] = 1.0
            else:
                terms['base_stacking'] = 0.0
        else:
            terms['base_stacking'] = 0.0
        
        return terms
    
    def score_complex(self, structure_file: str) -> Dict[str, float]:
        """Score protein-DNA complex"""
        print(f"Loading structure: {structure_file}")
        
        structure = self.load_structure(structure_file)
        if structure is None:
            return {}
        
        # Identify chains
        chains = self.identify_chains(structure)
        protein_chains = chains['protein']
        dna_chains = chains['dna']
        
        print(f"Identified chains - Protein: {protein_chains}, DNA: {dna_chains}")
        
        if not protein_chains or not dna_chains:
            print("Error: Could not identify both protein and DNA chains")
            return {}
        
        # Extract atoms
        protein_atoms = self.extract_atoms(structure, protein_chains)
        dna_atoms = self.extract_atoms(structure, dna_chains)
        
        print(f"Found {len(protein_atoms)} protein atoms and {len(dna_atoms)} DNA atoms")
        
        # Initialize scoring terms
        total_scores = {
            'gauss1': 0.0, 'gauss2': 0.0, 'repulsion': 0.0,
            'hydrophobic': 0.0, 'hydrogen_bonding': 0.0, 'electrostatic': 0.0,
            'major_groove': 0.0, 'minor_groove': 0.0,
            'base_stacking': 0.0, 'phosphate_contact': 0.0
        }
        
        interaction_count = 0
        for i, (protein_atom, protein_res, protein_chain) in tqdm(enumerate(protein_atoms), total=len(protein_atoms), ncols=80, desc="Calculating interactions"):
            for dna_atom, dna_res, dna_chain in dna_atoms:
                center_distance = np.linalg.norm(
                    protein_atom.get_coord() - dna_atom.get_coord())
                
                if center_distance > self.cutoff_distance:
                    continue
                
                interaction_count += 1
                
                # Calculate surface distance
                surface_distance = self.calculate_surface_distance(protein_atom, dna_atom)
                
                # Standard Vina terms
                vina_terms = self.calculate_vina_terms(surface_distance)
                for term, value in vina_terms.items():
                    total_scores[term] += value
                
                # Additional terms
                total_scores['hydrophobic'] += self.calculate_hydrophobic_term(
                    surface_distance, protein_atom, dna_atom, protein_res, dna_res)
                
                total_scores['hydrogen_bonding'] += self.calculate_hydrogen_bonding(
                    surface_distance, protein_atom, dna_atom)
                
                total_scores['electrostatic'] += self.calculate_electrostatic_term(
                    protein_atom, dna_atom, protein_res, dna_res, center_distance)
                
                # DNA-specific terms
                dna_residue = None
                for model in structure:
                    for chain in model:
                        if chain.id == dna_chain:
                            for residue in chain:
                                for atom in residue:
                                    if atom == dna_atom:
                                        dna_residue = residue
                                        break
                
                if dna_residue:
                    dna_terms = self.calculate_dna_specific_terms(
                        protein_atom, dna_atom, dna_residue, surface_distance)
                    for term, value in dna_terms.items():
                        total_scores[term] += value
        
        print(f"Calculated {interaction_count} interactions")
        
        # Calculate weighted scores
        weighted_scores = {}
        total_weighted = 0.0
        
        for term, raw_score in total_scores.items():
            if term in self.weights:
                weighted_score = self.weights[term] * raw_score
                weighted_scores[f'{term}_weighted'] = weighted_score
                total_weighted += weighted_score
        
        # Estimate rotatable bonds
        n_rot = len(protein_atoms) // 10
        
        # Apply rotatable bond penalty
        final_score = total_weighted / (1 + self.weights['n_rot'] * n_rot)
        
        # Compile results
        results = {
            'final_score': final_score,
            'total_weighted_score': total_weighted,
            'n_rotatable_bonds': n_rot,
            'interaction_count': interaction_count,
            'protein_atoms': len(protein_atoms),
            'dna_atoms': len(dna_atoms),
            **total_scores,
            **weighted_scores
        }
        
        return results
    
    def print_results(self, results: Dict[str, float]):
        """Print detailed scoring results"""
        print("\n" + "="*60)
        print("VINA-LIKE SCORING RESULTS")
        print("="*60)
        
        print(f"\nFinal Binding Score: {results['final_score']:.3f} kcal/mol")
        print(f"Total Weighted Score: {results['total_weighted_score']:.3f}")
        print(f"Rotatable Bonds: {results['n_rotatable_bonds']}")
        print(f"Interactions: {results['interaction_count']}")
        
        print(f"\nStructure Information:")
        print(f"  Protein atoms: {results['protein_atoms']}")
        print(f"  DNA atoms: {results['dna_atoms']}")
        
        print(f"\nEnergy Components (Raw Scores):")
        energy_terms = ['gauss1', 'gauss2', 'repulsion', 'hydrophobic', 
                       'hydrogen_bonding', 'electrostatic']
        for term in energy_terms:
            if term in results:
                print(f"  {term:15s}: {results[term]:8.3f}")
        
        if self.include_dna_specific:
            print(f"\nDNA-Specific Terms:")
            dna_terms = ['major_groove', 'minor_groove', 'base_stacking', 'phosphate_contact']
            for term in dna_terms:
                if term in results:
                    print(f"  {term:15s}: {results[term]:8.3f}")
        
        print(f"\nWeighted Contributions:")
        for term in energy_terms + (['major_groove', 'minor_groove', 'base_stacking', 'phosphate_contact'] if self.include_dna_specific else []):
            weighted_key = f"{term}_weighted"
            if weighted_key in results:
                print(f"  {term:15s}: {results[weighted_key]:8.3f}")


def main():
    """Main function for command-line usage"""
    if len(sys.argv) != 2:
        print("Usage: python vina_dna_scorer_standalone.py <structure_file.pdb>")
        print("\nStandalone Vina-like Scoring for Protein-DNA Complexes")
        print("Features:")
        print("  - No MDAnalysis dependency")
        print("  - Only requires BioPython and NumPy")
        print("  - AMBER-based partial charges")
        print("  - DNA-specific interaction terms")
        print("  - Detailed energy breakdown")
        sys.exit(1)
    
    structure_file = sys.argv[1]
    
    if not os.path.exists(structure_file):
        print(f"Error: Structure file not found: {structure_file}")
        sys.exit(1)
    
    # Initialize scorer
    scorer = StandaloneVinaScorer(
        cutoff_distance=8.0,
        include_electrostatics=True,
        dielectric_constant=4.0,
        include_dna_specific=True
    )
    
    # Score the complex
    results = scorer.score_complex(structure_file)
    
    if results:
        scorer.print_results(results)
    else:
        print("Error: Could not score the complex")
        sys.exit(1)


if __name__ == "__main__":
    main()

