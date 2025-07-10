## Technical Details

### Scoring Function Mathematics

#### Core Vina Function

The total score is calculated as:

```
Score = Σ w_i * f_i(d_ij) / (1 + w_rot * N_rot)
```

Where:
- `w_i`: Weight for energy term i
- `f_i(d_ij)`: Energy function for atom pair at surface distance d_ij
- `w_rot`: Rotatable bond weight (0.0585)
- `N_rot`: Number of rotatable bonds

#### Surface Distance Calculation

Surface distance between atoms i and j:

```
d_ij = r_ij - R_i - R_j
```

Where:
- `r_ij`: Center-to-center distance
- `R_i`, `R_j`: Van der Waals radii

#### Energy Function Definitions

**Gaussian Terms:**
```
gauss1(d) = exp(-(d/0.5)²)
gauss2(d) = exp(-((d-3)/2)²)
```

**Repulsion Term:**
```
repulsion(d) = d² if d < 0, else 0
```

**Hydrophobic Term:**
```
hydrophobic(d) = 1 if d < 0.5 Å
                = 0 if d > 1.5 Å
                = linear interpolation otherwise
```

**Hydrogen Bonding:**
```
hbond(d) = 1 if d < -0.7 Å
         = 0 if d > 0 Å
         = linear interpolation otherwise
```

#### Electrostatic Calculations

Coulombic interactions with distance-dependent dielectric:

```
E_elec = (332.0636 * q_i * q_j) / (ε_eff * r_ij)
```

Where:
- `q_i`, `q_j`: Partial charges (from AMBER force field)
- `ε_eff`: Effective dielectric constant
- `r_ij`: Interatomic distance

The effective dielectric constant varies with distance:
```
ε_eff = max(1.0, ε_0 * r_ij / 5.0)
```

This prevents singularities at short distances while maintaining realistic electrostatic interactions.

#### DNA-Specific Terms

**Major Groove Interactions:**
```
E_major = 1.5 * exp(-d²/4) if atom in major groove
```

**Minor Groove Interactions:**
```
E_minor = 1.2 * exp(-d²/4) if atom in minor groove
```

**Phosphate Backbone Contacts:**
```
E_phosphate = exp(-r²/8) if r < 5.0 Å and atom is phosphate
```

**Base Stacking:**
```
E_stacking = 1.0 if 3.0 Å ≤ d ≤ 4.5 Å and atoms are aromatic
```

### Algorithm Implementation

#### Chain Identification Algorithm

1. **Residue Analysis**: Count nucleotide vs. amino acid residues per chain
2. **Classification**: Chains with >50% nucleotides classified as DNA
3. **Validation**: Check for minimum chain lengths and valid residue types

#### Interaction Calculation

1. **Neighbor Search**: Use spatial partitioning for efficient distance calculations
2. **Cutoff Application**: Skip atom pairs beyond cutoff distance
3. **Energy Evaluation**: Calculate all energy terms for valid pairs
4. **Summation**: Aggregate weighted energy contributions

#### Performance Optimizations

- **Vectorized Operations**: Use NumPy for efficient array operations
- **Distance Cutoffs**: Early termination for distant atom pairs
- **Sparse Calculations**: Only compute non-zero energy terms
- **Memory Management**: Efficient data structures for large complexes