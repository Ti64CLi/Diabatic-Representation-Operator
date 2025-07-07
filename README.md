# Operator Representation Module

## English

This module provides a symbolic and symmetry-adapted construction of operator matrices between quantum states transforming under irreducible representations (irreps) of the \$C\_{nv}\$ point groups.

### Key Functionalities

#### 1. `operator`

```python
operator(n, opsymmetry, s1, s2, nvarsym, max_order) -> np.ndarray[tuple[Operator, Operator]]
```

**Description**: Constructs the matrix representation of an operator with symmetry `opsymmetry`, acting between two quantum states of symmetries `s1` and `s2`, using all symmetry-allowed invariants and monomials up to a given `max_order`.

- `n`: integer (\$C\_{nv}\$ group order)
- `opsymmetry`: `Symmetry` object (A1, A2, B1, B2, or E)
- `s1`, `s2`: `Symmetry` objects representing the state symmetries
- `nvarsym`: list[int], number of variables of each symmetry (indexed as [A1, A2, B1, B2, E\_1, E\_2, ...])
- `max_order`: maximum polynomial order

**Returns**: A 2x2 or 1x1 array of tuples `(A_x, A_y)`, where each `A_x`, `A_y` is an `Operator` object representing the X and Y components of the operator matrix.

---

#### 2. `A_x`

```python
A_x(n, opsymmetry, s1, s2, max_order) -> Operator
```

**Description**: Generates the \$X\$-component of the symbolic operator matrix using selection rules dictated by the \$C\_{nv}\$ symmetry.

Returns an `Operator` object representing the symmetry-allowed expansion in \$X\$ for the given transition.

---

#### 3. `A_y`

```python
A_y(n, opsymmetry, s1, s2, max_order) -> Operator
```

**Description**: Generates the \$Y\$-component of the symbolic operator matrix using the \$C\_{nv}\$ symmetry.

Returns an `Operator` object representing the symmetry-allowed expansion in \$Y\$.

---

#### 4. `Operator.compile`

```python
Operator.compile() -> str
```

**Description**: Compiles an `Operator` into a CSV-like text format containing:

- Operator and state symmetries
- Variable counts per irrep
- Definitions of pseudovariables (combinations of fundamental ones)
- Monomial expansions for X and Y components

This textual format can be used for interfacing with external computational tools or storing symbolic representations.

---

## Français

Ce module permet la construction symbolique d'opérateurs agissant entre états quantiques, en tenant compte des symétries ponctuelles \$C\_{nv}\$.

### Fonctionnalités principales

#### 1. `operator`

```python
operator(n, opsymmetry, s1, s2, nvarsym, max_order) -> np.ndarray[tuple[Operator, Operator]]
```

**Description** : Construit la matrice de l'opérateur de symétrie `opsymmetry` entre deux états de symétries `s1` et `s2`, à l'aide de tous les invariants et monômes compatibles jusqu'à un ordre maximal `max_order`.

- `n` : entier correspondant au groupe \$C\_{nv}\$
- `opsymmetry` : objet `Symmetry` (A1, A2, B1, B2 ou E)
- `s1`, `s2` : objets `Symmetry` décrivant la symétrie des états
- `nvarsym` : liste d'entiers donnant le nombre de variables par irrep ([A1, A2, B1, B2, E\_1, E\_2, ...])
- `max_order` : ordre maximal du développement

**Retourne** : Un tableau 2x2 (ou 1x1 si les états sont identiques) de couples `(A_x, A_y)`, où `A_x` et `A_y` sont des objets `Operator` contenant les matrices d'opérateurs en \$X\$ et \$Y\$.

---

#### 2. `A_x`

```python
A_x(n, opsymmetry, s1, s2, max_order) -> Operator
```

**Description** : Génère la composante \$X\$ de la matrice de l'opérateur à partir des règles de sélection imposées par le groupe \$C\_{nv}\$.

Retourne un objet `Operator` contenant l'expansion symbolique selon \$X\$.

---

#### 3. `A_y`

```python
A_y(n, opsymmetry, s1, s2, max_order) -> Operator
```

**Description** : Génère la composante \$Y\$ de la matrice de l'opérateur à partir des contraintes de symétrie.

Retourne un objet `Operator` contenant l'expansion symbolique selon \$Y\$.

---

#### 4. `Operator.compile`

```python
Operator.compile() -> str
```

**Description** : Compile un objet `Operator` en un format texte de type CSV contenant :

- Les symétries de l'opérateur et des états
- Le nombre de variables par représentation irréductible
- La définition des pseudo-variables
- Les composantes X et Y en monômes

Ce format est adapté pour l'exportation, la sauvegarde ou une lecture par des outils externes.

---

### Authors / Auteurs

- Elie DUMONT
- David WILLIAMS
- Alexandra VIEL

---

### License / Licence

MIT

---

