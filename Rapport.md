# Rapport TP5

Dans les exercices suivants, le but sera d'implémenter des algorithmes pour le calcul scientifique, notamment en algèbre linéaire, afin de comprendre et
comparer les performances et les limites de ces derniers.
Nous testerons et validerons les algorithmes écris en calculant leur erreur avant.

Pour rappel, les algorithmes implémentés ont été vus en cours et sont :
- CSR/CSC
- LU tridiagonal
- Jacobi & Gauss-Seidel
- Richardson

Tous les graphiques ont été réalisés à partir des valeurs calculées par nos implémentations et écrites dans des fichiers de sortie au format `.dat` (accessibles 
dans le dossier `graphs`).
Ces données ont ensuite été traitées grâce à l'outil [Gnuplot](http://www.gnuplot.info/) (dont les scripts sont disponibles dans le dossier `plot`).

## Méthode directe et stockage bande
### Exercice 3
1. Pour déclarer et allouer une matrice en C pour utiliser BLAS et LAPACK, on utilise la fonction "malloc" pour allouer de l'espace mémoire pour la matrice, puis on utilise un pointeur pour y accéder. 
```
AB  = (double  *) malloc(sizeof(double)*lab*la);
```
2. LAPACK_COL_MAJOR est une constante utilisée pour spécifier l'ordre de stockage des données dans une matrice en utilisant la bibliothèque de calcul scientifique LAPACK. La valeur de cette constante est 1 et indique que les données sont stockées en column-major, c'est-à-dire que les éléments d'une colonne sont stockés les uns à côté des autres en mémoire. Cela est en contraste avec l'ordre de stockage row-major, où les éléments d'une ligne sont stockés les uns à côté des autres.<br><br>
3. La dimension principale (leading dimension) généralement notée ld est utilisée pour spécifier l'espacement entre les éléments consécutifs d'une matrice en mémoire. Elle est souvent utilisée dans les bibliothèques de calcul scientifique LAPACK pour permettre aux matrices de disposer de suffisamment d'espace pour stocker des sous-matrices ou des extensions.<br><br>
4. La fonction dgbmv (Double-precision General Banded Matrix-Vector product) est une fonction qui permet de faire le produit matrice-vecteur entre une matrice bandée générale GB et un vecteur.
Elle calcule $(y = alpha * A * x + beta * y)$,  où $A$ est une matrice bandée générale avec $kl$ diagonales inférieures et $ku$ diagonales supérieures, $x$ est le vecteur entrée, $y$ est le vecteur de sortie, $alpha$ et $beta$ sont des scalaires. La fonction dgbmv de LAPACK implémente la méthode du produit matrice-vecteur traditionnel pour calculer le produit matrice-vecteur entre une matrice bandée générale et un vecteur.<br><br>
5. La fonction dgbtrf (Double-precision General Banded Matrix LU factorization) est une fonction qui permet de calculer la factorisation LU d'une matrice bandée générale. Elle permet de décomposer une matrice bandée générale $A$ en deux matrices triangulaires $L$ et $U$ de telle sorte que $A = L * U$. Cette fonction implémente une méthode de pivot de Gauss pour obtenir la factorisation LU. Elle utilise les propriétés spécifiques de la matrice bandée générale pour optimiser les calculs en n'utilisant que les éléments non-nuls de la matrice. La fonction dgbtrf remplit la matrice $A$ avec les matrices $L$ et $U$. <br><br>
6. La fonction dgbtrs (Double-precision General Banded Matrix triangular Solver) est une fonction de la bibliothèque LAPACK qui permet de résoudre un système linéaire $Ax = b$ avec une matrice bandée générale A. Elle utilise la factorisation LU de la matrice A obtenue par la fonction dgbtrf pour résoudre le système linéaire. Cette fonction implémente les matrices L et U obtenues par la factorisation LU de la matrice A pour résoudre le système linéaire en utilisant des substitutions successives. Elle utilise également les permutations de lignes effectuées par la fonction dgbtrf pour rendre la matrice A échelonnée.<br><br>
7. La fonction dgbsv (Double-precision General Banded Matrix Solver) est une fonction de la bibliothèque LAPACK qui permet de résoudre un système linéaire $Ax = B$ avec une matrice bandée générale A. Elle utilise la factorisation LU pour résoudre le système linéaire de manière efficace. Cette fonction implémente une méthode de pivot de Gauss pour obtenir la factorisation LU de la matrice A. Elle utilise également les permutations de lignes effectuées pour rendre la matrice A échelonnée pour résoudre le système linéaire en utilisant des substitutions successives. La fonction remplit la matrice B avec la solution x du système.
8. Pour calculer la norme du résidu relatif d'un système linéaire $Ax = B$, on utilise les fonctions BLAS pour calculer la norme de la différence entre le vecteur B et le produit de la matrice A par la
solution x. La norme du résidu relatif est définie comme :
$$\frac{||b - Ax||}{||b||}$$
Pour calculer cette valeur, on utilise les fonctions BLAS suivantes : 
- dgbmv (ou dgemv si A est dense) pour calculer le produit matrice-vecteur Ax.
- daxpy pour calculer la différence B – Ax.
- dnrm2 pour calculer la norme de la différence B – Ax.
- dnrm2 pour calculer la norme de B.

### Exercice 4
3. Pour valider cette de méthode de résolution par la fonction dgbmv on calcul l’erreur relative :
$$relres = \frac{||RHS - EX\_SOL ||}{EX\_SOL }$$
### Exercice 6
2. Pour valider cette de méthode de factorisation on calcul l’erreur relative suivante :
$$\frac{||A- L*U||}{||A||}$$
## Méthode de résolution itérative
l’implémentation des algorithmes itératives ainsi que leur données et graphes correspondants ce trouve dans le repo github.

