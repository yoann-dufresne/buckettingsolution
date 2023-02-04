# Bucketing de superkmers

Votre objectif est de créer un logiciel capable, depuis un fichier fasta, de calculer tous les k-mers d'un ensemble de séquences ADN, puis des les distribuer en "buckets" sous la forme de superkmers.
Tous les termes sont expliqués un plus loin dans ce README.

ATTENTION : Veillez à bien tout lire de ce fichier avant de vous lancer dans la programmation.


# Projet global

Nous vous demandons de produire un logiciel qui prend en ligne de commande:
- 1 fichier au format [FASTA](https://fr.wikipedia.org/wiki/FASTA_(format_de_fichier))
- Un chemin vers un répertoire où seront stockées les sorties
- Un entier k (inférieur ou égal à 31)
- Un entier m inférieur à k représentant la taille des minimiseurs

Le fichier d'entrée sera forcément au format fasta et les fichiers de sortie au format texte avec un superkmer par ligne.

Exemple de ligne de commande:

```bash
  #         fichier          répertoire out  k  m    
  ./monprog data/ecoli.fasta data/ecoli_out/ 31 13
```

## Attentes

- Un repo Github contenant votre projet dont le lien nous sera envoyé avant Lundi 13 Février.
- Un readme dans le Github expliquant comment compiler et exécuter votre code.
- Un code qui compile
- Un code qui s'exécute avec la ligne de commande que vous avez donné.
- Des fichiers de code commentés
- Un code propre et reprenable par quelqu'un d'autre.
- Langage de programmation: C/C++/Rust fortement recommandés.


## Conseil d'ordre d'étapes du projet

Nous vous conseillons de procéder par résolution des étapes suivantes dans l'ordre :  
- Télécharger un génome (voir ci-dessous, commencez par le plus petit)  
- Faites-vous une fonction qui vous permet de lire le fichier fasta de manière à ce que chaque appel à la fonction vous renvoie uniquement le caractère suivant de la séquence ADN en cours de lecture  
- Faites une fonction qui encode une chaine de caractères de taille k en une valeur entière (voir la partie "Chaines de caractère vs entiers").
- Créez un fonction qui à partir d'un kmer précédent (sous format d'entier) et une lecture du caractère suivant dans le fichier vous donne l'entier correspondant au kmer suivant.  
- Listez tous les m-mers dans l'ordre du fichier, pour en déduire des superkmers
- Codez une fonction qui ajoute un superkmer à la fin d'un fichier texte
- Bravo, vous avez fini le code


# kmers, minimiseurs et superkmers

## k-mer

Un k-mer est une chaine de caractère de taille k.
Dans notre cas, nous nous intéresserons à tous les kmers d'une ou plusieurs séquences.

Ensemble des 5-mers pour une séquence:
```
  seq  : TTAGGACAAA
  kmers: TTAGG
          TAGGA
           AGGAC
            GGACA
             GACAA
              ACAAA
```


## Minimiseur

Un minimiseur est le m-mer le plus petit (selon une relation d'ordre prédéfinie) au sein d'un kmer.
Pour calculer un minimiseur, il est donc nécessaire de calculer tous les m-mers d'un kmer.
Parmi tout ces m-mer il faut ensuite ne conserver que le plus petit selon une relation d'ordre.
Comme relation d'ordre pour ce projet, on peut choisir l'ordre alphabétique.

En reprenant le premier kmer de l'exemple précédent avec m=3 on a donc:
```
kmer : TTAGG
mmers: TTA
        TAG
         AGG <-- minimiseur
```

En temps normal, il faut se préoccuper de l'ensemble des m-mers ainsi que de leur reverse-complement pour déterminer le minimiseur d'un kmer.
Ici, nous vous demandons seulement de vérifier parmi les mmers forward.


## Un peu de contexte

Les minimiseurs sont principalement utilisés comme signatures de kmers.
On peut séparer les kmers en "buckets" où chaque bucket contient tous les kmers avec le même minimiseur.
Ensuite lorsqu'on cherche un kmer, on peut donc aller directement dans le bon bucket le chercher.
Cela permet d'éviter de stocker tous les kmers d'un gros génome dans une gigantesque structure de données.
À la place, on préfère un grand nombre de petites structures, 1 par minimiseur.
Pour ce projet, nous ne programmerons que la séparation des kmers en buckets et pas les petites structures de données.

## Superkmer

Remarque: 2 kmers successifs ont de fortes chances de partager les même minimiseur puisqu'ils partagent tous leurs m-mers sauf les 2 aux extrémités.

On appelle superkmer la compaction des kmers successifs partageant le même minimiseur (même == même valeur et même position).

Reprenons le même exemple qu'un peu plus haut et déterminons le minimiseur pour chaque kmer:
```
  seq  : TTAGGACAAA  | Minimiseur
  kmers: TTAGG       |   AGG  
          TAGGA      |   AGG
           AGGAC     |   AGG
            GGACA    |   ACA
             GACAA   |   ACA
              ACAAA  |   AAA
```

En compactant les kmers successifs avec le même minimiseur, on obtient la liste des 3 superkmers suivants:
```
TTAGGAC
GGACAA
ACAAA
```

## Chaines de caractère vs entiers

En informatique, comparer 2 chaines de caractères coûte autant d'opérations que le nombre de caractères à comparer.
Dans notre cas, on va devoir vérifier si un m-mer est plus petit qu'un autre.
Il va donc falloir comparer des chaines de caractères (ce qui risque d'être lent).

Une astuce pour remplacer cette comparaison par une comparaison d'entiers.
Pour cela vous devrez utiliser un encodage 2 bits (base 4) pour transformer une chaine.
Par exemple avec l'encodage A:0 C:1 T:2 G:3 le kmer CTT vaut $1 * 4^2 + 2 * 4^1 + 2 * 4^0 = 26$ (01 10 10 en binaire).

L'avantage d'utiliser ce type d'encodage est qu'ils sont "glissants".
A partir d'un kmer encodé, vous pouvez "retirer" les 2 bits de poids fort (le nucléotide le plus à gauche), décaler l'entier de 2 bits vers la gauche puis ajouter les 2 bits du nouveau nucléotide à droite pour obtenir le kmer suivant (le mmer dans ce projet).


# Datasets : Séquences ADN

## Téléchargement

Pour télécharger un génome il faut:  
1 - Cliquer sur le lien correspondant au génome (voir ci-après)  
2 - Cliquer en haut à droite sur "send to"  
3 - Sélectionner "File"  
4 - Modifier le format pour mettre "FASTA"  
5 - Cliquer sur "Create File"  

- Virus (Sars-Cov-2): https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2  
- Bactérie (e.coli): https://www.ncbi.nlm.nih.gov/nuccore/NZ_LN874954.1  
- Mammifère (Opossums gris à queue courte): https://www.ncbi.nlm.nih.gov/nuccore/NC_008801.1

## Lecture du fichier .fasta

Les fichiers fasta ci-dessus sont tous formatés de la même manière.
La première ligne est une entête dont le premier caractère est un Chevron fermant.
Les lignes suivantes représentent la séquence ADN.
Les retours à la ligne sont uniquement pour faciliter la lecture d'un Humain et ne représentent rien.
Vous devez donc les ignorer lors de la lecture du fichier.  

Les séquences sont composées principalement des 4 lettres A, C, G, T.
Parfois, les meilleurs machines actuelles n'ont pas réussit à parfaitement lire certains nucléotides.
Ces nucléotides sont remplacés par des N (dans la bactérie et l'humain).
Pour ce projets, ignorez les N.
Ainsi, si la séquence est ACTTNNNNATNGCT considérez à la place que c'est ACTTATGCT.
TTA est donc un 3-mer valide ici.
Votre lecteur de fichier dois passer par dessus ces caractères sans les retourner.
ATTENTION : Lorsque nous testerons votre code, nous utiliserons des séquences avec des N. Votre programme ne dois pas planter.

