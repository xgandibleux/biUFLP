# Branche didactique / Luxor

Ce contenu est destiné à une branche git dédiée (ex. `didactic-luxor`), séparée
du code public principal (branche `main`).

## Ce qui diffère de la branche `main`

Deux fichiers supplémentaires dans `src/`, absents de `main` :

- **`src/main_didactic.jl`** — point d'entrée minimaliste, restreint à
  l'instance didactique (`data/dataDidactic/didactic1.txt` ou `didactic2.txt`).
  Pas de mode expérimentation, pas de sélection de jeu de données, pas de
  bascule `vOptSolver`/`backendGR` — juste : charger l'instance didactique,
  résoudre, tracer avec Luxor. Pensé pour une démonstration/présentation.

- **`src/graphicluxor.jl`** — le backend graphique Luxor, rendu autonome
  (il définit sa propre constante `CanvasSize`, qui vivait auparavant dans
  `datastru.jl` sur la branche principale). Ce backend n'a de sens que sur
  l'instance didactique : il utilise des index de boîtes codés en dur
  (`ib1=3`, `ib=3`) et une échelle de canevas fixe adaptée aux coûts de
  l'exemple pédagogique (< 100), pas aux jeux de données F/H.

Tous les autres fichiers (`datastru.jl`, `parser.jl`, `biUFLP.jl`, `reduce.jl`,
`mopRoutines.jl`, `pavingBnB.jl`, `reducePaving.jl`,
`labelingBourrinOrdrevct.jl`, `main.jl`, `graphicpyplot.jl`) sont identiques
à la branche `main` — le socle algorithmique n'est pas dupliqué ni divergent.

## Comment lancer

```
cd src
julia main_didactic.jl
```

Produit un fichier `figUFLPbox.png` (via Luxor) illustrant le pavage de
l'instance didactique.
