# Code pour analyser and evaluer les performances de DeepSC sur des données au data-tier RECO
Les données RECO sont celles en sortie de l'algorithme de reconstruction

## Reco dumper
Ceci est le procédé général pour produire dans un format pratique les données venant de l'algorithme de reconstruction. Il permet d'obtenir les données pour **Moustache** ou pour **DeepSC** avec les stratégies A, B, C ou D pour des données produites par ParticleGun pour électrons/photons. Ce dumper transforme des données stockées dans des arbres (.root) en données stockées sous tableaux panda. Les notebooks analysent les données stockées dans ces tableaux. Le script réalisant le dumper est `reco_dumper.py`, les autres permettent d'automatiser le procédé.

Les différents fichiers/dossiers ont le rôle suivant :
- notebooks : scripts d'analyse des données
- simScore : fichiers produits par [truth definition](https://github.com/valsdav/DeepSuperCluster/tree/master/NtuplesProduction/input_dataset_truth)
