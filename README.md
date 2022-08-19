# Code pour analyser and evaluer les performances de DeepSC sur des données au data-tier RECO
Les données RECO sont celles en sortie de l'algorithme de reconstruction

## Reco dumper
Ceci est le procédé général pour produire dans un format pratique les données venant de l'algorithme de reconstruction. Il permet d'obtenir les données pour **Moustache** ou pour **DeepSC** avec les stratégies A, B, C ou D pour des données produites par ParticleGun pour électrons/photons. Ce dumper transforme des données stockées dans des arbres (.root) en données stockées sous tableaux panda. Les notebooks analysent les données stockées dans ces tableaux. Le script réalisant le dumper est `reco_dumper.py`, les autres permettent d'automatiser le procédé.

Les différents fichiers/dossiers ont le rôle suivant :
- notebooks : scripts d'analyse des données
- simScore : fichiers produits par [truth definition](https://github.com/valsdav/DeepSuperCluster/tree/master/NtuplesProduction/input_dataset_truth), ils contiennent les "sim fraction thresholds" (*Expliquer ce que c'est*)
- `reco_dumper.py` est le script du dumper
- `run_reco_dumper.py` est le script qui permet de lancer le dumper avec les bonnes options (voir exemple dans la suite)
- `condor_reco_dumper.py` permet de lancer le dumper en utilisant HTCondor (permet de paralléliser les tâches efficacement) pour les calculs longs
- `Moustache.C` et `calo_association.py` font quoi ????
- `join_datasets.py` est un script a lancer après `condor_reco_dumper.py` pour agréger l'ensemble des données en un seul fichier panda
