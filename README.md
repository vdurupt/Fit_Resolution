# Code pour analyser and evaluer les performances de DeepSC sur des données au data-tier RECO
Les données RECO sont celles en sortie de l'algorithme de reconstruction.

## Reco dumper
Ceci est le procédé général pour produire dans un format pratique les données venant de l'algorithme de reconstruction. Il permet d'obtenir les données pour **Moustache** ou pour **DeepSC** avec les stratégies A, B, C ou D pour des données produites par ParticleGun pour électrons/photons. Ce dumper transforme des données stockées dans des arbres (.root) en données stockées sous tableaux panda. Les notebooks analysent les données stockées dans ces tableaux. Le script réalisant le dumper est `reco_dumper.py`, les autres permettent d'automatiser le procédé.

Les différents fichiers/dossiers ont le rôle suivant :
- notebooks : scripts d'analyse des données
- simScore : fichiers produits par [truth definition](https://github.com/valsdav/DeepSuperCluster/tree/master/NtuplesProduction/input_dataset_truth), ils contiennent les "sim fraction thresholds" (*Expliquer ce que c'est*)
- `reco_dumper.py` est le script du dumper
- `run_reco_dumper.py` est le script qui permet de lancer le dumper avec les bonnes options (voir exemple dans la suite)
- `condor_reco_dumper.py` permet de lancer le dumper en utilisant HTCondor (permet de paralléliser les tâches efficacement) pour les calculs longs
- `Moustache.C` et `calo_association.py` sont des scripts appelés par les scripts de plus haut niveau
- `join_datasets.py` est un script a lancer après `condor_reco_dumper.py` pour agréger l'ensemble des données en un seul fichier panda

Le dumper peut fonctionner de deux manières différentes :
- **loop-on-calo** :  information is saved for each PFCluster seed matched to a CaloParticle.
- **loop-on-reco** : the dumper works on the final SC or reco (electron/photon) collection and saves information for each reco object

Pour chaque objet (photon/électron), les informations suivantes sont sauvegardées :
- informations pour le calomatching (matcher les traces du tracker avec les clusters)
- informations pour le genmatching by DeltaR (matcher les caloparticules aux particules générées)
- position et énergie de chaque seed et objet
- nombre de clusters sélectionnés
- Gen energie (énergie de la particule générée) et sim energie (énergie reconstituée par la simulation utilisée tracker/ECAL) pour les pour les caloparticules/genparticules matchées
- information sur le pileup (PU)
- Number of event and run to be able to perform matching in the case the reconstruction is done two times with different algos.

##Comment faire tourner le dumper
Pour faire tourner le dumper avec condor, il faut éxécuter `condor_reco_dumper.py` :
`python3 condor_reco_dumper.py` avec les options suivantes (dans l'ordre) :
`usage: condor_reco_dumper.py [-h] -i INPUTDIR -nfg NFILE_GROUP -o OUTPUTDIR -a ASSOC_STRATEGY [--wp-file WP_FILE] -q QUEUE [-e EOS] [-c] [--redo] [-d] [--loop-on-calo] [-s SC_COLLECTION] [-r RECO_COLLECTION] [-cf CONDOR_FOLDER]`

Les options sont :
`  -h, --help            show this help message and exit
  -i INPUTDIR, --inputdir INPUTDIR
                        Inputdir
  -nfg NFILE_GROUP, --nfile-group NFILE_GROUP
                        How many files per numpy file
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Outputdir
  -a ASSOC_STRATEGY, --assoc-strategy ASSOC_STRATEGY
                        Association strategy
  --wp-file WP_FILE     File with sim fraction thresholds
  -q QUEUE, --queue QUEUE
                        Condor queue
  -e EOS, --eos EOS     EOS instance user/cms
  -c, --compress        Compress output
  --redo                Redo all files
  -d, --debug           debug
  --loop-on-calo        If true, loop only on calo-seeds, not on all the SC
  -s SC_COLLECTION, --sc-collection SC_COLLECTION
                        SuperCluster collection
  -r RECO_COLLECTION, --reco-collection RECO_COLLECTION
                        Reco collection (none/electron/photon)
  -cf CONDOR_FOLDER, --condor-folder CONDOR_FOLDER
                        Condor folder`
