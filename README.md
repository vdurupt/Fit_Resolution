# Code pour analyser and evaluer les performances de DeepSC par rapport à Moustache (sur des données au data-tier RECO)
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

## Comment faire tourner le dumper
Pour faire tourner le dumper avec condor, il faut éxécuter `condor_reco_dumper.py` :
`python3 condor_reco_dumper.py` 

Les options sont :
```  -h, --help            pour afficher les options

  -i INPUTDIR, --inputdir INPUTDIR
                        répertoire des données de simulation
                        
  -nfg NFILE_GROUP, --nfile-group NFILE_GROUP
                        Nombre de fichiers par fichier numpy
                        
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        répertoire de sortie du script
                        
  -a ASSOC_STRATEGY, --assoc-strategy ASSOC_STRATEGY
                        Stratégie d'association (sim_fraction voir [truth definition](https://github.com/valsdav/DeepSuperCluster/tree/master/NtuplesProduction/input_dataset_truth))
                        
  --wp-file WP_FILE     Fichier avec les limites de sim_fraction
  
  -q QUEUE, --queue QUEUE
                        File d'attente de condor (donne la priorité du calcul, voir [flavors](https://batchdocs.web.cern.ch/local/submit.html))
                        
  -e EOS, --eos EOS     EOS instance user/cms
  
  -c, --compress        Compresse les fichier d'output (à utiliser si on veut concaténer les résultats)
  
  --redo                Redo all files
  
  -d, --debug           debug
  
  --loop-on-calo        Si cette option est vraie, le script boucle seulement sur les calo-seeds et pas sur l'ensemble des SC
  
  -s SC_COLLECTION, --sc-collection SC_COLLECTION
                        SuperCluster collection
                        
  -r RECO_COLLECTION, --reco-collection RECO_COLLECTION
                        Reco collection (none/electron/photon)
                        
  -cf CONDOR_FOLDER, --condor-folder CONDOR_FOLDER
                        Dossier condor où sont stockées outputs de condor (condor_ndjson par défaut)
```


Pour analyser les photons (en utilisant **DeepSC** avec la stratégie de collection A) en bouclant sur les objets reconstruits (pour étudier la reconstruction par exemple), j'ai utilisé la commande suivante :

``` python3 condor_reco_dumper.py -i /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourGammasGunPt1-100_pythia8_StdMixing_Flat55To75_14TeV_123X_mcRun3_2021_realistic_v11_UL18_pfRechitThres_Dumper_SCRegression_PhoRegression_DeepSC_AlgoA_125X_bugFix -o /eos/user/v/vdurupt/reco_comparison_corrected/photons/pho_UL18_123X_DeepSC_AlgoA/ -a sim_fraction --wp-file simScore_Minima_PhotonsOnly_updated.root -nfg 40 -q espresso --compress --reco-collection photon ```

Attention, il faut que le répertoire de sortie `-o` soit dans EOS.
Une fois cette commande effectuée, le répertoire `CONDOR_FOLDER` est créé avec les différents fichiers nécessaires à la demande de calcul pour condor. Il faut donc aller dans ce répertoire et soumettre le travail à condor par la commande : `condor_submit condor_job.txt`

Pour vérifier l'état d'avancement du travail, on peut utiliser `condor_q` et pour enlever un travail de la liste d'attente, la commande est `condor_rm JOBID`.
Une fois que ces calculs sont terminés, les outputs sont disponibles dans le répertoire spécifié dans l'option `-o` de `condor_reco_dumper.py`. Il reste une étape pour le tableau panda plus facilement utilisable : la concaténation des output.

Pour cela la commande est la suivante :
`python3 join_datasets.py -i /eos/user/v/vdurupt/reco_comparison_corrected/photons/pho_UL18_123X_DeepSC_AlgoA/ -o /eos/user/v/vdurupt/reco_comparison_corrected/supercluster_regression/photons/pho_UL18_123X_DeepSC_AlgoA_{type}.h5py`
Avec `-i` le répertoire où sont stockés les données précédentes.

## Analyse DeepSC vs. Moustache
Les données (sous le format `.h5py` de panda) sont ensuite utilisées pour comparer DeepSC à Moustache ou pour comparer les différentes stratégies de collection de DeepSC. Les notebooks utilisés pour réaliser ces analyses se trouvent dans le répertoire [notebooks](https://github.com/vdurupt/Fit_Resolution/notebooks).
