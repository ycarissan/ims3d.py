# Sur votre machine
## Créer un environnement pour le programme ims3d

Dans un terminal, récupérer le programme ims3d :
```
git clone https://github.com/ycarissan/ims3d.py.git
echo "export IMS3D_PATH=$PWD/ims3d.py" >> .bashrc
export IMS3D_PATH=$PWD/ims3d.py
```
Récupérer le programme de gestion des environnements virtuels conda :
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate

echo "source ~/miniconda3/bin/activate" >> ~/.bashrc
```
Générer l'environnement pour ims3d :
```
cd ${IMS3D_PATH}
conda env create -f conda/ims3d_conda_env.yml
```
Activer l'environnement (il faudra le faire dans chaque nouveau terminal dans lequel vous voudrez utiliser le programme ims3d) : 
```
conda activate ims3d_env
```
## Monter les disques du mésocentre
On crée un répertoire dans lequel les fichiers du mesocentre vont apparaitre:
```
mkdir $HOME/_mesocentre
```
On monte les fichiers
```
sshfs utilisateur@mesocentre.univ-amu.fr: _mesocentre
```
## Générer les fichiers d'entrée pour ims3d
Admettons que vous souhaitiez faire vos calculs pour la molécule benzene dans le répertoire `$HOME/calculs/benzene` au mésocentre.
Ce répertoire apparaît sur `$HOME/_mesocentre/calculs/benzene` une fois les disques du mésocentre montés.
Allons-y
```
cd $HOME/calculs/benzene
ls
   benzene.xyz
```
Lançons la génération des fichiers d'entrée
```
conda activate ims3d_env # Si ce n'est pas déjà fait dans le terminal
python ${IMS3D_PATH}/ims3d.py -r 1 -b 400 -f orca benzene.xyz
```
La grille est créée et les fichiers d'entrée sont générés.
Rendons-nous sur le cluster du mésocentre :
```
ssh utilisateur@mesocentre.univ-amu.fr
```

# Sur le cluster
On active la possibilité d'utiliser les programmes du labo (à ne faire qu'une seule fois):
```
echo "export PATH=${PATH}:/home/ycarissan/CTOM_prog/bin" >> ~/.bashrc
source ~/.bashrc
```
On se rend dans le bon répertoire
```
cd calculs/benzene
```
On récupère un fichier pour lancer orca sur tous les fichier d'entrée normalement nommés `input_batchXXXXX.inp`, XXXXX étant un nombre entier codé sur 5 caractères.
```
get_orcafile_ims
```
Pour lancer les calculs 
```
sbatch submit_ims.job "1-$(ls *.txt|wc -l)"
```
Maintenant il faut attendre que les calculs soient terminés. Pour cela on vérifie que la queue de calcul est vide des jobs qu'on vient de lancer:
```
squeue --me
```
S'il reste des jobs, on attend. Un jour il n'y en a plus.
# De retour sur votre machine
Une fois tous les calculs terminés, on doit récupérer les données :
```
cd _mesocentre/calculs/benzene
python ${IMS3D_PATH}/ims3d_harv.py -f orca input_batch00000.out
```
Un fichier non vide nommé `ims.dat` doit exister. Alors on génère le fichier `test.vtk`
```
python ${IMS3D_PATH}/ims3d_view.py
```
Le programme crashe mais ce n'est pas grave, le fichier est créé.
On peut le visualiser avec `paraview` : https://www.paraview.org/.
