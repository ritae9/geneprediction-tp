# geneprediction-tp
# TP Prédiction de gènes
  
Les gènes correspondent à une sous-séquence des transcripts traduites en protéines par le ribosome. Ils comportent un cadre de lecture consistant en des triplets consécutifs depuis un codon d’initiation ( 'AUG', 'UUG', 'CUG', 'AUU' ou 'GUG') et un codon stop (UAA', 'UAG', ou 'UGA'). Ces codons sont dans le même cadre de lecture ! 
On retrouve en amont du codon d’initiation un motif permettant l'initiation de la traduction via la fixation de la sous-unité 16S de l’ARN ribosomique : AGGAGGUAA appelée séquence de Shine-Dalagarno  [Shine et Dalgarno 1974]. Ce motif n'est pas nécessairement dans le même cadre de lecture que le codon d'initiation et peut être incomplet.

Dans le dossier geneprediction-tp/data/, on a :
listeria.fna : génome de Listeria monocytogenes
proteins_665_300274.csv : fiche NCBI des protéines de listeria
prodigal.csv: position des gènes prédites par prodigal
position.csv: position des gènes vérifiée expérimentalement

On a édité le programme Python3 nommé gpred.py dans le dossier gpred/.  Il prendra en argument:
- i fichier fasta avec des séquences nucléotidiques (.fna)
-g la longueur minimale des gènes
-s la distance maximale du codon initiateur où chercher la séquence de Shine-Dalgarno
-d le gap minimal entre deux gènes (Séquence de Shine-Dalgarno indépendant)
-p fichier de sortie des positions brutes des gènes
-o fichier fasta donnant les séquences fasta des gènes prédits

Notre programme utilisera impérativement ces arguments.

Lecture du fichier de séquences

read_fasta  prend un seul argument correspondant au fichier fasta (fasta_file) et retourne une séquence sous la forme d’un string (sans retour chariot). Par souci de simplicité, le fichier fasta ne contiendra qu’une seule molécule d’ADN. read_fasta s’assurera que les caractères de cette string sont bien en majuscule.

Recherche des codons d’initiation, stop et de la séquence de Shine-Dalgarno 

La recherche des motifs sera effectuée à l’aide de regex compilées générés par la fonction re.compile. Nous remplacerons l‘uracile par la thymine car nous effectuons cette recherche directement sur le génome ! 
regexp_object = re.compile(“myregex”)
Deux fonctions s’appliquant aux regex_object vont nous intéresser:
regexp_object.search(sequence, [start, [stop]]) : Cherche la première occurrence d’une regex entre une position start et stop de la séquence et retourne un objet: match_object
regexp_object.finditer(sequence, [start,[stop]]): retourne un itérateur de match_object non-chevauchant


find_start prend 4 arguments
start_regex: regex object permettant d’identifier un codon start
sequence: Séquence du génome
start : position de début de la recherche
stop : position de fin de la recherche
find_start retourne la position de la première occurrence d’un codon d'initiation dans la zone de recherche et None si rien n’est identifié. S

find_stop 
stop_regex: regex object permettant d’identifier un codon stop
sequence: Séquence du génome
start : position de début de la recherche
find_stop retourne le premier codon stop se trouvant dans le même cadre de lecture que le codon d'initiation et None si rien n’est identifié

has_shine_dalgarno
shine_regex: regex object permettant d’identifier une séquence de Shine-Dalagarno
sequence: Séquence du génome
start: position de début de la recherche
max_shine_dalgarno_distance : Position relative de la séquence de Shine-Dalgarno par rapport au codon d’initiation (valeur obtenue en argument du programme)

has_shine_dalgarno  retourne True si un motif de Shine-Dalgarno a été identifié entre max_shine_dalgarno_distance et -6 nucléotides en amont du codon d’initiation et False sinon. La séquence de Shine-Dalgano ne se trouve pas nécessairement dans le même cadre de lecture que le codon d’initiation.

Algorithme de recherche

 predict_genes 
sequence : Séquence du génome
start_regex:  regex object permettant d’identifier un codon start
stop_regex:  regex object permettant d’identifier un codon stop
shine_regex: regex object permettant d’identifier une séquence de Shine-Dalagarno
min_gene_len: longueur minimale d’un gène (valeur obtenue en argument du programme)
max_shine_dalgarno_distance: Position relative de la séquence de Shine-Dalgarno par rapport au codon d’initiation (valeur obtenue en argument du programme)
min_gap: Distance relative minimale entre 2 gènes (valeur obtenue en argument du programme)


Notre programme principal fera appel aux fonctions précédentes nécessaires pour chercher les gènes dans le sens 5’ puis 3’, puis fera appel à reverse complement pour faire cette recherche dans le sens 3’ vers 5’. Aucune modification de vos fonctions ne doit être réalisée pour effectuer cette seconde partie.
La position des gènes prédits dans le sens 3’ vers 5’ devra être corrigée pour qu’elle apparaisse dans le sens 5’ 3’  ! Les positions ont une valeur croissante uniquement.
Vérifiez la qualité de vos prédictions sur jvenn (en copier coller les prodigal/position/predict_genes):
jvenn (inra.fr)

