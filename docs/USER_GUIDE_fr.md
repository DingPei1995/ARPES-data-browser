# ARPES viewer — Guide de l'utilisateur

*[English version](USER_GUIDE_en.md) · [中文版](USER_GUIDE_zh.md)*

*Dernière mise à jour : 28 octobre 2026 — pour ARPES viewer v1.2.0*

La version et sa date figurent dans la barre de titre et en bas de la fenêtre
principale. Signaler un bogue, poser une question, obtenir la dernière version :
<https://github.com/DingPei1995/ARPES-data-browser> (merci d'indiquer le numéro
de version en signalant un bogue).

Ce guide présente les principales fonctions du programme, captures d'écran à
l'appui, dans l'ordre où on les utilise d'habitude : charger les données, les
visualiser, les convertir en impulsion, les analyser et préparer des figures.
Le contexte technique de chaque étape (formules, résultats des tests, choix
de conception) se trouve dans le [README](../README.md) et le
[CHANGELOG](../CHANGELOG.md), en anglais.

Toutes les captures ont été prises sur le programme en fonctionnement, avec
les données synthétiques écrites par `devtools/demo_data.py` : vous pouvez
donc suivre le guide sans temps de faisceau (voir
[Essayer avec des données de démonstration](#2-essayer-avec-des-données-de-démonstration)).

Les noms des boutons et des menus sont donnés tels qu'ils apparaissent dans
le programme, c'est-à-dire en anglais.

---

## Sommaire

1. [Installation et lancement](#1-installation-et-lancement)
2. [Essayer avec des données de démonstration](#2-essayer-avec-des-données-de-démonstration)
3. [La fenêtre principale](#3-la-fenêtre-principale)
4. [Charger des données](#4-charger-des-données)
5. [Ce que les fenêtres de visualisation ont en commun](#5-ce-que-les-fenêtres-de-visualisation-ont-en-commun)
6. [Les cartes : le contour et ses deux coupes](#6-les-cartes--le-contour-et-ses-deux-coupes)
7. [Les coupes : du spectre brut aux paramètres de bande](#7-les-coupes--du-spectre-brut-aux-paramètres-de-bande)
8. [D'angles en impulsion pour une carte, et la zone de Brillouin](#8-dangles-en-impulsion-pour-une-carte-et-la-zone-de-brillouin)
9. [Coupes arbitraires et correction de la surface de Fermi](#9-coupes-arbitraires-et-correction-de-la-surface-de-fermi)
10. [Traitement : lissage, dérivées, courbure](#10-traitement--lissage-dérivées-courbure)
11. [Empilements de courbes et vues 3D](#11-empilements-de-courbes-et-vues-3d)
12. [Arithmétique de coupes : dichroïsme et division par une référence](#12-arithmétique-de-coupes--dichroïsme-et-division-par-une-référence)
13. [Balayages en énergie de photon (kz)](#13-balayages-en-énergie-de-photon-kz)
14. [Supprimer la grille du détecteur](#14-supprimer-la-grille-du-détecteur)
15. [Courbes : EDC, MDC et données résolues en spin](#15-courbes--edc-mdc-et-données-résolues-en-spin)
16. [Balayages résolus spatialement (nano-ARPES)](#16-balayages-résolus-spatialement-nano-arpes)
17. [Figures et exportation](#17-figures-et-exportation)
18. [Sauvegarde automatique, mémoire et récupération](#18-sauvegarde-automatique-mémoire-et-récupération)
19. [Conseils et dépannage](#19-conseils-et-dépannage)

---

## 1. Installation et lancement

Il faut Python 3.9 ou plus récent, et cinq paquets :

```bash
pip install PyQt5 pyqtgraph numpy h5py scipy
```

Puis, depuis le dossier `ARPES_viewer` :

```bash
python ARPES_viewer.py
```

Le programme lit :

| Ligne de lumière / source | Fichiers |
|---|---|
| SOLEIL **ANTARES** | NeXus/HDF5 `*.nxs` (coupes, cartes au déflecteur, balayages spatiaux) |
| SOLEIL **CASSIOPEE** | fichiers texte Scienta SES — spectres isolés, et dossiers numérotés assemblés en carte ou en balayage en énergie de photon |
| **Station de spin** de CASSIOPEE (MBS A-1) | `*.krx` et `*.txt` — coupes, cartes au déflecteur et EDC résolues en spin |
| Ce programme | son propre format `*.nxs` |

## 2. Essayer avec des données de démonstration

```bash
python devtools/demo_data.py            # écrit ./demo_data/*.nxs
```

Six fichiers de données synthétiques sont écrits : une carte de surface de
Fermi, une paire de coupes en polarisation LH et LV, un balayage en énergie
de photon, une EDC résolue en spin, un petit balayage en espace réel et une
carte sur laquelle est imprimée une grille de détecteur. Ils proviennent
d'une structure de bandes modèle (une poche d'électrons hexagonale et une
bande de trous en Γ, des poches de trous en K) avec la coupure de Fermi, la
résolution de l'analyseur et un bruit de Poisson, et s'ouvrent exactement
comme des données mesurées.

Pour refaire toutes les captures de ce guide après avoir modifié une fenêtre :

```bash
python devtools/guide_screenshots.py    # écrit docs/images/*.png
```

---

## 3. La fenêtre principale

![La fenêtre principale](images/main_window_info.png)

La fenêtre principale est un **navigateur** : une ligne par *jeu de
données*, et non par fichier (un fichier `.nxs` contient souvent plusieurs
mesures, chacune sur sa propre ligne). Le type est indiqué après le nom :
`[Cut]`, `[Map]`, `[kz map]`, `[k-map]`, `[SPEM]`, `[EDC]`, `[Spin EDC]`…

Elle a trois parties :

- **Le panneau d'outils** (à gauche) — tous les outils d'analyse et de
  traitement, sous forme de boutons, en deux onglets. **Basic** (affiché au
  démarrage) contient ce par quoi passe presque chaque jeu de données :
  ajustement du front de Fermi, suppression de la grille du détecteur,
  soustraction du fond, les trois conversions en k, la vérification des
  niveaux de cœur, l'alignement des niveaux de Fermi d'une carte kz et la zone
  de Brillouin. **Advanced** contient le reste : composeur de figures,
  lissage et filtres qui affinent les bandes, symétrisation, empilements de
  courbes, ajustements MDC et EDC, arithmétique de coupes, figure en série de
  tranches et les trois vues 3D. Dans chaque onglet, les outils sont groupés
  selon les données qu'ils acceptent (*2D & 3D — cuts, maps*, *2D only —
  cuts*, *3D only — maps* …).
- **La liste** (à droite), avec ses boutons.
- **Le cadre Memory** (en bas à gauche) et la ligne de version (en bas à
  droite).

![L'onglet Advanced du panneau d'outils](images/main_window_advanced.png)

**Les boutons d'outils suivent la sélection.** Un bouton ne s'allume que si
les lignes sélectionnées dans la liste lui conviennent — une conversion en k
(cut) demande une coupe en angle, l'arithmétique de coupes exactement deux
coupes, une vue 3D une carte. Un bouton grisé explique pourquoi au survol.

| Action | Comment |
|---|---|
| Ouvrir un jeu de données dans sa propre fenêtre | **double-clic** sur la ligne (ou clic droit → **Open**) |
| Lancer un outil sur les lignes sélectionnées | son bouton dans le **panneau d'outils** |
| Empiler plusieurs coupes en une carte, ou raccorder une carte / carte kz mesurée en plusieurs morceaux | **Combine...** sous la liste |
| Rogner à un intervalle, regrouper des points, normaliser chaque spectre sur lui-même | **Truncate...**, **Compress...**, **Self-normalise...** sous la liste |
| Ouvrir la coupe le long de la fente ou de l'axe balayé d'une carte | **clic droit** sur la ligne (ou sur une sélection de lignes) |
| Lire les métadonnées, écrire le journal d'expérience, renommer, enregistrer, retirer | le cadre **Datasets** du panneau d'outils, au-dessus de Memory |
| Voir et libérer la mémoire, fermer toutes les fenêtres, ouvrir le dossier de sauvegarde | le cadre **Memory**, en bas à gauche |
| Voir quelle version est utilisée, signaler un bogue, obtenir la dernière version | la barre de titre et la dernière ligne de la fenêtre (lien GitHub) |

Les quatre boutons sous la liste agissent sur les lignes sélectionnées
au-dessus et créent de nouvelles lignes — les originaux ne sont jamais
modifiés.

**Le menu contextuel** ne contient que les façons d'ouvrir les lignes :
**Open**, **Open slit cut** et **Open scan-axis cut** (pour une carte). Une
entrée qui ne convient pas à la sélection est grisée, et son infobulle en
donne la raison.

**Le cadre Datasets** (panneau d'outils, entre les outils et Memory) contient
le reste de ce qui se fait aux lignes : **Experiment info**, **Auto Exp.
log...** (voir [le journal d'expérience](#le-journal-dexpérience-auto-exp-log)),
**Rename...**, **Save datasets...**, **Remove from list** et **Session
log...**. Ses boutons se grisent comme les outils.

### Où se trouvent les choses maintenant (v1.2.0)

Pour qui a utilisé une version antérieure : les outils ont quitté le menu
contextuel et les menus *Functions* des fenêtres pour le panneau d'outils ;
chacun a ainsi une seule entrée, toujours visible.

| Maintenant | Avant |
|---|---|
| Panneau → Basic → **Fermi edge fit...** (coupes et cartes) | fenêtre *Functions → Fermi edge fit* |
| Panneau → Basic → **Remove detector grid...** | fenêtre *Functions → Remove detector grid* |
| Panneau → Basic → **k conversion (cut)** / **(map)**, **kz conversion** | fenêtre *Functions → Convert to momentum* |
| Panneau → Basic → **Align Fermi levels (kz)...** | fenêtre *Functions → Align Fermi levels (kz map)* |
| Panneau → Basic → **Brillouin zone...** | carte k *Functions → Brillouin zone* |
| Panneau → Basic → **Core-level check...** | clic droit → Core-level check |
| Panneau → Basic/Advanced → **Background subtraction**, **Smooth**, **Second derivative**, **Curvature**, **Symmetrise** (une fenêtre chacun) | **Smooth / derivative...** (une fenêtre à onglets) |
| **Truncate...**, **Compress...**, **Self-normalise...** sous la liste | **Crop / bin / normalise...** |
| Panneau → Advanced → **MDC fit (k-cut only)**, **EDC fit** | fenêtre *Functions → MDC / EDC fit* |
| Panneau → Advanced → **Cut arithmetic (two cuts)** | clic droit / fenêtre *Functions → Cut arithmetic* |
| Panneau → Advanced → **Slice series figure...** | carte *Functions → Multi-slice figure* |
| Panneau → Advanced → **3D slices**, **Notched cube**, **Volume projection** | clic droit → 3D view (une fenêtre, trois modes) |
| Panneau → Advanced → **Figure composer...** | clic droit → Open in figure composer |
| Clic droit → **Open scan-axis cut**, bouton **Scan-axis cut** | Open deflector cut, Deflector cut |
| Cadre Memory → **Free memory**, **Close all viewers**, **Open auto-save folder** | sous la liste |
| Cadre Datasets → **Experiment info** | clic droit → Show information |
| Cadre Datasets → **Rename...**, **Save datasets...**, **Remove from list**, **Session log...** | clic droit |
| **Combine...** sous la liste (seulement) | aussi clic droit → Combine N datasets |

**Experiment info** ouvre aussi toutes les métadonnées du jeu de données
dans une fenêtre à part, avec un champ de filtre et **Export metadata
(.csv)** :

![Informations sur le jeu de données](images/information_window.png)

### Le journal d'expérience (Auto Exp. log)

**Datasets → Auto Exp. log...** écrit le journal du temps de faisceau -- un
fichier CSV, une ligne par mesure -- à partir des métadonnées des fichiers
eux-mêmes.

![La fenêtre Auto Exp. log](images/exp_log_dialog.png)

1. **Log file.** Choisissez le dossier et tapez le nom du fichier. Si un
   journal de ce nom existe déjà, les nouvelles lignes sont écrites **après
   sa dernière ligne** et rien au-dessus n'est modifié ; sinon il est créé,
   avec les colonnes pour en-tête. (Si les colonnes ont changé depuis, une
   ligne marquée `# columns` avec les nouveaux titres vient d'abord.)
2. **Columns.** *Fields of:* choisit un jeu de données de la liste -- du
   type que vous consignez le plus -- et **Read fields** liste tout ce qu'il
   contient, avec sa valeur : nom du fichier, type, heure de début et
   intervalles balayés tirés des données, puis chaque champ de métadonnées
   et chaque moteur. **Cochez** un champ pour en faire une colonne. Dans le
   tableau de droite, double-cliquez un titre pour le renommer, et rangez
   les colonnes dans l'ordre voulu avec **Top / Up / Down / Bottom**. **Add
   empty column** ajoute une colonne à remplir à la main (commentaires).
   *Ranges as:* écrit un intervalle `start/stop/step` ou `start:step:stop`.
   Les colonnes sont retenues pour la session suivante ; **Save columns
   as...** / **Load columns...** en gardent un jeu par ligne de lumière.
3. **If scanned.** Ce qu'écrit une colonne de position ou d'angle dépend de
   la mesure :
   - une **coupe** : l'échantillon n'a pas bougé -- la valeur du moteur ;
   - un balayage **SPEM** : l'axe de platine balayé, `start/stop/step`
     (*SPEM: scanned X* / *Y* pour les colonnes de la platine) ;
   - une **carte** prise au manipulateur plutôt qu'au déflecteur : l'angle
     en `start/stop/step` (*map: scanned angle*), et une position de
     l'échantillon corrigée avec l'angle en `début..fin` ;
   - une **carte kz** : l'intervalle d'énergie de photon (*kz map: photon
     energy*).

   *auto* (par défaut) reconnaît le moteur balayé à son nom, et une position
   qui a bougé à ses valeurs. Toute autre mesure reçoit la valeur du champ.
4. **Record.** *The selected rows* ou *Every dataset in the list*, puis
   **Write records**. *Skip datasets already in the log* omet les fichiers
   dont le nom figure déjà dans la colonne des noms. **Keep recording**
   écrit la ligne de chaque fichier ajouté ensuite à la liste, même fenêtre
   fermée : pendant le temps de faisceau, ajoutez chaque fichier mesuré et
   le journal se tient à jour.

Si le journal est ouvert dans Excel, Windows refuse de l'écrire : fermez-le
dans Excel et recommencez.

**Clear list** demande confirmation, et propose d'enregistrer ce qui a été
calculé pendant la session sans avoir encore été sauvegardé. Les fichiers sur
le disque ne sont jamais supprimés.

### Combiner des jeux de données : empiler et raccorder

**Combine...** (sous la liste) fait deux choses différentes, et demande
laquelle :

- **Stack (empiler)** — *coupes seulement*. Chaque coupe devient une tranche
  d'une nouvelle carte, le long d'un axe que vous nommez : série en
  température, en tension de grille, une coupe par énergie de photon. Le
  résultat s'examine comme n'importe quelle carte.
- **Join (raccorder)** — *coupes, cartes, cartes kz ou cartes k, toutes du
  même type*. Les morceaux d'une mesure prise en plusieurs fois sont remis
  bout à bout le long de l'un de leurs propres axes : une carte dont la plage
  du déflecteur a été balayée en deux parties, une carte kz dont la plage
  d'énergie de photon a été répartie sur plusieurs acquisitions, une coupe
  dont la fenêtre d'énergie a été prise en deux moitiés.

**Empiler des coupes**

1. Sélectionnez les coupes dans la liste (Ctrl/Maj-clic) et cliquez sur
   **Combine...**. L'empilement exige les mêmes axes point par point
   (même nombre de points, mêmes valeurs, mêmes unités) ; sinon *Stack* est
   grisé avec la première différence, et *Join* est proposé à la place.
2. Il demande s'il faut les combiner **dans l'ordre de la liste**. *Yes* garde
   cet ordre ; *No* permet de le choisir.
3. Choisissez **Stack**, puis définissez le nouvel axe : **New axis is**
   (température, tension de grille, énergie de photon, délai, position,
   autre — ce choix remplit *Name* et *Unit*, modifiables), **Start** et
   **Step** (le pas). Chaque coupe est affichée avec la valeur qu'elle reçoit.
   Si vous avez choisi de réordonner, **faites glisser** les lignes (ou
   utilisez **Move up** / **Move down** / **Reverse**). Cliquez sur
   **Combine**.

![Combiner des coupes](images/combine_cuts_dialog.png)

Une série en énergie de photon de coupes en angle devient une `kz map`, que
*Align Fermi levels (kz map)* et *Convert to momentum (kz, k∥)* acceptent
(voir [§13](#13-balayages-en-énergie-de-photon-kz)) ; tout autre axe donne une
`map` dont le premier axe est marqué comme n'étant pas un angle, si bien
qu'aucune conversion en k n'est proposée.

**Raccorder les morceaux d'une mesure**

1. Sélectionnez les morceaux (toutes des cartes, ou toutes des cartes kz, ou
   toutes des coupes …) et cliquez sur **Combine...**. Les cartes et
   cartes kz passent directement à *Join*.
2. **Join along** : l'axe selon lequel la mesure a été découpée. Il part de
   l'axe dont l'étendue diffère d'un morceau à l'autre — l'axe déflecteur
   d'une carte, l'énergie de photon d'une carte kz, ou l'énergie de coupes
   prises en deux fenêtres d'énergie — et il est donc en général déjà juste.
3. **Where pieces overlap** : *Average them* (moyenne), ou garder le morceau
   le plus tôt / le plus tard dans l'ordre affiché (faites glisser les lignes
   pour le changer).
4. La boîte de dialogue donne l'intervalle de chaque morceau et le
   **Result** — l'intervalle raccordé et le nombre de points — avec des notes
   sur ce qu'il a fallu faire. Cliquez sur **Combine**.

![Raccorder une carte mesurée en deux parties](images/combine_dialog.png)

Traitement des axes :

- Des morceaux sur une même grille (même pas, intervalles qui se touchent ou
  se recouvrent d'un nombre entier de points — le cas habituel) sont
  **copiés tels quels**, sans interpolation.
- Un morceau de pas différent ou décalé d'un demi-pas est **interpolé** sur
  l'axe commun (au pas du premier morceau), et la boîte de dialogue le dit.
- Si les **autres axes** diffèrent — une fenêtre en énergie décalée de
  quelques meV entre deux acquisitions, un nombre différent de canaux de
  fente — tous les morceaux sont rééchantillonnés sur l'intervalle que
  **tous** couvrent, et la boîte de dialogue le dit.
- Un intervalle couvert par aucun morceau est un **trou**, laissé à zéro et
  signalé.
- Une carte kz raccordée reste une `kz map` : lancez *Align Fermi levels (kz
  map)* sur le résultat plutôt que sur les morceaux (les ajustements par
  énergie de photon des morceaux ne sont pas repris).

---

## 4. Charger des données

**Load data...** ouvre la fenêtre de chargement. **Quick add...** la saute et
ajoute les fichiers choisis avec la détection automatique.

![La fenêtre de chargement](images/loader_dialog.png)

1. **Files** — choisir un ou plusieurs fichiers (*Select file(s)...*). Le
   filtre de fichiers est construit à partir de ce que les lecteurs installés
   savent ouvrir.
2. **Reader** — détecté automatiquement, et le résultat est affiché
   (*Detected: 3 × This program's own format*). À ne changer que si un
   fichier est mal lu.
3. **Datasets found** — le contenu de chaque fichier, lu à partir de sa seule
   *structure*. C'est instantané même pour un balayage de plusieurs Go ; rien
   n'est lu avant d'ouvrir une ligne.
4. **Axes** — *Axis order* réordonne une fois pour toutes les dimensions du
   tableau, dès le chargement. *First axis of a map is* indique ce qui a été
   balayé : un angle de déflecteur (par défaut), une **énergie de photon**, une
   température, une tension de grille, un délai… L'unité reste ainsi correcte
   partout ensuite, et aucune conversion en k n'est proposée pour un axe qui
   n'est pas un angle.

**Dossiers CASSIOPEE.** Sélectionnez n'importe quel fichier d'une série
numérotée (`nom_1_ROI1_.txt`, `nom_2_ROI1_.txt`, …) : tout le dossier est
listé une seule fois, comme un cube assemblé. Qu'il devienne une carte en
angle ou une carte kz est déterminé à partir des fichiers de paramètres `_i`
placés à côté des spectres. La ligne porte l'énergie de photon et la
polarisation lues dans ces fichiers — `FineMap 117eV LH`, ou
`KzScan 40-120eV LH` pour un balayage en énergie de photon — pour distinguer
les répétitions d'un même balayage à une autre énergie ou polarisation.

**Fichiers ANTARES en énergie de photon (PhD).** Un fichier `.nxs` PhD (une
entrée `Phd_NNNN` par énergie de photon) est listé comme une seule `kz map` :
les images sont triées par énergie de photon et chaque axe en énergie
cinétique est rapporté à E_F avec son propre hν (travail de sortie 4,5 eV,
que la conversion k_z utilise ensuite). Chaque image s'ouvre encore seule
comme une coupe.

### Types de données

| Type | Axes | Ce que c'est |
|---|---|---|
| `cut` | angle (ou k) × E | un spectre |
| `map` | angle × angle × E | un balayage au déflecteur ou en angle polaire |
| `k_map` | kx × ky × E | une carte convertie en impulsion |
| `kz_map` | hν × angle × E | un balayage en énergie de photon |
| `kz_map_k` | k_z × k∥ × E | ce balayage converti en impulsion |
| `spem_4d`, `spem_1d` | x (× y) × angle × E | balayages en espace réel |
| `edc`, `mdc`, `spin_edc` | E (ou k) × canaux | courbes, avec canaux σ facultatifs |

---

## 5. Ce que les fenêtres de visualisation ont en commun

Chaque jeu de données s'ouvre dans sa propre fenêtre : on peut placer côte à
côte une carte, ses coupes et un autre fichier.

![Une fenêtre de coupe avec le curseur de lecture](images/cut_viewer.png)

- **La ligne du haut** : le menu **Functions** (tout ce que la fenêtre sait
  faire), **Plot tools** (lignes guides et étiquettes sur ce graphique — voir
  plus bas), la palette de couleurs (propre à chaque fenêtre) avec **Flip**, et
  les réglages de vue — étendue des axes, **Auto**, **Inv** (inverser),
  **Reset**, **Grid**.
- **Au-dessus de chaque image** : la barre **Levels** (faire glisser les deux
  poignées, ou saisir des centiles **Min%/Max%** ; **Clip** les applique) et
  **Gamma**. Un min/max que vous fixez est conservé quand la tranche se
  déplace (l'énergie d'un contour, l'angle d'une coupe, le curseur d'une
  carte spatiale), si bien que toutes les positions partagent la même échelle
  de couleurs ; ramenez les poignées aux deux extrémités, ou réglez 0 %–100 %,
  pour suivre de nouveau la plage propre à chaque image.
- **Sous chaque image** : **Smooth** (affichage seulement, les données ne
  changent pas) et **Lock x/y** (rapport d'aspect).
- **Curseur de lecture** — clic droit sur l'image → *Readout cursor*. Il
  donne la valeur sous le curseur et, sur une coupe, trace l'**EDC** (à
  droite) et la **MDC** (en dessous) passant par ce point. **EDC ±** /
  **MDC ±** intègrent sur une fenêtre, affichée comme une bande ombrée sur
  l'image. Saisissez une position dans la ligne **Cursor** et appuyez sur
  Entrée pour l'y placer exactement.
- **EDC → list** / **MDC → list** ajoutent la courbe affichée à la liste
  principale comme jeu de données à part entière (voir
  [§15](#15-courbes--edc-mdc-et-données-résolues-en-spin)).
- **Clic droit sur l'image → Export** enregistre ce que montre le panneau —
  données (`.nxs`), image, ou axes vectoriels (`.eps`) ; voir
  [§17](#17-figures-et-exportation).
- **Functions → Slice → Save slice to the main list / Freeze slice in a new
  window** conserve une copie de la tranche affichée. La fenêtre figée (et
  *Open this position in a new window* sur une carte spatiale) s'ouvre avec le
  min/max et le gamma du panneau, pour une comparaison directe. Elle a ses
  propres **Plot tools**, et un bouton **Figure composer...** pour en faire
  une figure ([§17](#17-figures-et-exportation)).

Le menu **Functions** ne contient plus que ce qui agit sur la fenêtre
elle-même — les outils qui prennent un jeu de données sont dans le panneau
d'outils de la fenêtre principale ([§3](#3-la-fenêtre-principale)) :

| Fenêtre | Functions |
|---|---|
| carte / carte kz / carte k (contour) | Arbitrary cut ; Slice : enregistrer / figer la tranche |
| coupe | Straighten bands (FS correction) ; Slice : enregistrer / figer |
| coupe selon la fente d'une carte | Straighten bands (FS correction) pour tout le cube ; Show the map |
| courbe | Curve fit, Spin analysis ; Crop, Bin, Normalise, Background, Shift axis, Add √N errors |

Une entrée qui s'applique mais ne peut pas encore être lancée est grisée, et
son infobulle en donne la raison.

**Case Brillouin zone.** Une carte k ou une carte kz convertie sur laquelle
une zone a été tracée ([§8.2](#82-zone-de-brillouin)) la garde : sa fenêtre a
une case **Brillouin zone** sur la ligne du haut pour l'afficher ou la
masquer.

## 6. Les cartes : le contour et ses deux coupes

Une carte s'ouvre sur son **contour à énergie constante**. Le curseur du haut
choisit l'énergie (**Pos** en eV, ou **Ind** en indice ; **±** intègre sur une
fenêtre, en **mean** ou **sum**). **Scan-axis cut** et **Slit cut** ouvrent
les deux coupes orthogonales, chacune dans sa fenêtre. (L'axe de balayage
est ce qui a été balayé pour faire la carte — en général l'angle du
déflecteur, mais aussi une rotation du manipulateur ou, pour une carte kz,
l'énergie de photon — d'où le nom du bouton.)

![Contour à énergie constante d'une carte](images/map_contour.png)

| Coupe selon la fente (angle de fente vs énergie) | Coupe selon le déflecteur (angle de déflecteur vs énergie) |
|---|---|
| ![Slit cut](images/map_slit_cut.png) | ![Scan-axis cut](images/map_scan_axis_cut.png) |

**Un seul curseur pour les trois fenêtres.** Activer le curseur de lecture
dans l'une l'active dans les autres, et le déplacer n'importe où recoupe les
autres en ce nouveau point. Chaque axe garde sa couleur partout :
**déflecteur en rouge, fente en vert, énergie en bleu**. Une ligne de curseur
a la couleur de l'axe qu'elle maintient constant, et l'EDC ou la MDC qu'elle
alimente est tracée de la même couleur. Les largeurs d'intégration (les
champs ±) appartiennent à chaque fenêtre et ne sont jamais partagées.

On peut aussi ouvrir directement une coupe d'une carte depuis la liste (clic
droit → **Open slit cut** / **Open scan-axis cut**) ; le contour reste alors
caché jusqu'à ce qu'on le demande (*Functions → Show the map*).

---

## 7. Les coupes : du spectre brut aux paramètres de bande

L'ordre habituel pour une coupe est : **ajustement du front de Fermi → décaler
l'axe en énergie → conversion en impulsion → ajustement MDC/EDC**. Chaque
étape est un bouton du panneau d'outils : sélectionnez la coupe dans la
liste et appuyez dessus.

### 7.1 Ajustement du front de Fermi

Panneau → Basic → **Fermi edge fit...** ajuste un front de Fermi (Fermi–Dirac
× densité d'états linéaire, convoluée par une résolution gaussienne) sur
l'EDC sommée dans un cadre. Il ouvre la fenêtre de la coupe avec la fenêtre
d'ajustement ; les deux vont ensemble, et fermer l'une ferme l'autre.

![Fenêtre du niveau de Fermi](images/fermi_dialog.png)

1. Appuyez sur **Show the selection box** et faites glisser le cadre sur le
   front, sur la coupe (ou saisissez les intervalles d'angle et d'énergie).
   L'intervalle d'angle fixe les canaux sommés ; l'intervalle d'énergie est
   la fenêtre d'ajustement.

   ![Le cadre d'ajustement sur la coupe](images/fermi_cut_box.png)

2. **Estimate**, puis **Fit**. Cochez **Hold** pour fixer un paramètre (la
   température est fixée par défaut à la valeur enregistrée). Chaque
   paramètre est donné avec son incertitude.
3. **Offset energy axis** ajoute à la liste une copie de la coupe en
   **E − E_F**.

Le coin *Advanced* divise par la coupure de Fermi (pour les gaps et les
échantillons de référence) ou ajuste E_F canal par canal.

**Sur une carte.** Sélectionnez une carte (en angle ou en k) et appuyez sur
le même bouton : la **coupe selon la fente** de la carte s'ouvre avec la
fenêtre d'ajustement, et son curseur choisit la coupe ajustée. Après un
ajustement, **Offset energy axis** demande s'il faut décaler **toute la carte
par l'E_F de cette coupe** (*Yes*) ou ajuster **un E_F par coupe** (*No*).
Dans ce second cas, la liste **Fitted cuts** à droite se remplit au fur et à
mesure : déplacez le curseur, **Fit**, recommencez. Un second ajustement au
même endroit remplace le premier, et chaque E_F de la liste se modifie à la
main. Les coupes entre deux coupes ajustées sont décalées selon une droite
entre leurs voisines (au-delà de la première et de la dernière, selon
l'ajustement le plus proche), tracée sous la liste. **Offset energy axis by
these fits** ajoute la carte corrigée à la liste.

![Ajustement de Fermi sur une carte, un E_F par coupe](images/fermi_map_cuts.png)

### 7.2 Conversion en k d'une coupe

Panneau → Basic → **k conversion (cut)...** convertit l'axe angulaire en k∥
(Å⁻¹). Une **ligne verticale** est tracée sur la coupe : faites-la glisser
sur l'angle de fente de Γ. La coupe est alors supposée passer par Γ — ce
n'est pas exact si elle le manque, mais c'est le seul choix quand seule cette
coupe a été mesurée.

![La ligne Γ sur la coupe](images/cut_k_line.png)

![Conversion en k d'une coupe](images/cut_k_dialog.png)

- **Advanced: Γ off this cut** est pour une coupe qui ne passe pas par Γ :
  saisissez les angles déflecteur et fente de Γ, ou **From a converted
  map...** pour les reprendre d'une carte du même échantillon déjà convertie.
- Si la coupe est déjà en E − E_F, **Energy offset** est rempli avec le E_F
  soustrait (la note en dessous le signale) : la conversion utilise les
  vraies énergies cinétiques, et le résultat reste en E − E_F.
- **Réglages mémorisés.** La fenêtre part de la dernière conversion d'une
  coupe de même taille (même nombre de points sur chaque axe) — le cas
  habituel pour une série de coupes d'un même échantillon. Les conversions
  de cartes et kz font de même.

![La coupe en impulsion](images/cut_in_k.png)

### 7.3 Ajustement MDC et ajustement EDC

Deux outils de l'onglet Advanced ajustent des pics ligne par ligne :

- **MDC fit (k-cut only)...** ajuste les courbes de distribution en impulsion
  et en déduit la bande à partir des centres ajustés : vitesse de Fermi,
  masse effective, self-énergie. N'accepte qu'une coupe dont l'axe est déjà
  en Å⁻¹ — une vitesse en eV·deg n'aurait aucun sens — et son bouton est grisé
  sinon.
- **EDC fit...** ajuste les courbes de distribution en énergie, sur
  n'importe quelle coupe, éventuellement avec la coupure de Fermi dans le
  modèle. v_F et m\* de la bande ajustée demandent une coupe en impulsion.

![Panneau d'ajustement MDC](images/mdc_fit_panel.png)

**Guided fit (step by step)**, en haut du panneau, guide tout l'ajustement
pas à pas à côté des réglages, en disant à quoi sert chaque étape et en
encadrant la partie du panneau qu'elle utilise : ce que fait un ajustement
MDC (EDC), la forme de pic, l'intervalle (*Start here* / *End here* à la
ligne bleue), les bandes, les points de départ (le guide place la ligne en
cinq endroits de l'intervalle — cliquez à chaque fois le pic de la bande sur
la courbe), une ligne d'essai, la série et le résultat. *Next* attend que
l'étape ait ce qu'il faut à l'ajustement.

Les images et les réglages se partagent la fenêtre avec une poignée entre
eux : les réglages restent visibles ; la molette sur une image fait défiler
la fenêtre quand elle est trop courte (Ctrl + molette zoome l'image).

À la main :

1. Onglet **Bands** : **Add band** pour chaque branche ; sélectionnez une
   bande, appuyez sur **Place a seed** et cliquez sur le pic, sur la courbe (à
   droite). Quatre ou cinq points de départ le long d'une bande suffisent ;
   **Line at** change la ligne examinée.
2. **Try this line** montre l'ajustement de la ligne courante (en tirets :
   l'estimation initiale).
3. Onglet **Fit settings** : forme de pic (lorentzienne, gaussienne, Voigt),
   fond, pas, puis **Fit the series**. Les positions ajustées sont tracées
   sur l'image, avec une taille proportionnelle à l'amplitude.
4. **Dispersion...** ajuste la bande : **v_F** (droite), **m\*** (parabole) —
   et vérifie jusqu'où l'ajustement reste stable. **Self-energy...** donne
   Re Σ et Im Σ.

![Fenêtre de dispersion](images/dispersion_window.png)

**Save table...** écrit toutes les lignes ajustées ; **Band to list** ajoute
une bande à la liste principale.

---

### 7.4 Vérification des niveaux de cœur

Un balayage large des niveaux de cœur est lui aussi listé comme une coupe.
Sélectionnez-le et appuyez sur panneau → Basic → **Core-level check...** pour
savoir d'où viennent ses pics.

1. Saisir les éléments de l'échantillon — `NbSCl` ou `Nb, S, Cl` — puis
   **Next**. *Match within* (0,5 eV) est l'écart maximal entre un pic et
   l'énergie de liaison tabulée pour l'attribuer à l'échantillon.

   ![Éléments](images/core_level_elements.png)

2. Vérifier la référence en énergie. E_F est ajusté si la fenêtre contient
   le front de Fermi, sinon placé à hν − W. Corriger l'énergie de photon,
   le travail de sortie ou E_F s'ils sont faux (la ligne verte se déplace à
   la souris) ; un hν − E_F loin de 4–5 eV signale une énergie de photon
   erronée. **Analyse**.

   ![Référence en énergie](images/core_level_reference.png)

3. Les résultats : l'EDC en énergie de liaison, chaque pic marqué — vert
   *sample* (orbitale nommée), orange *uncertain* (un niveau de
   l'échantillon à moins de 2 eV : déplacement chimique ou charge — à vous
   de juger), rouge *not the sample* (éléments courants qui pourraient
   l'expliquer : contamination, lumière d'ordre supérieur, raies Auger).
   Cliquer une ligne affiche tous les candidats ; la colonne *Assignment*
   accepte votre propre attribution.

   ![Vérification des niveaux de cœur](images/core_level_check.png)

**Raies faibles et pics à la main.** Cochez **Log counts** pour voir les raies
faibles à côté d'une raie forte ; **Min height** (0,2 % de la plus forte par
défaut) fixe la hauteur minimale d'un pic trouvé. Pour ajouter un pic que la
recherche a manqué, appuyez sur **Click to add peak** puis cliquez sur le
graphe (un double-clic fonctionne toujours), ou saisissez son énergie de
liaison après *or at* et cliquez sur **Add** : sa position est affinée par un
ajustement local, il est marqué ◆ et nommé avec tous les niveaux tabulés à
±2 eV. **Remove selected** retire des pics ; **Reset peaks** rétablit les pics
trouvés et oublie ceux ajoutés à la main. Les attributions saisies survivent
à **Re-analyse**.

**Peak table...** liste tous les niveaux tabulés, filtrés par élément et par
plage d'énergie de liaison : chercher un pic, afficher un niveau sur le
graphe, ou l'ajouter comme pic avec ce niveau pour attribution. **Select
range** affiche une bande à faire glisser sur une région ; **Zoom window...**
l'ouvre (ou la vue courante) dans une fenêtre à part, avec un curseur, un
ajustement gaussien au curseur et la dérivée seconde négative divisée par
l'intensité — un épaulement ou une raie faible à côté d'une raie forte y est
aussi visible que la raie forte.

Si tous les pics sont décalés de la même quantité, les notes le disent et
**Apply suggested shift** corrige (échantillon chargé, ou E_F faux).
**Add EDC to list** ajoute l'EDC en énergie de liaison à la liste ;
**Copy report** copie le résultat en texte. Une mesure à une seconde
énergie de photon tranche ce qui reste douteux : les niveaux de cœur
restent fixes en énergie de liaison, les raies Auger et la lumière d'ordre
supérieur se déplacent.

## 8. D'angles en impulsion pour une carte, et la zone de Brillouin

### 8.1 Conversion en k d'une carte

Sélectionnez une carte en angle et appuyez sur panneau → Basic → **k
conversion (map)...**. La carte s'ouvre sur un contour nu — sans coupes ni
autres fonctions, juste l'image — à côté de la fenêtre de conversion ; les
deux se ferment ensemble. Faites glisser le curseur de lecture sur le point
qui doit devenir **k = (0, 0)** — la fenêtre le suit — puis **Convert**. Le
résultat est une nouvelle ligne `k-map` dans la liste.

| Le contour pour choisir Γ | La conversion |
|---|---|
| ![Contour de conversion](images/map_k_contour.png) | ![Conversion en k d'une carte](images/map_k_dialog.png) |

- **Γ deflector angle / Γ slit angle** : les angles du point qui devient
  k = (0, 0), remplis depuis le curseur.
- **Sample rotation** fait tourner la carte autour de la normale à la
  surface ; **Set rotation from contour...** la déduit d'une direction
  cliquée sur le contour (deux points le long d'une ligne de symétrie, par
  exemple).
- **Output grid** : nombre de points, ou résolution (Å⁻¹, meV). Les nombres
  de points partent de ceux de la mesure — autant de points kx que d'angles
  déflecteur, de points ky que d'angles de fente, de points E que d'énergies.
- **Energy offset** est ajouté à l'axe en énergie avant la conversion, pour
  un axe qui n'est pas une vraie énergie cinétique (rempli automatiquement
  pour une carte en E − E_F).
- Les réglages partent de la dernière conversion d'une carte de même
  taille : une deuxième carte d'une série n'a rien à ressaisir.

![La carte convertie en k](images/k_map_contour.png)

### 8.2 Zone de Brillouin

Sélectionnez une carte k (ou une carte kz convertie) et appuyez sur panneau
→ Basic → **Brillouin zone...**. Comme la conversion de carte, il ouvre un
contour à lui avec la fenêtre de la zone, et les deux se ferment ensemble.

![Fenêtre de la zone de Brillouin](images/brillouin_dialog.png)

1. D'abord le **groupe d'espace** (1–230) ; les paramètres de maille qu'il
   impose sont grisés et suivent les paramètres libres. Ou passez en **2D
   layer** pour un réseau de surface.
2. **Cut plane** : la normale à la surface (hkl) et un décalage le long de
   celle-ci — on trace la section de la zone 3D par ce plan, pour n'importe
   quelle surface et pas seulement (001). **Through point** donne les
   décalages pour lesquels le plan passe par un point de haute symétrie (Γ,
   Z, L, F ... pour un cristal rhomboédrique) : choisissez-en un plutôt que
   de calculer la valeur. **Show cut plane in 3D...** la prévisualise, avec
   les points nommés et, entourés, ceux que contient le plan ; **Label
   high-symmetry points** les écrit sur le contour, ce qui aide à régler
   l'azimut. La fenêtre s'ouvre sur le dernier cristal utilisé.
3. Zone **Conventional** ou **Irreducible** ; **Azimuth** (ou **Pick
   direction from contour...**) ; **Repeat across the field of view**.
4. **Moiré** ajoute une seconde couche, tournée ou désaccordée, et sa zone de
   moiré.
5. **Plot**.

La fenêtre est haute : son contenu défile, et elle se redimensionne comme
toutes les autres.

**La zone reste avec le jeu de données.** Une fois tracée, elle est une
propriété de la carte k : elle est enregistrée avec elle dans un fichier
`.nxs`, et sa fenêtre l'affiche avec la case **Brillouin zone** de la ligne
du haut. **Clear overlay** la retire.

![Carte k avec la zone affichée par sa case](images/k_map_with_bz.png)

Les exportations peuvent aussi la tracer : **Export → Image...** propose
**Include the Brillouin zone**, et **Axes only (.eps)...** propose **Draw the
Brillouin zone in the box**, pour obtenir la zone en traits vectoriels sur
l'image dans Illustrator.

![Exportation avec la zone](images/export_dialog_bz.png)

---

## 9. Coupes arbitraires et correction de la surface de Fermi

### 9.1 Coupe arbitraire

*Functions → Arbitrary cut...* sur n'importe quel contour : cliquez jusqu'à
six points (ou saisissez-les) ; la coupe les parcourt dans l'ordre — ici
Γ → M → K → Γ sur la carte k. **One dataset per segment** crée un jeu de
données par segment.

| Chemin sur le contour | Points |
|---|---|
| ![Chemin de la coupe arbitraire](images/arbitrary_cut_path.png) | ![Fenêtre de la coupe arbitraire](images/arbitrary_cut_dialog.png) |

Le résultat est une coupe en fonction de la distance le long du chemin ; les
coins du chemin sont marqués par des tirets verts, pour qu'un coude à cet
endroit ne soit pas pris pour de la physique.

![La coupe le long de Γ–M–K–Γ](images/arbitrary_cut_result.png)

### 9.2 Correction de la surface de Fermi (FS)

*Functions → Straighten bands (FS correction)...* sur une coupe (ou sur la coupe selon la fente
d'une carte) redresse une structure qui devrait être plate — en général le
front de Fermi, courbé en fonction de l'angle de l'analyseur. Cliquez le long
de la structure : un polynôme (d'ordre 2 par défaut) est ajusté et chaque
colonne d'angle est décalée en énergie. **From a reference (Au)...** prend
plutôt la courbe sur un spectre d'or présent dans la liste. Sur la coupe
selon la fente d'une carte, le même décalage est appliqué à tout le cube.

![Correction FS](images/fs_correction_dialog.png)

---

## 10. Traitement : lissage, dérivées, courbure

Chaque filtre a sa propre fenêtre, ouverte depuis le panneau d'outils avec
une ou plusieurs coupes (ou cartes) sélectionnées :

| Onglet | Bouton | Ce que cela fait |
|---|---|---|
| Basic | **Background subtraction...** | retire l'intensité qui n'est pas celle des bandes : fond inélastique de Shirley ou de Tougaard, un plancher commun à tous les angles, un polynôme lisse — ou divise par la coupure de Fermi |
| Advanced | **Smooth...** | lissage gaussien, Savitzky–Golay ou moyenne glissante, avec suppression des pics facultative |
| Advanced | **Second derivative...** | −d²I/dE² ou −d²I/dk² après un pré-lissage, pour affiner les bandes |
| Advanced | **Curvature...** | courbure 1D ou 2D (Zhang et al.), qui suit la bande mieux qu'une dérivée seconde |
| Advanced | **Symmetrise...** | miroir d'une coupe par rapport à une droite ; pour une carte, miroir ou rotation de ses plans à énergie constante |

Chaque résultat est un nouveau jeu de données qui garde la trace des étapes
qui l'ont produit.

![La fenêtre de courbure](images/process_dialog.png)

- **Source** et **Result** sont affichés côte à côte, sur la même échelle,
  pour voir quand un réglage a « mangé » la bande. Ils se déplacent et
  zooment ensemble.
- **Lisser d'abord.** Les fenêtres de dérivée seconde et de courbure
  s'ouvrent avec un conseil : toutes deux amplifient le bruit, et marchent
  donc mieux sur des données déjà lissées — lancez **Smooth...**,
  sélectionnez son résultat dans la liste et ouvrez le filtre sur celui-ci,
  ou augmentez les largeurs de pré-lissage dans la fenêtre même. Un résultat
  lissé porte l'étiquette **‹smoothed›** dans la liste (tout ce qui en est
  tiré aussi) ; sur ces données le conseil laisse place à une ligne « déjà
  lissé » et les largeurs de pré-lissage partent de 0.
- **Shape** : *Square pixels* (un point de données est carré — sans
  distorsion pour une coupe angle-énergie), *Equal axes* (une unité a la même
  longueur sur les deux axes — adapté à angle-angle ou k-k ; choisi
  automatiquement quand les deux axes ont la même unité), ou *Fill the box*.
- **Live preview** recalcule à chaque changement de valeur ; **Apply** ajoute
  le résultat à la liste.
- Valeurs de départ : **Symmetrise** fait un miroir par rapport à une droite
  verticale passant par angle (ou k) = 0, c.-à-d. k → −k ; la **Fermi
  cut-off** de *Background subtraction* part de l'E_F enregistré dans le
  fichier, sinon de la chute la plus raide de l'EDC intégrée en angle, avec
  la température et la résolution du fichier.

Pour une carte, la fenêtre travaille plan par plan et demande d'abord **sur
quels plans** — contours à énergie constante ou dispersions — car la courbure
d'un contour et celle d'une coupe sont des grandeurs différentes. **Preview
at** choisit le plan affiché ; **Apply to the whole cube** traite tous les
plans. Les opérations et les largeurs portent le nom des axes du plan choisi
(*Second derivative along Kinetic energy*…), et les largeurs de lissage
partent de 1/100 de chaque axe.

![Lissage d'une carte](images/volume_process_dialog.png)

**Truncate...**, **Compress...** et **Self-normalise...** (sous la liste)
changent la grille ou l'échelle d'un ou plusieurs jeux de données à la fois
— à condition qu'ils soient du même type et aient le même nombre de points
sur chaque axe : **Truncate** rogne à un intervalle, **Compress** regroupe
des points, **Self-normalise** divise chaque spectre par son propre total.

**Self-normalise** divise chaque morceau des données par son propre total (ou
son maximum) :

- **Normalise** choisit le morceau en termes d'ARPES et part du choix
  habituel : *each EDC* pour une coupe, *each cut (slit × energy at one scan
  step)* pour une carte, une carte k ou une carte kz — ce qui retire une
  dérive du flux ou du courant faisceau le long du balayage. Autres choix :
  chaque MDC, chaque EDC d'une carte, chaque carte à énergie constante, ou
  *custom*.
- **Summed over** montre les axes sur lesquels on somme ; les cases peuvent
  être changées à la main, et la phrase en dessous dit ce qu'elles font
  (« Each plane over Angle, slit × Kinetic energy is divided by its own
  total — one factor for every Angle, deflector »).
- **Sum only over** restreint la somme à une fenêtre sur l'un des axes
  sommés — l'énergie par défaut — par exemple les énergies d'un niveau de
  cœur.

![Auto-normalisation](images/data_operations.png)

---

## 11. Empilements de courbes et vues 3D

**Stack plot (EDC / MDC)...** (panneau → Advanced, une coupe sélectionnée) :
une cascade d'EDC ou de MDC, avec le nombre de courbes, le regroupement, la
normalisation, le décalage, le cisaillement et la coloration ; **Send to
figure composer** en fait un panneau de figure.

![Empilement de courbes](images/stack_plot.png)

Trois outils de l'onglet Advanced montrent une carte en 3D, chacun dans sa
propre fenêtre — faites glisser l'image pour la faire tourner :

- **3D slices...** — les trois plans orthogonaux du cube, chacun placé où
  vous voulez (case, curseur ou position saisie pour chaque axe).
- **Notched cube...** — le cube avec un coin retiré : une carte à énergie
  constante et deux dispersions dans une seule image.
- **Volume projection...** — tout le volume projeté selon la direction de
  vue : intensité maximale, somme, composition alpha, ou isosurface ombrée.

| 3D slices | Notched cube |
|---|---|
| ![3D slices](images/volume_slices.png) | ![Notched cube](images/volume_view.png) |

- **Échantillonnage.** Chaque image rééchantillonne le cube ; chaque fenêtre
  s'ouvre donc en échantillonnage *Balanced*, qui suit la souris.
  L'échantillonnage courant est indiqué sous **Save the view...** ; **Change
  sampling...** choisit la région et *Fast*, *Balanced* ou *Fine*, puis
  redessine depuis le cube complet.

  ![Changer l'échantillonnage](images/volume_range_dialog.png)

- **Couleurs.** *Apply to* règle les couleurs de **toutes les tranches** à
  la fois ou d'une seule orientation (le plan en énergie dans une palette,
  les dispersions dans une autre) : palette (et *reversed*), niveaux
  (**Auto**, ou **Each its own** pour une échelle par tranche), gamma,
  **Opacity** et **Fade low values**, qui rend transparente l'intensité
  faible pour laisser voir les plans de derrière.
- **Axes.** **Draw the axes, with ticks and values**. *Placement* est
  automatique — sur les arêtes extérieures face à vous, qui suivent la vue
  quand elle tourne — ou choisi : pour chaque axe, l'arête du cube qu'il
  longe et s'il compte **from min** ou **from max**.
- **Camera** : azimut, élévation, rapport de profondeur ; **Look into it**
  fait face à un plan choisi. **Outline the faces** et un **Orientation
  tripod** sont facultatifs.
- **Save the view...** exporte l'image par la même fenêtre d'exportation que
  tout panneau (PNG, TIFF, JPEG à la taille et à la résolution choisies).
  **To figure composer...** ajoute la vue comme panneau d'une figure.

---

## 12. Arithmétique de coupes : dichroïsme et division par une référence

Deux coupes en donnent une troisième : dichroïsme linéaire (LH − LV),
dichroïsme circulaire, division par une référence d'or, ou simples
différences, rapports et sommes. Sélectionnez deux coupes dans la liste et
appuyez sur panneau → Advanced → **Cut arithmetic (two cuts)...** (la
première sélectionnée est A ; la fenêtre permet de les échanger). Le bouton
reste grisé tant qu'exactement deux coupes ne sont pas sélectionnées.

![Arithmétique de coupes](images/cut_arithmetic.png)

- **Le tableau du haut** compare ce que les deux mesures ont enregistré
  (polarisation, énergie de photon, angles, position de l'échantillon,
  température, énergie de passage, mode de lentille). Les lignes en orange
  diffèrent au-delà d'une tolérance — une carte de dichroïsme entre deux
  points de l'échantillon est une carte de *cette* différence-là.
- Préréglages **Purpose** : *Linear dichroism* (A − B, B mis à la même
  intensité totale), *Circular dichroism* ((A − B)/(A + B)), *Divide by a
  reference*, *Custom*. Un avertissement apparaît si les polarisations
  semblent inversées (**Swap A ↔ B**).
- **Scale B to A** par l'intensité totale ou celle d'une région ; **Hide
  where A+B below** masque le fond.
- A et B sont tracés sur une seule échelle de couleurs ; le résultat sur une
  échelle rouge-blanc-bleu symétrique autour de zéro. Avec des comptages
  bruts des deux côtés, l'incertitude de Poisson est propagée et résumée sous
  le graphique.
- **Result to list**, **Uncertainty to list**, **All three as a figure**.

---

## 13. Balayages en énergie de photon (kz)

Un balayage en énergie de photon est listé comme `kz map` ; son contour est
l'énergie de photon en fonction de l'angle, à une énergie donnée.

![Une carte kz](images/kz_map_contour.png)

### 13.1 Aligner les niveaux de Fermi d'une carte kz

Rien dans les fichiers n'indique où se trouve réellement le niveau de Fermi
de chaque spectre. Panneau → Basic → **Align Fermi levels (kz)...** le mesure
sur les spectres eux-mêmes ; la coupe selon la fente de la carte kz s'ouvre
avec lui, et les deux se ferment ensemble.

1. La coupe selon la fente s'ouvre avec un cadre. Faites glisser le cadre sur
   un front de Fermi — l'intervalle d'angle fixe les canaux sommés,
   l'intervalle d'énergie la fenêtre d'ajustement.

   ![Le cadre sur la coupe selon la fente](images/kz_processing_cut.png)

2. **Fit and assemble** : le même cadre (en indices du détecteur) est utilisé
   pour chaque spectre ; chacun est décalé pour placer son propre front à
   zéro ; le cube est rogné à l'intervalle couvert par tous ; et, par défaut,
   chaque spectre est **normalisé par son intensité totale**, ce qui élimine
   les variations du flux de photons.

![Aligner les niveaux de Fermi (kz map)](images/kz_processing_dialog.png)

E_F en fonction de l'énergie de photon est tracé. Une courbe régulière
signifie des ajustements fiables ; les points dont l'ajustement n'est pas
crédible sont signalés et remplacés par interpolation de leurs voisins.

### 13.2 Convertir une carte kz en impulsion

Panneau → Basic → **kz conversion...** convertit le balayage en (k_z, k∥). Il
part de la dernière conversion d'un balayage de même taille.

![kz → impulsion](images/kz_conversion_dialog.png)

- Le **potentiel interne V0** est l'inconnue : l'aperçu reconvertit une
  tranche en énergie à chaque modification, en quelques millisecondes.
- Avec le **réseau** (groupe d'espace, a, b, c), la **normale à la surface**
  et la direction dans le plan, les limites de zone sont tracées sur
  l'aperçu.
- **Scan V0** convertit pour une série de valeurs de V0 et indique où la
  période en k_z correspond au réseau — et avec quelle précision cela fixe V0.
- **Measure a k_z period...** : cliquez deux points correspondant à la même
  structure, distants d'une zone ; les plans du réseau réciproque compatibles
  sont listés.
- Le bandeau jaune avertit quand le front de Fermi est encore courbé en
  fonction de l'angle de l'analyseur — corrigez-le d'abord (Straighten bands (FS correction) sur
  la coupe selon la fente).
- **Convert the whole cube** ajoute une `kz map (k)` à la liste.

![La carte kz convertie](images/kz_map_in_k.png)

---

## 14. Supprimer la grille du détecteur

Les motifs du MCP (hexagonal) et de la grille (carré) s'impriment sur chaque
image. **À faire en premier** : seules les données encore sur les pixels du
détecteur peuvent être dégrillées ; le bouton refuse donc tout ce qui a déjà
été converti, corrigé ou lissé, et dit pourquoi.

Sélectionnez une carte (ou une carte kz) ou une coupe et appuyez sur panneau
→ Basic → **Remove detector grid...**. Pour une carte, aucune référence n'est
nécessaire — la grille reste sur les mêmes pixels dans chaque tranche, alors
que la photoémission se déplace. **Find and remove the grid** :

![Dégrillage](images/degrid_dialog.png)

- **Before / After** pour une tranche choisie, sur une même échelle de
  couleurs (cochez **Show −∂²I/∂E²** pour rendre toute grille évidente).
- **La grille trouvée**, ses **régions dans l'espace réciproque** (en rouge),
  et le **contraste et le décalage** ajustés pour chaque tranche.
- La puissance restant dans les régions de la grille (1 = plus de grille).
- **Also list the grid pattern** conserve un jeu de données `[grid]`, qu'une
  coupe prise plus tard avec le même mode de lentille et la même énergie de
  passage retrouvera et utilisera.

---

## 15. Courbes : EDC, MDC et données résolues en spin

**EDC → list** / **MDC → list** sous n'importe quelle image, ou un fichier
résolu en spin, donnent un jeu de données de courbes, qui s'ouvre dans la
**fenêtre de courbes** :

![Fenêtre de courbes](images/curve_viewer.png)

Les canaux peuvent être masqués, décalés en cascade (**Offset**), affichés
avec des **Error bars** (à partir d'un canal σ), et **Range** délimite une
fenêtre en x utilisée par les opérations. *Functions → Transform*
propose **Crop**, **Bin**, **Normalise**, **Subtract a background**
(constant, linéaire, Shirley), **Shift the axis** et **Add counting errors
(√N)** ; chacune crée un nouveau jeu de données et propage correctement σ.

### 15.1 Ajustement de courbes

*Functions → Curve fit...* : **Peaks on a background** (pics lorentziens,
gaussiens ou de Voigt, placés par clic ou par **Find peaks** ; fond nul /
constant / linéaire / Shirley ; résolution instrumentale ; **× Fermi edge**
facultatif), ou le **Fermi edge** lui-même. La pondération est le σ propre au
canal, √N, ou uniforme. Chaque paramètre est donné avec son erreur ; le
résidu est affiché en unités de σ.

![Ajustement de courbe](images/curve_fit_dialog.png)

**Fit to list** conserve les données, l'ajustement, le résidu et les
composantes ; **Copy results** copie un tableau séparé par des tabulations.
Un ajustement du front de Fermi permet **List a copy with E_F = 0**.

### 15.2 Analyse de spin

Une EDC de spin s'ouvre avec une courbe par canal de spin :

![EDC de spin](images/spin_edc_viewer.png)

*Functions → Spin analysis...* :

![Analyse de spin](images/spin_analysis.png)

- Le tableau lit l'étiquette de chaque canal (réglage du manipulateur,
  aimantation de la cible, axe et signe comptés comme « haut ») ; **Axis** et
  **Counts as** peuvent être modifiés.
- Le **Cross ratio** (par défaut quand les canaux le permettent) élimine
  exactement à la fois la transmission et la réflectivité de la cible ; une
  **seule paire** est disponible mais porte l'asymétrie instrumentale, qui est
  mesurée et indiquée (ε).
- **S_eff** (fonction de Sherman, 0,2 par défaut), **Bin**, fenêtre **Zero
  reference** facultative.
- Graphiques : canaux bruts, P avec ses barres d'erreur, et I↑ / I↓ avec
  leurs erreurs.
- **P to list**, **I↑ / I↓ to list**, **As a figure** (ouvre le composeur de figures).

---

## 16. Balayages résolus spatialement (nano-ARPES)

Un balayage en espace réel ouvre la **carte spatiale** (à gauche, l'intensité
totale en chaque position) avec le **spectre E-k au curseur** (à droite).
Faites glisser la cible rouge pour lire le spectre en un autre point. Un clic
droit sur un panneau donne les deux intégrations :

- sur la carte, **Integrate selection → E vs k panel** somme les spectres sur
  une région de l'échantillon ;
- sur le spectre, **Integrate selection → spatial panel** cartographie sur
  l'échantillon l'intensité d'une fenêtre (k, E) choisie — c'est ainsi qu'on
  repère un domaine ou le bord d'un flocon.

![Balayage spatial](images/spatial_scan.png)

---

## 17. Figures et exportation

### 17.1 Composeur de figures

Panneau → Advanced → **Figure composer...** **sans rien sélectionner** ouvre
une page vierge, comme une figure d'article : quatre panneaux vides (a),
(b), (c), (d) à la largeur de colonne d'une revue. Avec des coupes (ou des
courbes) sélectionnées, chacune devient un panneau (de même pour le bouton
**Figure composer...** d'une tranche figée ; un empilement de courbes ou une
vue 3D a son propre bouton **To figure composer**).

![Composeur de figures](images/figure_composer.png)

- **Mise en page.** Cliquez un panneau pour le sélectionner ; faites-le
  glisser pour le déplacer, tirez ses poignées pour le redimensionner (Maj
  sur un coin garde ses proportions). Les panneaux s'alignent sur les bords
  des autres et sur une grille de 0,5 mm (**Snap**) ; les flèches du clavier
  le décalent de 0,5 mm (Maj : 5 mm). *Arrange as grid* (onglet Panels à
  droite) les remet en lignes × colonnes. **Re-letter** les nomme (a), (b)…
  dans l'ordre de lecture.
- **Remplir un panneau.** Double-cliquez un panneau vide, ou **Load from
  list...** / **Load picture...** à gauche : un jeu de la liste principale,
  ou une image d'un fichier (PNG, JPEG, TIFF… — un schéma ou une figure
  faite ailleurs ; elle garde ses proportions et n'a pas d'axes).
- **Encarts.** **Zoom-in inset** : tracez une région d'un panneau — un
  petit panneau qui l'agrandit apparaît dans le coin, et la région est
  encadrée. Un **Panel** tracé dans un autre panneau est un encart à remplir
  avec un autre jeu. Les encarts n'ont pas de lettre et suivent leur panneau.
- **Dessiner.** **Arrow**, **Line**, **Box**, **Ellipse**, **Text** :
  faites glisser (ou cliquez, pour le texte) sur la page. Un dessin
  sélectionné se déplace, ses extrémités ou coins se tirent, et sa couleur,
  son épaisseur, son trait, son remplissage et son texte se règlent au-dessus
  des outils. Double-cliquez un texte pour le modifier ; Suppr efface la
  sélection ; clic droit pour le reste.
- **Axes et polices.** Les outils à côté de la page agissent sur le panneau
  sélectionné (ou sur tous) : bornes, graduations, étiquettes, bords
  numérotés, taille de police propre au panneau, carte de couleurs, barre de
  couleurs, lignes de référence.
- **Undo / Redo** (Ctrl+Z / Ctrl+Y) couvrent tout. **PDF**, **PNG**, **SVG**
  exportent, **Copy** copie dans le presse-papiers ; les panneaux vides ne
  sont jamais écrits. Les **Figure tools** du composeur ne sont pas les
  **Plot tools** d'une fenêtre, qui dessinent des guides sur le graphique
  affiché.

### 17.2 Figure en série de tranches

Sélectionnez une carte → panneau → Advanced → **Slice series figure...**
dispose une page de tranches 2D planes du cube en une figure pour un
article — contours à énergie constante à une série d'énergies, ou coupes à
une série d'angles ou d'impulsions — régulièrement espacées ou à des valeurs
choisies, avec intégration, une échelle de couleurs commune et des
étiquettes seulement sur les bords extérieurs. (À ne pas confondre avec
**3D slices**, qui trace trois plans dans une vue 3D.)

| Réglages | Résultat |
|---|---|
| ![Réglages de la série de tranches](images/slice_figure_dialog.png) | ![Figure en série de tranches](images/slice_figure.png) |

### 17.3 Exporter un panneau

Clic droit sur n'importe quelle image → **Export** → **Data (.nxs)...**,
**Image...** ou **Axes only (.eps)...**. L'exportation reprend toujours
exactement ce que montre le panneau (étendue comprise) : les données, comme
un jeu de données que ce programme sait relire ; une image avec la palette
appliquée (PNG, TIFF ou JPEG, à la taille et à la résolution choisies, avec
ou sans les axes, avec les guides des Plot tools s'il y en a, et avec la zone
de Brillouin si le jeu de données en a une) ; ou le cadre, les graduations et
les étiquettes en PostScript vectoriel, pour y placer l'image dans
Illustrator (avec, au choix, la zone de Brillouin tracée dans le cadre — voir
[§8.2](#82-zone-de-brillouin)).

![Exportation](images/export_dialog.png)

**Save datasets...** (panneau d'outils, cadre Datasets) écrit les jeux de
données sélectionnés dans un fichier `.nxs` de votre choix ; plusieurs
lignes sélectionnées vont dans un même fichier.

---

## 18. Sauvegarde automatique, mémoire et récupération

- **Sauvegarde automatique.** Chaque jeu de données calculé (une carte k, une
  coupe traitée, un résultat d'ajustement…) est écrit, sans compression, dans
  le dossier de la session sous `~/.arpes_viewer/sessions/` dès qu'il existe.
  Les balayages en espace réel ne le sont pas (ce sont les mesures
  elles-mêmes, et elles sont volumineuses).
- **Une fermeture normale supprime les copies automatiques.** Quand vous
  fermez la fenêtre principale et confirmez (tout ce que vous voulez garder
  est enregistré, ou *Close anyway*), les fenêtres ouvertes sont fermées et le
  dossier de la session est supprimé immédiatement. La copie d'un jeu de
  données disparaît aussi dès que vous l'enregistrez dans un fichier à vous.
- **Récupération après un plantage.** Ce n'est que si le programme ne s'est
  pas fermé normalement (plantage, arrêt forcé) que le dossier de la session
  est conservé : au démarrage suivant, ses jeux de données vous sont
  proposés, et il est supprimé au bout de trois jours.
- **Journal des opérations.** `~/.arpes_viewer/sessions/operations.log`, une
  ligne par écriture / enregistrement / suppression — pour retrouver ce qui
  était où.
- **Cadre Memory** (en bas du panneau d'outils). La mémoire du programme et celle du système, le nombre de
  fenêtres et de fichiers HDF5 ouverts, et la quantité de données calculées
  conservées. **Free memory** abandonne les copies en cache des données
  calculées (elles sont relues depuis leur fichier de sauvegarde au besoin) ;
  **Close all viewers** ferme toutes les fenêtres puis libère la mémoire.
- **Open auto-save folder** (dans le cadre Memory) ouvre le dossier des
  sessions dans l'explorateur de fichiers, pour supprimer à la main les
  anciennes copies de travail. Le dossier de la session en cours (nommé dans
  la barre d'état) est utilisé : n'y touchez pas.
- **Les opérations longues** (lecture d'une grande carte, conversion,
  dégrillage, assemblage d'un dossier) tournent hors du fil de l'interface,
  avec une barre de progression et un bouton Cancel qui fonctionne.
- **La fermeture** avec des données calculées non enregistrées demande
  d'abord confirmation, en nommant les jeux de données : **Save to a
  file...**, **Close anyway** (les copies automatiques sont supprimées) ou
  **Cancel** (rien ne se passe).

## 19. Conseils et dépannage

- **Un bouton d'outil grisé** — survolez-le : l'infobulle dit pourquoi (par
  ex. l'ajustement MDC exige une coupe déjà en Å⁻¹, l'arithmétique de coupes
  deux coupes). Sélectionnez les bonnes lignes dans la liste et il s'allume.
- **Un outil manque au menu Functions d'une fenêtre** — il est dans le
  panneau d'outils (voir *Où se trouvent les choses maintenant* au
  [§3](#3-la-fenêtre-principale)) ; regardez les deux onglets.
- **Une fenêtre est plus haute que l'écran** — elle s'ouvre au plus à la
  hauteur de l'écran et son contenu défile ; les boutons Apply / Close
  restent en bas. Toutes les fenêtres se redimensionnent.
- **Aucune conversion en k n'est proposée** — le premier axe de la carte
  n'est pas un angle (énergie de photon, température…). Vérifiez *First axis
  of a map is* dans la fenêtre de chargement.
- **Conversion en k d'une coupe en E − E_F** — le décalage en énergie est
  rempli pour vous à partir de l'ajustement de Fermi ; si vous avez décalé
  l'axe autrement, saisissez dans *Energy offset* l'énergie cinétique de E_F
  (énergie de photon − travail de sortie).
- **Un fichier est mal lu** — choisissez explicitement le lecteur dans la
  fenêtre de chargement, et voyez *Assumptions worth knowing about* dans le
  README (ordre des axes dans les fichiers SOLEIL).
- **Le curseur est lent sur une grande carte** — une carte lue à la demande
  dans son fichier se recoupe en environ 0,1 s par fenêtre, si bien que les
  autres fenêtres attendent que la souris s'arrête. Une fois la carte lue en
  mémoire, elle suit en direct.
- **Saisie des nombres** — toutes les cases numériques prennent un point
  décimal (`0.05`), quelle que soit la langue du système ; la virgule n'est
  pas acceptée.
- **La machine utilise la mémoire d'échange** — surveillez la barre Memory
  (ambre au-delà d'un tiers de la RAM, rouge au-delà de la moitié) ; utilisez
  **Free memory** ou **Close all viewers**.
