# Holographie des parcelles d’observateur

**La physique à partir des registres que les observateurs peuvent partager.**

[Article principal](https://philpapers.org/rec/MUEFOC) · [Livre](https://oph-book.floatingpragma.io/) · [Carnet](https://notebook.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a) · [Bibliothèque Lean](Lean/) · [Reproduire](REPRODUCE.md) · [Pragma Research](https://floatingpragma.io/) · [English](README.md)

[![Mandatory Suite](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml)
[![Lean CI](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml)
[![Runnable observer model](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml)

L’Holographie des parcelles d’observateur (OPH) propose que la physique est
ce que des observateurs bornés peuvent établir ensemble. Des observateurs
finis comparent leurs registres sur des frontières partagées, réparent leurs
désaccords et tiennent pour réel ce qui résiste à toutes les comparaisons.
Trois problèmes qui résistent à une description classique motivent le
programme : la mesure quantique, les corrélations qu’aucun modèle classique
local ne reproduit et l’origine d’un espace-temps apparent. OPH cherche une
explication commune des trois dans la cohérence des registres et des
opérations des observateurs, et la poursuit sous la forme d’un observateur
exécutable, de mathématiques vérifiées par machine et de prédictions fixées à
l’avance.

**L’Univers se lit lui-même jusqu’à l’existence.** OPH a une seule
hypothèse, l’hypothèse de clôture : l’Univers est une structure intemporelle
et autoréférentielle dont les observateurs retrouvent la spécification de
leur propre monde à partir de leurs registres et construisent un matériel qui
réalise cette même histoire, y compris ses constructeurs. Rien ne se trouve
hors de cette boucle ; la spécification retrouvée n’a donc besoin d’aucune
autre justification. Elle doit fonctionner et reproduire ce que mesurent les
observateurs. L’espace, le temps et la matière sont des descriptions
effectives des registres, et « notre Univers » est le registre public, ce qui
résiste à toute comparaison entre observateurs. La [base](docs/BASIS.md)
énonce l’hypothèse, la règle de l’absence d’extérieur, la condition d’accord
et la spécification exécutable ; le
[guide des termes](docs/INTERPRETATION_AND_TERMS.md) explique la vision
intemporelle et les horloges qu’elle contient.

## Points forts

- **Un observateur qui fonctionne.** Le plus petit observateur a douze ports
  et deux événements. Il écrit un registre à chaque événement, le conserve et
  ne laisse que ce registre choisir la réparation suivante, et le registre
  change ce qui est mesuré : l’occupation finale passe de 1/6 à 35/264. Un
  seul processus quantique fini réunit le flux électrique, la matière chargée
  et un champ scalaire sur un même registre, et un calcul exact certifie que
  sa rétroaction enregistrée modifie la réponse mesurée. Son protocole
  continu à mesure lisse exécute quatre lectures du même état perturbé,
  compile deux lois d’enregistrement quantiques de la réparation native et protège un qubit codé
  par une récupération commandée par les registres, sous des services déclarés.
  Deux observateurs dont les régions se recouvrent lisent leurs champs et
  leurs horloges dans une même histoire continue et s’accordent sur une même
  loi complète. Des supports joints forment des écrans allant jusqu’à
  245 760 ports, et sur un écran qui se raffine une seule histoire de
  réparation définit un champ continu aux fluctuations gaussiennes, dont les
  ondes se propagent à vitesse finie.
- **La théorie relativiste des champs comme description effective des
  observateurs.** Sur un produit radial déclaré au-dessus de cet écran,
  l’action d’onde native a une densité de signature lorentzienne, et pour tout
  ensemble fini de lectures d’horloge, de champ et de réponse des
  observateurs, une seule théorie des champs lorentzienne lisse et statique en
  dimension 1+3 les reproduit toutes, à toute tolérance choisie. Le
  [bilan intégré](model/CONSTRUCTION.md#integrated-effective-description)
  en énonce les limites : la métrique native exacte n’est pas lisse dans son
  atlas actuel, et la structure radiale est fournie. Les registres de courant
  du support dérivent pourtant un rayon : leur compteur tridimensionnel suit
  la loi de diffusion radiale standard depuis tout état initial. Les
  histoires de mesure complètes ne distinguent pas les structures radiales,
  si bien que relier ce rayon à la direction radiale du champ et à la
  distance physique reste ouvert.
- **La physique standard à partir d’un accord fini.** Des théorèmes finis
  mènent des lectures de l’observateur à des structures familières : des
  registres publics indépendants de l’ordre des réparations ; la règle de
  Born, qui fixe les probabilités quantiques, et la borne de Tsirelson sur les
  registres achevés ; un support tridimensionnel dans sa réponse lente ; l’algèbre de
  Lie de jauge du Modèle standard avec une génération sans anomalies ; un
  ordre causal de dimension (1+3) dont le temps propre vient de comptages de
  registres ; et l’équation d’Einstein issue d’un bilan nul réalisé.
- **Des mathématiques vérifiables.** La [bibliothèque Lean](Lean/) contient
  plus de 12600 théorèmes et lemmes publics sans preuve admise, et chaque reçu
  de calcul conservé est rejoué en intégration continue.
- **Des tests fixés à l’avance.** Les prédictions et leurs règles de rejet
  entrent dans l’[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md)
  avant toute donnée de comparaison, avec des horodatages attestés sur la
  chaîne de blocs Bitcoin ; la prédiction de dispersion pour l’orbite des ports
  primitifs en fait partie. Le [registre des postdictions](docs/POSTDICTION_LEDGER.md)
  et le [programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md)
  séparent les comparaisons avec des mesures des tests des mathématiques.

**L’accord avec l’ensemble des statistiques de mesure observées n’est pas
établi.** Les théorèmes valent sous des prémisses énoncées, et l’observateur
exécutable est un candidat jugé selon des critères d’acceptation déclarés.
Les sections ci-dessous indiquent quelles prémisses la construction fournit
et lesquelles restent des entrées.

## Lancer l’observateur

Depuis la racine du dépôt, avec Python 3.11 ou plus récent :

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

L’exécution affiche une ligne : le risque de Brier de trois prévisions du
registre suivant et l’empreinte du résultat.

```text
{"brier": {"action_known_prior": 0.0728, "native_record": 0.0563, "within_action_shuffle": 0.0884}, "reference_matches": 10000, ...}
```

Plus la valeur est basse, mieux c’est. La prévision qui utilise le registre
conservé bat celle qui ne connaît que l’action choisie et celle dont les
registres sont mélangés au sein d’une même action ; les valeurs exactes sont
395/7392, 140/2013 et 12875/150304. Le [répertoire du modèle](model/README.md)
contient l’implémentation, ses tests et ses résultats conservés. Ses
expériences conjointes comprennent l’[observateur cohérent sur écran](model/coherent_screen_observer/README.md),
qui fait traverser à un qubit logique et à une référence intacte des tics
d’horloge natifs et une rétroaction commandée par les registres sur des
écrans de 15 360 et 245 760 ports, et le
[processus commun en interaction](model/common_source/README.md), dans lequel
le flux électrique quantique, des fermions chargés et un champ scalaire
partagent un même registre fini avec des instruments complets, des registres
conservés et une rétroaction native ; l’action du champ et la préparation de
ce processus sont des entrées fournies. Le
[guide des lois de réparation](model/CANDIDATES.md) montre quelles parties de
la loi de réparation sont imposées et lesquelles restent libres : la règle
classique découle de quatre exigences, ses lois d’enregistrement quantiques
forment une seule famille à cinq multiplicateurs de cohérence, une seule
expérience d’interférence lit la cohérence qui sépare la loi exécutée C de la
loi de comparaison E, et neuf mesures compilées identifient la table entière. Le
[paquet de prédictions](model/predictions/README.md) fixe cette expérience et
sa règle de rejet par comptage avant tout comptage indépendant, et le
[guide de construction](model/CONSTRUCTION.md) énonce les critères
d’acceptation et les connexions qui restent à construire.

Un seul support à douze ports relie les observations spatiales et
d’horloge : les registres de courant retrouvent la géométrie locale, et des
lois exactes de corrélation temporelle reconstruisent l’espace de Hilbert
cyclique d’un détecteur et ses modes de transfert visibles. La réparation fixe
une action de fluctuation du courant et une loi de probabilité reliant les
passages incidents aux tics de mouvement. Des exécutions conservées testent
cette conversion contre une horloge hôte. Sur des écrans sources plus grands,
les prédictions de champ et de rétroaction valent aussi à occupation
conservée fixe, sans occupations de Bernoulli préparées indépendamment. Deux
observateurs qui partagent un port d’un même support lisent leurs champs,
leurs horloges locales et leur rétroaction dans une même histoire continue,
et l’exécution native directe comme une description effective donnent la
même probabilité à chaque registre conjoint ; une version construite à partir
d’impulsions sources compilées reproduit cette histoire sur les 512 adresses
de code. Sur les écrans sphériques qui se raffinent, une seule histoire de
réparation donne un [opérateur de champ continu](model/physical_space_time/OCCUPATION_FIELD.md)
dont les champs gaussiens, la décroissance et les ondes vivent sur l’écran
bidimensionnel lui-même. Le [bilan de récupération](model/OBSERVER_RELATIVE_RECOVERY.md)
sépare ces résultats finis d’un espace-temps relativiste étendu et d’une
identification de laboratoire calibrée. Une
[relecture de données publiques de Bell](model/predictions/BELL_RECORDS.md)
teste la limite de l’interprétation de la branche diagonale comme modèle
complet de réponse locale ; elle ne teste pas la construction quantique
complète. L’[audit des données publiques](docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md)
compare l’inclinaison scalaire conditionnelle aux chaînes postérieures d’ACT DR6
et à une tranche de vraisemblance ACT, et teste un rattachement direct des comptages d’horloge natifs aux rapports
publiés d’horloges optiques.

Le [bilan empirique](docs/EMPIRICAL_STATUS.md) distingue les comparaisons
avec les mesures de la sélection de la règle de réparation native. La réparation
locale de l’histoire reproduit les lectures de l’appareil oscillateur fourni ;
les données conservées ne déterminent pas la loi d’enregistrement native.

## Fondements mathématiques

- **L’accord et les registres stables.** Sous des conditions explicites de
  terminaison et de cohérence, les réparations locales atteignent un même
  résultat public quel que soit leur ordre. D’autres résultats contrôlent
  quels registres survivent au raffinement.
  [Article sur le consensus](paper/reality_as_consensus_protocol.pdf).
- **Les observations quantiques.** Dans l’algèbre quantique finie énoncée, des
  résultats vérifiés couvrent les mesures répétables, les registres classiques
  que les opérations ultérieures ne peuvent perturber, l’incertitude et la
  borne de corrélation de Bell, avec un état qui l’atteint.
  [Algèbres finies d’événements quantiques](extra/machine_checked_finite_event_algebras.pdf).
- **Trois directions spatiales.** Le support fixe à douze ports et sa règle de
  réparation moyenne donnent une réponse lente exactement tridimensionnelle de
  géométrie euclidienne, que retrouvent les registres de courant d’un observateur.
  [Construction du support](code/a5_closure/README.md).
- **La symétrie de jauge.** Sous des prémisses énoncées de réponse réversible
  complète et de transport entre observateurs, la construction impose
  l’algèbre de Lie de jauge du Modèle standard, et une représentation de
  matière spécifiée donne une génération sans anomalies.
  [Dérivation de jauge](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **La thermodynamique.** Une seule construction finie de réparation donne les
  quatre lois, les relations de fluctuation et la borne d’effacement de
  Landauer sous des prémisses énoncées de référence, de calibration de
  l’énergie et de basse température.
  [Exposé thermodynamique](paper/observers_are_all_you_need.pdf).
- **La capacité d’information.** Des résultats exacts comptent les registres
  qui restent distinguables après les opérations permises, bornent ce nombre
  par le support disponible et le suivent lors du raffinement.
  [Capacité des registres](code/capacity_readback/README.md).
- **La matière et les champs.** Une action fournie de scalaire chargé et de
  Maxwell, avec sa géométrie et sa quantification, donne un hamiltonien en
  interaction après réduction de jauge, et un raffinement scalaire contrôlé
  relie les lectures classiques d’un détecteur à une équation de champ non
  linéaire. [Continuum de matière](paper/observers_are_all_you_need.pdf).

Ce sont des résultats mathématiques conditionnels. Faire tenir leurs prémisses
pour les observations d’un même modèle commun est le problème central de la
construction. La [carte des prémisses](model/PREMISES.md) indique quelles
prémisses la spécification exécutable fournit et lesquelles elle prend comme
entrées.

## Articles

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  est la synthèse principale : la théorie en un seul exposé, des résultats
  choisis en détail, le tableau de bord de la construction et l’hypothèse de
  clôture. Citez sa version publiée, la
  [fiche PhilPapers](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer for Observer Patch Holography](paper/executable_observer_construction.pdf)**
  énonce la spécification, ses critères d’acceptation, ses composants
  certifiés et la carte des prémisses.
- **La bibliothèque de théorèmes** contient les preuves longues, du
  [consensus](paper/reality_as_consensus_protocol.pdf) et des
  [observateurs et registres](paper/observers_are_all_you_need.pdf) jusqu’à
  l’[espace-temps et la dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).
  L’[index des articles](paper/README.md) et l’[ensemble des articles](docs/PAPER_SET.md)
  donnent chaque article et son rôle.
- **Les articles compagnons** portent des résultats ciblés. Une relation entre
  masses de leptons, la relation de Koide, tient exactement sous une prémisse
  d’équilibre déclarée. Le programme de capacité demande si la capacité
  publique attribuée à l’Univers s’accorde avec celle reconstruite depuis
  l’intérieur. Ces calculs de constantes, comme celui de la structure fine,
  ont un statut diagnostique.

## Où se trouvent les choses

| Chemin | Contenu |
| --- | --- |
| [model/](model/README.md) | L’observateur exécutable : implémentation, tests, reçus conservés et enregistrements de protocole |
| [flagship/](flagship/) | L’article de synthèse |
| [paper/](paper/README.md) | L’article de construction et la bibliothèque de théorèmes, avec les articles archivés sous [paper/archive/](paper/archive/README.md) |
| [extra/](extra/README.md), [cosmology/](cosmology/) | Les articles compagnons |
| [Lean/](Lean/) | La bibliothèque vérifiée par machine |
| [code/](code/README.md) | Les calculs exacts, les vérificateurs et les reçus derrière les articles |
| [evidence/](evidence/README.md) | Les paquets de preuves conservés et les enregistrements de garde |
| [claims/](claims/README.md) | Les énoncés établis et leurs dépendances |
| [tracking/](tracking/) | Les registres des prémisses, observations, prédictions et instruments |
| [docs/](docs/README.md) | La base, le guide des termes, les politiques et les registres |
| [book/](book/), [essays/](essays/) | Le livre et les essais |

## Reproduire et vérifier

Le [guide de reproduction](REPRODUCE.md) lance la suite de vérification depuis
un clone propre. Le vérificateur du modèle rejoue chaque reçu conservé :

```sh
python -m pip install -r model/requirements.txt
python -m model.verify
```

La [bibliothèque Lean](Lean/) se compile avec Lake, l’[index de la documentation](docs/README.md)
donne les registres et les politiques, et le [livre](book/) offre un chemin
d’accès plus long. Les manuels sur [learn.floatingpragma.io](https://learn.floatingpragma.io/)
et la simulation interactive sur
[simulation.floatingpragma.io](https://simulation.floatingpragma.io/) sont
les introductions les plus accessibles.

## Contribuer

La tâche consiste à faire qu’un seul processus d’observation et de réparation
fournisse la physique effective : relier l’horloge de l’observateur, son
interprétation spatiale et ses mesures de champ sur les mêmes histoires, puis
comparer les probabilités de résultats, les corrélations et les réponses aux
interventions avec des observations indépendantes. Une équation de champ
fournie reste une entrée même lorsque sa simulation fonctionne. Les
contributions doivent résoudre une connexion nommée ou révéler un défaut du
modèle commun ; le [guide de construction](model/CONSTRUCTION.md), le
[guide des lois de réparation](model/CANDIDATES.md) et les
[tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
décrivent ce travail. Les publications d’articles suivent la
[politique de publication](model/CONSTRUCTION.md#research-programme-and-paper-release-hold) ;
les PDF de la branche principale sont des aperçus des sources, et les versions
déposées conservent leurs octets d’origine.

## Licence

Le logiciel, y compris Lean et le code du modèle, est publié sous
[Apache-2.0](code/LICENSE). Les articles, la documentation, les figures, les
données et les registres sont publiés sous [CC BY-NC-SA 4.0](LICENSE), avec les
conditions amont du Particle Data Group pour les valeurs de particules
empaquetées. La [carte des licences](LICENSE) donne les conditions par
répertoire.

<details>
<summary>Projets associés et financement</summary>

[Cadence](https://github.com/muellerberndt/cadence) est le projet associé de
machines apprenantes. Ses éléments de preuve sur les capacités sont distincts
des éléments de preuve physiques d’OPH.
[Pragma Research](https://floatingpragma.io/) réunit les projets.

Le memecoin officiel TEO (Theory Of Everything) d’OPH a été déployé par la
communauté Base et se trouve sur Base à l’adresse
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([graphique](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Ses frais de transaction sont versés à Bernhard Mueller pour financer le
projet. Le jeton n’a aucun rôle dans la théorie. Les autres jetons qui
utilisent le nom OPH ou TEO n’ont aucun lien avec le projet. Rien ici ne
constitue un conseil financier.

</details>

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->
