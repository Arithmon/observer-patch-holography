# Holographie des parcelles d’observateur

**La physique à partir des registres que les observateurs peuvent partager.**

[Article principal](https://philpapers.org/rec/MUEFOC) · [Observateur exécutable](model/README.md) · [Bibliothèque Lean](Lean/) · [Pragma Research](https://floatingpragma.io/research/) · [English](README.md)

L’Holographie des parcelles d’observateur (OPH) propose que la physique est
ce que des observateurs bornés peuvent établir ensemble. Des observateurs
finis comparent leurs registres sur des frontières partagées, réparent leurs
désaccords et tiennent pour réel ce qui résiste à toutes les comparaisons. Le
projet construit cette image sous la forme d’un observateur exécutable,
démontre ce qui en découle par des mathématiques vérifiées par machine, et
confronte les résultats à des prédictions fixées à l’avance.

**L’Univers se lit lui-même jusqu’à l’existence.** OPH a une seule
hypothèse : l’Univers est une structure intemporelle et autoréférentielle
dont les observateurs retrouvent la spécification de leur propre monde à
partir de leurs registres et construisent un matériel qui réalise cette même
histoire, y compris ses constructeurs. Rien ne se trouve hors de cette
boucle ; la spécification retrouvée n’a donc besoin d’aucune autre
justification. Elle doit fonctionner et reproduire ce que mesurent les
observateurs. L’espace, le temps et la matière sont des descriptions
effectives des registres, et « notre Univers » est le registre public, ce qui
résiste à toute comparaison entre observateurs. La [base](docs/BASIS.md)
énonce l’hypothèse, sa règle et la spécification exécutable ; le
[guide des termes](docs/INTERPRETATION_AND_TERMS.md) explique la vision
intemporelle et les horloges qu’elle contient.

## Points forts

- **Un observateur qui fonctionne.** Un support à douze ports muni d’une loi
  de réparation native écrit un registre à chaque événement, le conserve et
  ne laisse que ce registre commander la réparation suivante. Un seul
  processus quantique fini réunit le flux électrique, la matière chargée et
  un champ scalaire sur un même registre, et un calcul exact certifie que sa
  rétroaction enregistrée modifie la réponse mesurée.
- **La physique standard à partir d’un accord fini.** Des théorèmes finis
  mènent des lectures de l’observateur à des structures familières : des
  registres publics indépendants de l’ordre des réparations, la règle de Born
  et la borne de Tsirelson sur les registres achevés, un support
  intrinsèquement tridimensionnel, l’algèbre de Lie de jauge du Modèle
  standard avec une génération sans anomalies, un précurseur causal de
  dimension (1+3) dont le temps propre vient de comptages de registres, et
  l’équation d’Einstein issue d’un bilan nul réalisé.
- **Des mathématiques vérifiables.** La [bibliothèque Lean](Lean/) contient
  plus de 12300 théorèmes et lemmes publics sans preuve admise, et les
  résultats calculés se rejouent à partir de reçus conservés.
- **Des tests fixés à l’avance.** Les prédictions et leurs règles de rejet
  sont enregistrées avant les données de comparaison. Le
  [registre des postdictions](docs/POSTDICTION_LEDGER.md) et le
  [programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md) séparent
  les comparaisons avec des mesures des tests des mathématiques.

**L’accord avec l’ensemble des statistiques de mesure observées n’est pas
établi.** Les théorèmes valent sous des prémisses énoncées, et l’observateur
exécutable est un candidat jugé selon des critères d’acceptation déclarés.
Les sections ci-dessous indiquent quelles prémisses la construction fournit
et lesquelles restent des entrées.

## L’observateur exécutable

Le [répertoire du modèle](model/README.md) contient l’implémentation, ses
tests et ses résultats conservés. L’observateur lui-même est un programme
court : douze ports, une loi de réparation qui écrit une pièce et un courant
à chaque événement, des registres conservés et une rétroaction que seuls ces
registres commandent. Lancez-le depuis la racine du dépôt avec Python 3.11 ou
plus récent :

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

Son expérience conjointe la plus forte est le
[processus commun en interaction](model/common_source/README.md) : le flux
électrique quantique, des fermions chargés et un champ scalaire partagent un
même registre fini, avec des instruments complets, des registres conservés et
une rétroaction native. Une borne exacte certifie que la rétroaction modifie
la réponse mesurée, y compris dans la limite d’occupation infinie, et les
mêmes histoires fournissent des espaces de Hilbert d’observateur, des
horloges modulaires et des lectures géométriques dans trois directions.
L’action du champ et la préparation sont des entrées fournies. Le
[guide de construction](model/CONSTRUCTION.md) énonce les critères
d’acceptation et les connexions qui restent à construire, et le
[paquet de prédictions](model/predictions/README.md) fixe deux prédictions
concurrentes de détecteurs et leur test de rejet avant tout comptage
indépendant.

## Fondements mathématiques

- **L’accord et les registres stables.** Sous des conditions explicites de
  terminaison et de cohérence, les réparations locales atteignent un même
  résultat public quel que soit leur ordre. D’autres résultats contrôlent
  quels registres survivent au raffinement.
  [Article sur le consensus](paper/reality_as_consensus_protocol.pdf).
- **Les observations quantiques.** Dans l’algèbre quantique finie énoncée, des
  résultats vérifiés couvrent les mesures répétables, les registres classiques
  protégés, l’incertitude et la borne de corrélation de Bell, avec un état qui
  l’atteint. [Algèbres finies d’événements quantiques](extra/machine_checked_finite_event_algebras.pdf).
- **Trois directions spatiales.** Le support fixe à douze ports et sa règle de
  réparation moyenne donnent une réponse lente exactement tridimensionnelle de
  géométrie euclidienne, la brique géométrique de l’espace d’un observateur.
  [Construction du support](code/a5_closure/README.md).
- **La symétrie de jauge.** Avec une réponse réversible complète et des
  prémisses de transport entre observateurs, la construction impose l’algèbre
  de Lie de jauge du Modèle standard, et une représentation de matière
  spécifiée donne une génération sans anomalies.
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
  clôture. Sa version publiée se trouve sur
  [PhilPapers](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer](paper/executable_observer_construction.pdf)**
  énonce la spécification, ses critères d’acceptation, ses composants
  certifiés et la carte des prémisses.
- **La bibliothèque de théorèmes** contient les preuves longues, du
  [consensus](paper/reality_as_consensus_protocol.pdf) et des
  [observateurs et registres](paper/observers_are_all_you_need.pdf) jusqu’à
  l’[espace-temps et la dynamique d’Einstein](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).
  L’[index des articles](paper/README.md) et l’[ensemble des articles](docs/PAPER_SET.md)
  donnent chaque article et son rôle.
- **Les articles compagnons** portent des résultats ciblés : la relation de
  Koide tient exactement sous une prémisse d’équilibre déclarée, et le
  programme de capacité demande si la capacité publique attribuée à l’Univers
  s’accorde avec celle reconstruite depuis l’intérieur. Ces calculs de
  constantes, comme celui de la structure fine, ont un statut diagnostique.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Explorer et reproduire

Lisez l’[article principal](https://philpapers.org/rec/MUEFOC) pour l’exposé
complet, le [livre](book/) pour un chemin d’accès plus long et
l’[index de la documentation](docs/README.md) pour les registres et les
politiques. Le [guide de reproduction](REPRODUCE.md) rejoue chaque reçu
conservé et lance la suite de vérification depuis un clone propre.

<details>
<summary>Travaux associés, financement et licences</summary>

[Cadence](https://github.com/muellerberndt/cadence) est le projet associé de
machines apprenantes. Ses éléments de preuve sur les capacités sont distincts
des éléments de preuve physiques d’OPH.
[Pragma Research](https://floatingpragma.io/) réunit les projets.

Le memecoin officiel TEO (Theory Of Everything) d’OPH est sur Base à l’adresse
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([graphique](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Bernhard Mueller peut réclamer ses frais de transaction, qui aident à financer
la recherche. Le jeton n’a aucun rôle dans la théorie. Les autres jetons qui
utilisent le nom OPH ou TEO n’ont aucun lien avec le projet. Rien ici ne
constitue un conseil financier.

Le logiciel, y compris Lean et le code du modèle, est publié sous
[Apache-2.0](code/LICENSE). Les articles, la documentation, les figures, les
données et les registres sont publiés sous [CC BY-NC-SA 4.0](LICENSE), avec les
conditions amont du Particle Data Group pour les valeurs de particules
empaquetées. La [carte des licences](LICENSE) donne les conditions par
répertoire.

</details>

## Axe de recherche

**Faire qu’un seul processus d’observation et de réparation fournisse la
physique effective.** La tâche consiste à relier l’horloge de l’observateur, son
interprétation spatiale et ses mesures de champ sur les mêmes histoires, puis à
comparer les probabilités de résultats, les corrélations et les réponses aux
interventions avec des observations indépendantes. Une équation de champ
fournie reste une entrée même lorsque sa simulation fonctionne, et les choix de
construction n’ont besoin d’aucun théorème d’unicité, seulement de conséquences
qui résistent aux tests. Les contributions doivent résoudre une connexion
nommée ou révéler un défaut du modèle commun. Le
[guide de construction](model/CONSTRUCTION.md), le
[registre des candidats](model/CANDIDATES.md) et les
[tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
décrivent ce travail.

Les publications d’articles suivent la
[politique de publication](model/CONSTRUCTION.md#research-programme-and-paper-release-hold).
Les PDF de la branche principale sont des aperçus des sources ; les versions
déposées conservent leurs octets d’origine.
