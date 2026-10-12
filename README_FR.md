# Holographie des parcelles d’observateur

**La physique à partir des registres que les observateurs peuvent partager.**

[Article de physique](https://philpapers.org/rec/MUEFOC) · [Pragma Research](https://floatingpragma.io/research/) · [English](README.md)

L’Holographie des parcelles d’observateur (OPH) étudie si la physique familière
peut être reconstruite à partir d’observateurs bornés qui parviennent à un
accord. Une parcelle d’observateur possède un état local, une frontière munie
de ports, une capacité de relecture, des registres et une rétroaction qui
répare les désaccords. Dans ce modèle, les faits deviennent publics lorsque
les registres résistent à la comparaison entre parcelles qui se recouvrent.

L’univers proposé est une même structure, décrite comme une histoire cohérente
intemporelle et comme le calcul local interne qui la réalise. Ses registres
et ses réponses aux interventions constituent ce fonctionnement, sans ordinateur
extérieur ni horloge globale privilégiée. L’espace, la durée physique et la
matière sont les lectures émergentes recherchées. Les constructions actuelles
fournissent encore certaines géométries, actions et correspondances d’horloges
comme prémisses. Le [guide des termes](docs/INTERPRETATION_AND_TERMS.md)
distingue cette proposition des résultats conditionnels établis.

Ce dépôt contient le travail scientifique derrière cette idée : articles,
preuves vérifiées par machine et modèles exécutables. Il relie l’accord fini
aux probabilités quantiques, à la géométrie spatiale et à la structure de
symétrie des forces connues, avec une voie conditionnelle vers la dynamique
d’Einstein. Chaque construction précise ses hypothèses. Établir une
réalisation physique commune exige des mesures qui distinguent les modèles
proposés.

## État de la recherche

Nous disposons de **résultats mathématiques conditionnels et de candidats
finis concrets pour la loi de réparation**. Ces candidats reproduisent et
relient certaines structures requises ; des simulations permettent d’en
distinguer plusieurs. **Aucun modèle commun n’est établi comme physiquement
correct.** Le [registre des candidats](model/CANDIDATES.md) et le
[suivi des modèles](model/CONSTRUCTION.md) distinguent les résultats construits,
les hypothèses fournies et les liens physiques à établir.

**La validité physique d’OPH est une hypothèse à tester.** La cohérence
logique, les théorèmes conditionnels et les simulations réussies ne prouvent
pas que la nature suit OPH ou une loi de réparation particulière. Nous
conservons les lois alternatives et révisons ou rejetons les interprétations
physiques que les tests mettent en défaut.

Un test physique exige une identification opérationnelle des états, des
registres, du temps et des paramètres, une observable discriminante, un bilan
de calibration et d’incertitude, puis une prédiction et un critère de rejet
fixés avant l’examen de données indépendantes. La comparaison inclut la
physique établie et les modèles alternatifs. Le
[programme de tests physiques](https://github.com/FloatingPragma/observer-patch-holography/issues/1025)
et le [registre des prédictions figées](docs/FROZEN_PREDICTION_LADDER.md)
suivent ce travail. Reproduire une valeur utilisée pour choisir le modèle
constitue un ajustement, pas un test indépendant.

La [recherche d'obstructions à l'échelle du programme](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
examine en parallèle les contradictions, les impossibilités de réalisation
commune, les entrées physiques indépendantes inévitables et les exclusions
empiriques. Chaque résultat précise la formulation ou la classe de modèles
qu'il rejette. Un résultat négatif décisif est utile ; ne pas en trouver ne
confirme pas OPH. Cette tâche distingue aussi une réfutation scientifique d'une
décision pratique d'arrêter les recherches.

## Priorité de recherche : un même processus de réparation

**Trouver la règle de réparation commune la plus simple qui satisfait les
contraintes connues.**

**Le programme n’est achevé que lorsqu’une telle règle commune a été trouvée
et que sa compatibilité avec les contraintes requises est établie dans un même
modèle.** Les preuves partielles, simulations et améliorations du compilateur
sont des étapes ; aucune ne suffit à clore le programme. L’identification
physique exige aussi des observations discriminantes.

La priorité est une loi microscopique explicite dont l’évolution fournit les
entrées des théorèmes existants dans un même modèle. Des résultats conditionnels
peuvent être corrects sans que leurs prémisses physiques soient réalisées
ensemble. Nous partons de la construction compatible la plus simple,
réutilisons les résultats démontrés et confrontons chaque hypothèse ajoutée
à des alternatives. La géométrie du support, la réponse de jauge réversible
et les horloges opérationnelles locales doivent être des lectures compatibles de ce même
processus. Tout mécanisme ajouté doit répondre à un lien manquant démontré,
et toute hypothèse physique ajoutée doit être déclarée. L’échec d’un candidat
purement scalaire ne prouve pas qu’une loi commune simple soit impossible.

Le [registre des lois candidates](model/CANDIDATES.md) compare leurs
propriétés testées et leurs hypothèses. Le [suivi des modèles](model/CONSTRUCTION.md), la
[carte des prémisses](model/PREMISES.md) et les
[tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
décrivent ce travail. Un lien testé ou un candidat rejeté constitue un progrès ;
le nombre de théorèmes ne mesure pas l’achèvement physique.

Les publications d’articles et leur restructuration générale sont suspendues
jusqu’à une étape cohérente de construction du modèle. Audits, corrections
et vérifications de reproductibilité continuent selon la
[politique de recherche](model/CONSTRUCTION.md#research-programme-and-paper-release-hold).

Les PDF de la branche principale sont des versions de travail et peuvent
conserver l’identifiant de publication existant après correction. Les fichiers
associés à une version Git ou à un dépôt conservent leurs octets d’origine ;
une version de travail n’est pas une nouvelle publication.

## Commencer ici

| Pour explorer… | Commencer par… |
| --- | --- |
| Le modèle en construction et ses liens physiques à établir | [Suivi des modèles](model/CONSTRUCTION.md) |
| Les lois candidates, leurs hypothèses et leurs tests distinctifs | [Registre des candidats](model/CANDIDATES.md) |
| L’argument physique existant | [*Finite Observer Consensus as a Reconstruction Principle*](https://philpapers.org/rec/MUEFOC) |
| Les preuves et les calculs reproductibles | [Bibliothèque Lean](Lean/) et [guide de reproduction](REPRODUCE.md) |
| Les travaux associés sur les machines apprenantes | [Cadence](https://floatingpragma.io/cadence/) et [démonstrations interactives](https://floatingpragma.io/demos/) |
| Les explications pour le grand public | [Blog Pragma Research](https://blog.floatingpragma.io/) |

## Fondements mathématiques

Le programme de recherche relie ces domaines de la physique par une même
architecture d’observateurs, avec des hypothèses distinctes pour chaque résultat.

- **L’accord entre observateurs.** Sous les conditions de terminaison et de
  cohérence énoncées, les réparations locales aboutissent à un registre public
  protégé, indépendamment de leur ordre. Les théorèmes finis décrivent
  l’accord, la stabilité et le raffinement.
- **Les probabilités quantiques.** Les algèbres finies d’événements portent la
  règle de Born, le conditionnement de Lüders et des résultats exacts de
  Tsirelson, dont une construction de Bell qui atteint 2√2. Ces résultats
  mathématiques concernent la branche d’algèbre d’événements déclarée.
- **La géométrie et les forces connues.** Une construction de réponse à douze
  ports fournit une lecture spatiale tridimensionnelle. La réponse réversible
  complète et le transport interne entre observateurs imposent le type de Lie
  de jauge du Modèle standard sous les hypothèses énoncées. Une représentation
  de matière fournie donne une génération à quinze états sans anomalies ; sa
  réalisation physique et le groupe de jauge global exigent une structure
  supplémentaire.
- **L’espace-temps et la gravité.** Des lois de source et de lecture déclarées
  permettent des constructions contrôlées de géométrie causale. L’implication
  vers l’équation d’Einstein suppose une réalisation physique commune avec les
  données énoncées de contrainte mécanique, d’entropie, de continuum et
  d’échelle.
- **La matière classique et quantique.** Une action fournie de champ scalaire
  chargé couplé au champ de Maxwell permet un mouvement continu non linéaire
  contrôlé dans son secteur réel et un espace d’états quantiques en interaction
  sur un maillage fixé. Sélectionner cette action à partir des histoires des
  observateurs et l’identifier à la matière physique sont des exigences
  distinctes.
- **Les constantes comme problèmes de point fixe.** La relation de Koide tient
  exactement sous une prémisse d’équilibre déclarée. Le calcul de structure
  fine certifie une racine d’une application de fermeture déclarée et possède
  un statut diagnostique. Le programme de capacité demande si la capacité
  publique attribuée à l’univers s’accorde avec celle reconstruite depuis
  l’intérieur. Relier ces constructions aux constantes mesurées exige une
  identification physique.

L’[article principal](https://philpapers.org/rec/MUEFOC) développe ces liens
et leur portée précise. L’[index des articles](paper/) donne les exposés
spécialisés, notamment sur la thermodynamique, la géométrie des champs et les
constructions de points fixes pour les constantes.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Preuves et éléments de vérification

La [bibliothèque Lean](Lean/) contient plus de 12300 théorèmes et lemmes
publics sans preuve admise. Des certificats exacts et des simulations
reproductibles accompagnent les arguments mathématiques. Le
[simulateur de physique](https://github.com/muellerberndt/oph-physics-sim)
associé fournit une dynamique exécutable des observateurs et les éléments
de vérification conservés.

La [référence des axiomes](docs/AXIOM_REFERENCE.md) énonce les trois axiomes
fondamentaux ; le [registre des prémisses](docs/registers/PREMISE_REGISTER_V3.md) et la
[carte des prémisses](model/PREMISES.md) consignent les prémisses
supplémentaires et les hypothèses propres à chaque résultat, y compris les
identifications physiques et les données empiriques.
Le [registre des postdictions](docs/POSTDICTION_LEDGER.md)
consigne les comparaisons avec des valeurs mesurées et l’origine de leurs
entrées ; l’[échelle des prédictions gelées](docs/FROZEN_PREDICTION_LADDER.md)
consigne les tests dont les conditions doivent être fixées avant l’examen
des données de comparaison. Le [programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md)
donne les observations qui réfuteraient des affirmations précises.

### Reproduire le noyau fini

Après avoir installé les dépendances décrites dans [REPRODUCE.md](REPRODUCE.md),
vérifiez le graphe des affirmations et une sélection de résultats sur
l’algèbre finie, la capacité des registres et le consensus :

```bash
python3 tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

Le guide de reproduction décrit les vérifications plus larges et la portée
de chaque famille de preuves.

## Carte de reconstruction

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="Chaîne de reconstruction OPH" width="92%">
  </a>
</p>

<p align="center"><sub>La carte de reconstruction OPH relie registres d’observateurs, géométrie spatiale à trois dimensions, ordre causal, horloges, champs et états quantiques. Chaque flèche représente un lien mathématique ; les articles précisent les hypothèses qui les réunissent dans une description physique effective.</sub></p>

## Guide du dépôt

| Chemin | Contenu |
| --- | --- |
| [`flagship/`](flagship/) | Article principal autonome de physique, source TeX et PDF |
| [`paper/`](paper/) | Articles principaux et index des publications |
| [`Lean/`](Lean/) | Développement mathématique vérifié par machine |
| [`code/`](code/) et [`evidence/`](evidence/) | Modèles exécutables, certificats et éléments de reproduction |
| [`extra/`](extra/) et [`cosmology/`](cosmology/) | Articles de recherche mathématique et physique spécialisés |
| [`book/`](book/) | *Reverse Engineering Reality*, source et livre téléchargeable |
| [`docs/`](docs/) | Références canoniques : axiomes, prédictions gelées, postdictions et critères de réfutation |
| [`docs/registers/`](docs/registers/) | Registres générés des prémisses, observations, instruments et sélections |
| [`docs/research/`](docs/research/) | Notes de recherche en Markdown, chacune vérifiée par un module de code |
| [`docs/policies/`](docs/policies/) | Politiques de données, de preuves matérielles et de rédaction |

## Contribuer

Une contribution corrige un défaut existant ou résout une incertitude précise
du modèle commun et de ses tests. Réutilisez les résultats démontrés ; ajoutez
une preuve lorsqu’une obligation du modèle l’exige. L’expansion autonome des
théorèmes finis et les objectifs de comptage sont hors programme. Commencez par
les [tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
et le [guide de reproduction](REPRODUCE.md). Les journaux exploratoires restent
locaux ; les décisions publiques indiquent leurs hypothèses et leurs preuves.
Les notes vérifiées vont dans [`docs/research/`](docs/research/), et les articles
respectent la suspension décrite ci-dessus.

## Des parcelles d’observateur aux machines apprenantes

[Cadence](https://github.com/muellerberndt/cadence) rend l’idée de stabilisation
exécutable sous forme d’architecture d’apprentissage. Ses parcelles logicielles
bornées portent un état local et une capacité de relecture, avec une
rétroaction qui répare les erreurs de prédiction. Son cerveau par
défaut, le Système 1, apprend par l’expérience comme les cerveaux animaux, et
non par rétropropagation. Des observateurs optionnels, le Système 2, ajoutent
une rétroaction récursive au sein de la même stabilisation.

L’[article Cadence](https://philpapers.org/rec/MUECAP-2) décrit l’architecture
et les expériences. Son évaluation repose sur le comportement
d’apprentissage et les ressources mesurées.
[Pragma Research](https://floatingpragma.io/) relie ce travail à l’IA incarnée.

## Licence

Le dépôt utilise des licences séparées par type d’artefact. Tout le logiciel, y compris la bibliothèque Lean, [`code/`](code), [`tools/`](tools) et les schémas de [`schemas/`](schemas), est publié sous [Apache-2.0](code/LICENSE). Les articles, le livre, la documentation, les figures, les données, les registres générés de [`tracking/`](tracking) et les données de particules empaquetées de [`pdg_data/`](pdg_data) sont publiés sous [CC BY-NC-SA 4.0](LICENSE). Les valeurs tabulées de `pdg_data/` conservent les conditions amont du Particle Data Group. Le fichier [LICENSE](LICENSE) donne la carte par répertoire.
