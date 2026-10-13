# Holographie des parcelles d’observateur

**La physique à partir des registres que les observateurs peuvent partager.**

[Article de physique](https://philpapers.org/rec/MUEFOC) · [Pragma Research](https://floatingpragma.io/research/) · [English](README.md)

L’Holographie des parcelles d’observateur (OPH) demande si des observateurs
bornés, qui comparent leurs registres et réparent leurs désaccords, peuvent
reproduire les expériences que nous observons. Le but est un observateur
exécutable dont les mesures et les interprétations partagées concordent avec
les nôtres.

**Dans l’image proposée par OPH, la réalité est fondamentalement
computationnelle et l’observation est première.** La réalité est une structure
intemporelle et autoréférentielle d’observations, de registres et de leurs
relations. Chaque observateur donne un sens à ses observations, et le sens
partagé doit résister aux comparaisons avec les autres observateurs. Il n’existe
pas de continuum spatio-temporel fondamental, indépendant des observateurs, qui
contiendrait ces expériences. L’espace, le temps, la matière et les objets
classiques en sont des interprétations effectives. Le continuum apparent
appartient à cette description ; la « réalité objective » désigne la structure
opérationnelle stable sur laquelle les observateurs s’accordent.

Les observateurs vivent le changement depuis l’intérieur de cette histoire
complète. Sa cohérence intemporelle et son processus interne et local
d’observation et de réparation sont deux descriptions de la même structure
proposée. L’horloge d’un observateur est une lecture interne à cette structure,
sans temps extérieur dans lequel toute l’histoire devrait se dérouler. Le
[guide des termes](docs/INTERPRETATION_AND_TERMS.md) explique cette vision.

**L’Univers se lit lui-même jusqu’à l’existence.** OPH part d’une seule
hypothèse : ses observateurs retrouvent la spécification de leur propre monde à
partir de leurs registres et construisent un matériel dont l’exécution réalise
cette même histoire, y compris ses constructeurs. Comme rien ne se trouve hors
de cette structure, ce que les observateurs retrouvent et mettent en œuvre n’a
besoin d’aucune autre justification. Le support, l’algèbre de l’observateur et
la loi de réparation sont justifiés parce qu’ils sont retrouvés, qu’ils
fonctionnent et qu’ils reproduisent ce que mesurent les observateurs. Rien de
non retrouvé ne peut être supposé, et le registre public, ce qui résiste à
toute comparaison entre observateurs, est ce que signifie « notre Univers ».
La spécification exécutable est le point où le récit entre dans cette boucle ;
la boucle elle-même n’a pas de premier pas. La [base](docs/BASIS.md) l’énonce
avec précision.

Cette approche radicale, qui place les observateurs en premier, s’attaque à des
problèmes fondamentaux qui dépassent l’image classique : la mesure quantique,
les corrélations incompatibles avec une description classique locale et
l’origine d’un espace-temps apparent. OPH cherche une explication commune dans
la cohérence des registres et des opérations des observateurs. Expliquer ces
phénomènes par une seule construction est l’ambition du programme ; les
résultats ci-dessous en établissent des parties précises.

Le projet dispose de résultats mathématiques finis vérifiés et d’observateurs
exécutables qui en relient plusieurs dans un même processus. **L’accord avec
l’ensemble des statistiques de mesure observées n’est pas établi.** Les travaux
ci-dessous donnent des éléments concrets de cette construction et testent
comment ils s’assemblent.

## Observateur exécutable

Le [répertoire du modèle](model/README.md) contient l’implémentation, les tests
et les résultats reproductibles. Son expérience la plus intégrée place un champ
scalaire, des fermions chargés et un flux électrique sur un même registre
quantique fini. Une observation produit un registre, ce registre commande une
rétroaction, et une observation ultérieure commande une réparation locale.
Chaque résultat est conservé, et une borne d’erreur exacte certifie que la
rétroaction change la réponse mesurée. Pour la même action de champ déclarée,
un calcul exact à occupation infinie maintient le contraste de rétroaction loin
de zéro avec le détecteur de demi-droite en coordonnée ; l’action et la
préparation restent des entrées fournies.

Les mêmes histoires fournissent des espaces de Hilbert d’observateurs et des
sondes d’horloge modulaire lisibles. Les courants natifs donnent des relations
exactes de conservation et de transport, avec des lectures géométriques liées
au support à douze ports. Un détecteur lisse a des erreurs contrôlées lorsque
la coupure d’occupation à deux sites augmente. Ce sont des conséquences
exécutées et démontrées de l’appareil spécifié ; son hamiltonien de champ, sa
préparation et certains choix d’horloge sont des entrées fournies. Le
[processus commun](model/common_source/README.md) donne les résultats et leurs
limites, y compris la distinction entre relaxation de la réparation et ondes
dans le temps d’un observateur.

Essayez l’observateur fini depuis la racine du dépôt avec Python 3.11 ou plus
récent :

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

Le [paquet de prédictions](model/predictions/README.md) donne des prédictions
concurrentes pour des détecteurs et un test de rejet ; aucun comptage
expérimental indépendant n’est fourni. Le modèle contient aussi un test de
source fini qui distingue deux instruments de réparation par leurs lois de
registres conservés, avec une préparation et une lecture cohérentes
spécifiées. L’accord d’une simulation avec une loi de champ fournie n’établit
pas qu’OPH dérive les statistiques d’observation correspondantes.

## Fondements mathématiques

- **L’accord et les registres stables.** Sous des conditions explicites de
  terminaison et de cohérence, les réparations locales atteignent un résultat
  public commun, indépendant de leur ordre. D’autres résultats contrôlent la
  conservation des registres et leurs changements lors du raffinement.
  [Article sur le consensus](paper/reality_as_consensus_protocol.pdf).
- **Les observations quantiques.** Dans l’algèbre quantique finie énoncée, des
  résultats vérifiés couvrent les mesures répétables, les registres classiques
  protégés, l’incertitude et la borne de corrélation de Bell, avec un état qui
  l’atteint. [Algèbres finies d’événements quantiques](extra/machine_checked_finite_event_algebras.pdf).
- **Trois directions spatiales.** Le support fixe à douze ports et sa règle de
  réparation moyenne donnent une réponse lente exactement tridimensionnelle et
  sa géométrie euclidienne. C’est une brique géométrique pour l’interprétation
  spatiale d’un observateur. [Construction du support](code/a5_closure/README.md).
- **La symétrie de jauge.** Avec une réponse réversible complète et des
  hypothèses de transport entre observateurs, la construction impose l’algèbre
  de Lie de jauge du Modèle standard. Une représentation de matière spécifiée
  donne une génération sans anomalies.
  [Dérivation de jauge](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **La thermodynamique.** Une construction finie commune de réparation donne
  les quatre lois, les relations de fluctuation et la borne d’effacement de
  Landauer sous des hypothèses énoncées de référence, de calibration de
  l’énergie et de basse température.
  [Exposé thermodynamique](paper/observers_are_all_you_need.pdf).
- **La capacité d’information.** Des résultats exacts comptent les registres
  qui restent distinguables après les opérations permises, bornent ce nombre
  par le support disponible et le suivent lors du raffinement.
  [Capacité des registres](code/capacity_readback/README.md).
- **La matière et les champs.** Une action fournie de scalaire chargé couplé au
  champ de Maxwell, une géométrie et une quantification donnent un hamiltonien
  en interaction après réduction de jauge. Un raffinement scalaire contrôlé
  relie les lectures classiques d’un détecteur à une équation de champ non
  linéaire. [Continuum de matière](paper/observers_are_all_you_need.pdf).

Ce sont des résultats mathématiques conditionnels. Faire tenir leurs prémisses
pour les observations d’un même modèle commun est le problème central de la
construction. La [base](docs/BASIS.md) explique comment les prémisses des
théorèmes se rapportent à la spécification exécutable, et la
[carte modèle–théorèmes](model/PREMISES.md) indique quelles prémisses la
construction fournit. 

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
  donnent chaque article et son rôle. La relation de Koide tient exactement
  sous une prémisse d’équilibre déclarée, et le programme de capacité demande
  si la capacité publique attribuée à l’univers s’accorde avec celle
  reconstruite depuis l’intérieur. Ces calculs de constantes, comme celui de
  la structure fine, ont un statut diagnostique.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Explorer et reproduire

Commencez par l’[article principal](https://philpapers.org/rec/MUEFOC),
l’[index des articles](paper/), le [livre](book/) ou le
[guide de reproduction](REPRODUCE.md). La [bibliothèque Lean](Lean/) contient
plus de 12300 théorèmes et lemmes publics sans preuve admise. Le
[registre des postdictions](docs/POSTDICTION_LEDGER.md) et le
[programme de falsification](docs/OPH_FALSIFICATION_PROGRAM.md) distinguent les
comparaisons avec des mesures des tests des constructions mathématiques.

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
interventions avec des observations indépendantes. Une limite continue
contrôlée signifie que ces observations admettent une description continue de
plus en plus précise. Elle n’affirme pas qu’un continuum spatio-temporel existe
sous elles. La voie conditionnelle vers la dynamique d’Einstein a le même but :
que les mesures des observateurs concordent avec la description gravitationnelle
effective.

Utilisez le support et les règles existants pour fournir les entrées par une
preuve ou une exécution. Une équation de champ fournie reste une entrée même
lorsque sa simulation fonctionne. Les choix de construction n’ont besoin
d’aucun théorème d’unicité, mais leurs conséquences doivent résister aux tests.
Les contributions doivent résoudre une connexion nommée ou révéler un défaut du
modèle commun. Le [guide de construction](model/CONSTRUCTION.md), le
[registre des candidats](model/CANDIDATES.md) et les
[tâches de recherche](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
décrivent ce travail.

Les publications d’articles suivent la
[politique de recherche](model/CONSTRUCTION.md#research-programme-and-paper-release-hold).
Les PDF de la branche principale sont des aperçus des sources ; les versions
déposées conservent leurs octets d’origine.
