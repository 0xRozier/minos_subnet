# Minos SN107 — Observations et plan pour le top 20

*Rédigé le 2026-10-11, mis à jour le 2026-10-13 à 13 h 45 UTC (§3.3, §3.5, §3.6, §5, §8). Les chiffres viennent de runs réels : GATK 4.5.0.0 →
hap.py/vcfeval → score v2, le même pipeline que les validateurs, sur des rounds
réseau récents téléchargés via `/verification`.*

## 1. L'objectif, chiffré

- Paiement (`/scoring/network-config`) : `winner_weight=0.8`, `dust_top_n=20`,
  `dust_decay=0.9`. Seuls les 20 premiers d'un round sont payés.
- **Le top 20 est collé au top 1.** Sur les 8 rounds mesurés, l'écart top 1 →
  rang 20 va de 0,01 à 0,2 point (ex. 84,93 vs 84,92). Un grand cluster de
  miners produit des VCF quasi identiques. Entrer dans le top 20, c'est donc
  **égaler ce cluster**, pas s'en approcher.
- Rotation actuelle : **chr12 40 %, chr13 40 %**, 2,2 % pour chacun des autres
  chromosomes. Ce qui compte, c'est la performance sur chr12/chr13.

## 2. Ce que mesure le score (v2)

```
score = 100 × (0,70 × core + 0,30 × germline)     (0 si le gate ti/tv ou het/hom échoue)
core     = F1 par classe pondéré : indel_1bp 0,40 · snp_het 0,18 · indel_2_3bp 0,16
           · indel_8bp 0,14 · indel_4_7bp 0,10 · snp_hom 0,02
germline = exp(−fp_per_target / 8)
fp_per_target = tous les FP de la région / nombre de mutations synthétiques
```

Sur les vrais rounds, il y a 500 à 750 mutations synthétiques, donc
`fp_per_target` vaut 1,4 à 5. En practice, il n'y a que 90 à 200 mutations :
le germline y est écrasé et les scores practice ne sont **pas** comparables
au réseau.

## 3. Observations

### 3.1 La plateforme ne fait pas tourner ce que tu as mesuré

`/scoring/parameter-ranges` impose des bornes plus étroites que
`templates/tool_params.py`. Une valeur hors bornes est **remplacée
silencieusement par le défaut**. Dans tes configs :

| Config | Valeur testée | Ce qui tourne en live |
|---|---|---|
| `soft_thr65`, `soft_base` | `dont_use_soft_clipped_bases=true` | `false` (seule valeur autorisée) |
| `soft_thr65`, `soft_base`, `ei_best` | `max_num_haplotypes_in_population=201` | 128 |
| `ei_best` | `active_probability_threshold=0.0003` | 0,002 (min 0,001) |
| `ei_balanced` | `max_reads_per_alignment_start=864` | 50 (max 300) |

Conséquence : le « seul levier confirmé » du journal (soft clips, +0,0015)
n'existe pas en live, et les scores du journal pour ces configs sont
optimistes. `python scripts/optimize_config.py check --conf …` détecte ces
cas.

### 3.2 Pièges qui font perdre un round entier

| Valeur | Effet |
|---|---|
| `base_quality_score_threshold` < 6 | GATK refuse et plante, le round vaut 0 (la plateforme accepte 0–50) |
| `base_quality_score_threshold=35` (dans `ei_best`, `ei_balanced`) | 0 variant sur chr14/16/17 (constaté dans le journal) |
| `emit_ref_confidence=GVCF` | ~900 000 lignes, hap.py échoue, le round vaut 0 (testé) |

### 3.3 Mesure sur 8 vrais rounds chr12/chr13 (version live)

Rounds du 10/10 16:42 au 11/10 02:18. Les 2 derniers (chr12:6–14M,
chr13:37–45M) sont mis de côté pour la validation finale.

Toutes les configs sont mesurées **telles que les validateurs les exécutent**,
c'est-à-dire avec les valeurs hors bornes remises à leur défaut.

| Config | Moyenne (8 rounds) | Écart moyen vs top 1 |
|---|---|---|
| **`soft_thr65`** (ta prod) | **80,42** | **−1,82** |
| `journal_best` (ma reconstruction de `soft_thr65`) | 80,42 | −1,82 |
| `soft_thr65` + `contamination_fraction_to_filter=0.05` | 80,35 | −1,89 |
| `soft_base` (seuil 42,7) | 80,24 | −2,00 |
| `ei_best` | 79,53 | −2,71 |
| `gatk_defaults_plus` | 77,13 | −5,11 |
| `ei_balanced` | 76,57 | −5,67 |
| Défauts GATK | 75,99 | −6,25 |

- `soft_thr65` est ta meilleure config. Ma reconstruction à partir du journal
  donne exactement le même score : les écarts d'arrondi (adaptive pruning,
  indel_het, pruning_lod) n'ont aucun effet.
- **La config dédiée à chr12 n'apporte rien.** Sur les 4 rounds chr12,
  `soft_thr65` fait 78,93 contre 78,77 pour `soft_base` (sur chr13 : 81,90 contre
  81,72). Le gain de `soft_base` sur chr12 dans le journal venait d'un seul
  round, mesuré avec les soft clips activés.
- `gatk_defaults_plus` faisait 0,8946 sur un seul round chr19 dans le journal ;
  ici elle n'est qu'à 77,13. C'est un autre exemple de mesure trompeuse sur
  un seul round.

Détail de `soft_thr65` :

| Round | Score | Top 1 | Écart | Core | Germline | FP/cible |
|---|---|---|---|---|---|---|
| chr12:33–45M | 76,80 | 78,14 | −1,34 | 0,868 | 0,535 | 5,00 |
| chr13:97–111M | 77,37 | 79,30 | −1,93 | 0,857 | 0,581 | 4,35 |
| chr12:57–71M | 77,28 | 79,53 | −2,25 | 0,836 | 0,626 | 3,74 |
| chr13:97–110M | 83,20 | 84,93 | −1,73 | 0,868 | 0,747 | 2,33 |
| chr12:47–53M | 82,75 | 84,98 | −2,23 | 0,840 | 0,798 | 1,80 |
| chr13:42–57M | 80,06 | 82,22 | −2,16 | 0,843 | 0,702 | 2,83 |
| chr12:6–14M | 78,89 | 80,50 | −1,61 | 0,846 | 0,656 | 3,37 |
| chr13:37–45M | 86,97 | 88,32 | −1,35 | 0,882 | 0,841 | 1,39 |

Ce qu'on en tire :

1. **L'écart est constant (−1,3 à −2,3 sur chaque round).** Ce n'est pas la
   faute d'une région difficile : c'est une différence systématique de config.
2. **Ta config est déjà bonne** : +4,4 points sur les défauts, et elle divise
   par deux les FP.
3. **Le germline est le terme le plus faible partout** (0,54 à 0,84). Pour
   combler 1,8 point, il suffirait par exemple de baisser `fp_per_target`
   d'environ 25 % à core égal, ou de gagner environ 0,026 de core à FP égaux.
4. **Les FP qui pèsent dans le core** (comptes TP/FN/FP de `journal_best`) :
   - `indel_8bp` a une précision d'environ 50 % : autant de FP que de TP sur
     certains rounds (145 FP pour 187 TP sur chr12:33–45M). Avec un poids de
     0,14, c'est la classe la plus mal calibrée.
   - `snp_het` a entre 336 et 1605 FP par round, alors que les `snp_hom` en ont
     presque 0. Le bruit SNP est hétérozygote.
   - `indel_1bp` (poids 0,40) a un F1 de 0,85 à 0,91, avec à peu près deux
     fois plus de FP que de FN.
5. **Un essai sur un seul échantillon trompe.** La contamination à 0,05 donnait
   +0,2 sur un échantillon practice, puis −0,07 sur les 8 vrais rounds. Il en va
   de même pour les configs par chromosome du journal : leurs écarts
   (≈0,002) sont plus petits que la variation d'un round à l'autre sur un même
   chromosome.

### 3.4 Ce que le journal a déjà établi (et qui reste vrai)

- ~800 configs testées, sur les 22 paramètres. Le paysage est plat : aucun
  paramètre isolé ne bouge le score de plus de quelques millièmes.
- Leviers réels : `min_pruning` (3–5), `standard_min_confidence_threshold_for_calling`
  (60–75 sur chr13), `heterozygosity` (0,007), `min_dangling_branch_length` (≈20),
  `min_mapping_quality_score` (≈50), `pcr_indel_model=AGGRESSIVE`.
- Paramètres sans effet : `max_alternate_alleles`, `max_num_haplotypes`,
  `indel_heterozygosity`, `adaptive_pruning`, `pruning_lod` (2,3–5),
  `max_reads_per_alignment_start`, `recover_all_dangling_branches`.

### 3.5 Phase 2 — variations ciblées autour de `soft_thr65`

Moyenne sur les 6 rounds de tuning. `soft_thr65` y fait 79,577 : c'est plus
bas que 80,42 parce que les 2 rounds mis de côté, plus faciles, ne sont pas
comptés.

| Variation (le reste = `soft_thr65`) | Moyenne | Écart |
|---|---|---|
| `min_mapping_quality_score=58` | **79,634** | **+0,057** |
| seuil 75 | 79,578 | +0,001 |
| `min_pruning=3` | 79,561 | −0,016 |
| MQ 45 | 79,510 | −0,067 |
| seuil 55 | 79,497 | −0,080 |
| `min_pruning=5` | 79,473 | −0,104 |
| seuil 42,7 (`soft_base`) | 79,416 | −0,161 |
| `heterozygosity=0.004` | 79,568 | −0,009 |
| `min_dangling_branch_length=15` | abandonnée après 2 rounds (77,0) | nettement pire |

- **`soft_thr65` est au sommet d'un plateau.** Aucune retouche d'un paramètre
  ne gagne plus de 0,06. Le seuil de confiance est neutre entre 55 et 75.
- MQ 58 fait légèrement mieux : elle gagne sur 4 rounds sur 6 (+0,23 au mieux)
  et coupe quelques faux positifs. C'est un gain marginal, pas suffisant pour
  combler l'écart de ~1,8 point.
- Conclusion provisoire : le chemin vers le top 20 ne passe pas par une petite
  retouche de `soft_thr65`. La suite de la recherche (combinaisons proposées
  par Optuna) et la phase 3 (hypothèse 4.2) prennent le relais.

### 3.6 Phase 2 — propositions d'Optuna (point à 22 h 30 UTC)

25 essais lancés : 13 complets, 11 abandonnés en route (pruning) parce qu'ils
décrochaient dès les premiers rounds.

**Nouveau meilleur : essai 15, 79,855 (+0,28 vs `soft_thr65`).**

| Paramètre | `soft_thr65` | Essai 15 |
|---|---|---|
| `min_pruning` | 4 | **2** |
| `standard_min_confidence_threshold_for_calling` | 65 | 71,7 |
| `heterozygosity` | 0,007 | **0,0010** |
| `min_dangling_branch_length` | 20 | 19 |
| `min_mapping_quality_score` | 51 | 54 |
| `min_base_quality_score` | 18 | 26 |
| `base_quality_score_threshold` | 18 | 24 |
| `phred_scaled_global_read_mismapping_rate` | 28 | 38 |
| `pair_hmm_gap_continuation_penalty` | 15 | 22 |
| `pruning_lod_threshold` | 2,89 | 5,03 |

Round par round : +0,07, +0,01, +0,52, −0,23, +0,73, +0,57. Il gagne sur 5
rounds sur 6, surtout par le core (jusqu'à +0,006) plus que par les FP.
Gain encore loin des ~1,8 point qui manquent, et pas encore validé sur les 2
rounds mis de côté : c'est un candidat, pas encore une nouvelle prod.

**Point à 2 h 10 UTC (39/40 essais) :** deux nouveaux meilleurs, de la même
famille que l'essai 15.

| Essai | Moyenne | vs `soft_thr65` | Par round | Différences clés |
|---|---|---|---|---|
| **35** | **79,906** | **+0,33** | +0,35, +0,41, +0,90, −0,42, +0,22, +0,53 | pruning 2, seuil 73, het 0,0021, dang 15, MQ 59, mbqs 30, bqst 14, phred 30, gap 17, plod 6,3 |
| 34 | 79,863 | +0,29 | −0,08, +0,12, +0,61, −0,15, +0,43, +0,79 | pruning 2, seuil 78, het 0,0014, dang 16, MQ 46, mbqs 27, bqst 28, phred 35, gap 26, plod 3,2 |
| 15 | 79,855 | +0,28 | +0,07, +0,01, +0,52, −0,23, +0,73, +0,57 | voir tableau ci-dessus |

Les trois ont en commun `min_pruning=2`, un seuil de 72 à 78,
`heterozygosity` autour de 0,001–0,002 et `min_base_quality_score` vers 26–30.
Tous perdent sur le même round (chr13:97–110M). Après le dernier essai, les
5 meilleurs sont rejoués sur les 2 rounds mis de côté.

**Résultat final (2 h 45 UTC) — validation sur les 2 rounds mis de côté.**
Les 5 meilleures configs de la recherche, comparées à `soft_thr65` sur
exactement les mêmes rounds :

| Essai | Réglage (6 rounds) | Mis de côté (2 rounds) | Les 8 rounds |
|---|---|---|---|
| 15 | +0,28 | **−0,08** | +0,19 |
| 34 | +0,29 | **−0,26** | +0,15 |
| 35 | +0,33 | **−0,42** | +0,14 |
| 38 | +0,20 | **−0,05** | +0,14 |
| 26 | +0,09 | +0,30 | +0,14 |

**Aucune config ne passe la règle de décision.** Les configs qui gagnent le
plus sur les rounds de réglage perdent sur les rounds mis de côté. Leur gain
de +0,3 était donc du bruit de sélection : en essayant 40 configs, certaines
ont forcément de la chance sur 6 rounds. `soft_thr65` reste la prod. La
meilleure de la recherche est gardée comme simple candidate dans
`configs/candidates/search_trial15.conf`.

Importance des paramètres selon Optuna : `base_quality_score_threshold` 25 %,
MQ 15 %, `min_dangling_branch_length` 12 %, `pair_hmm_gap_continuation_penalty`
11 %, `heterozygosity` 10 %, `min_base_quality_score` 10 %. Mais tous ces
effets restent à l'intérieur de ±0,4 point. L'écart avec le top 1 reste de
−1,3 à −2,0 sur chaque round, pour toutes les configs.

**Conclusion de la phase 2.** Elle confirme le plateau du journal, cette fois
sur de vrais rounds chr12/chr13, en version live et avec validation : le
réglage fin de ces 11 paramètres GATK ne comble pas l'écart. La suite passe
par la phase 3 (hypothèse 4.2) et par l'analyse des FP en excès.

Tendance observée pendant la recherche : Optuna remonte vers des configs avec `min_pruning` bas (2–3),
MQ haut (54–60) et `pruning_lod_threshold` haut (5–9). Le pruning du graphe
se fait alors par la vraisemblance plutôt que par le nombre de reads, ce qui
récupère des vrais indels sans laisser passer plus de FP.

### 3.7 Phase 3 — ablation vers les défauts (4 rounds, terminée)

`soft_thr65` (live) avec un paramètre remis à son défaut ; moyenne sur
chr12:33–45M, chr13:97–111M, chr12:57–71M, chr13:97–110M. Référence : 78,66.

| Remis à défaut | Moyenne | Écart |
|---|---|---|
| Groupe « sans effet » (adaptive pruning, indel_het, max_alt, max_reads, pruning_lod) | 78,66 | 0,00 |
| `min_assembly_region_size` 98 → 50 | 78,66 | 0,00 |
| `pair_hmm_gap_continuation_penalty` 15 → 10 | 78,66 | 0,00 |
| `max_assembly_region_size` 394 → 300 | 78,64 | −0,02 |
| `min_pruning` 4 → 2 | 78,56 | −0,10 |
| `heterozygosity` 0,007 → 0,001 | 78,55 | −0,11 |
| `min_dangling_branch_length` 20 → 4 | 78,55 | −0,11 |
| `min_base_quality_score` 18 → 10 | 78,46 | −0,20 |
| `phred_scaled_global_read_mismapping_rate` 28 → 45 | 78,46 | −0,20 |
| `pcr_indel_model` AGGRESSIVE → CONSERVATIVE | 78,12 | −0,54 |
| `min_mapping_quality_score` 51 → 20 | 77,27 | **−1,39** |
| Défauts GATK complets | 73,14 | −5,52 |

- **Hypothèse 4.2 rejetée** : aucun paramètre remis à défaut ne fait mieux.
- 8 paramètres de `soft_thr65` sont neutres (les 5 du groupe,
  `min_assembly_region_size`, `pair_hmm_gap_continuation_penalty`, et quasiment
  `max_assembly_region_size`). Les leviers utiles : MQ (de loin), puis
  `pcr_indel_model=AGGRESSIVE`, `min_base_quality_score`, `phred`.
- MQ 51 vs 20 : même core (−0,0008), mais +0,72 FP/cible. MQ agit
  uniquement sur les faux positifs.

### 3.8 Vérifications sur les leaderboards publics (12/10)

1. **Notre pipeline local est exact.** Sur chacun des 8 rounds, un miner du
   réseau a exactement notre score « défauts GATK » (ex. 65,760, 73,773,
   72,604…). Même chose pour `soft_thr65` sur 4 rounds (probablement ton
   miner). L'écart mesuré est donc réel, pas un artefact local.
2. **Le top 20 est 100 % GATK** sur les 8 rounds. Meilleur bcftools : −2 à
   −7 points, meilleur DeepVariant : −5 à −13 points. Changer d'outil n'est
   pas la solution.
3. **Le classement réel de `soft_thr65` est 55e–68e sur ~70.** Une
   soixantaine de miners font mieux, probablement avec une même config
   partagée. Le repo n'en recommande aucune.
4. **L'écart avec le top vient uniquement des faux positifs.** D'après
   `/scoring/detailed/{hotkey}` (2 rounds encore en historique) :

| Round | | Score | Core | Germline | FP région (SNP + indel) |
|---|---|---|---|---|---|
| chr13:97–111M | top | 79,30 | 0,8582 | 0,641 | 1 740 (824 + 916) |
| | `soft_thr65` | 77,37 | 0,8565 | 0,581 | 2 127 (~1 000 + ~1 120) |
| chr12:57–71M | top | 79,53 | 0,8374 | 0,697 | 1 426 (663 + 763) |
| | `soft_thr65` | 77,28 | 0,8356 | 0,626 | 1 850 |

   Même core (à 0,002 près), mais environ 20 % de faux positifs en moins,
   répartis entre SNP et indels. **La cible est maintenant précise : couper
   ~20 % des FP de la région sans toucher au core.** Ça contredit le
   diagnostic « gap de core » du journal, qui comparait des régions
   différentes.

### 3.9 Phases 3b–3c — leviers « faux positifs » : criblage puis plan factoriel

Sur les 2 rounds où les FP du top sont connus (chr13:97–111M, chr12:57–71M).
Référence `soft_thr65` : 77,32.

**Criblage un par un** : `max_reads_per_alignment_start=25`,
`contamination_fraction_to_filter=0.03`, `indel_heterozygosity=3e-5`,
`max_alternate_alleles=2`, `recover_all_dangling_branches=true` → 77,32 (aucun
effet, VCF quasi identique). `active_probability_threshold=0.01` → 77,31,
`assembly_region_padding=25` → 77,36. `active_probability_threshold=0.05` →
**59,33** (catastrophique : à ne jamais utiliser).

**Plan factoriel 2^(7-3) résolution IV** (16 configs, `scripts/make_factorial.py`).
Effets principaux (niveau haut − niveau bas) :

| Levier | Score | Core | FP/cible |
|---|---|---|---|
| `assembly_region_padding` 100 → 25 | +0,114 | +0,0013 | −0,009 |
| `active_probability_threshold` 0,002 → 0,01 | +0,070 | +0,0006 | −0,011 |
| les 5 autres | 0,00 à −0,006 | 0 | ~0 |

**Interaction** : les 4 configs qui combinent padding 25 **et** apt 0,01
(fac06/07/14/15) font **77,50, soit +0,18**, alors que chacun seul ne fait
presque rien. C'est exactement le genre d'effet qu'un criblage un par un ne
voit pas. (Contraste BC = +0,085, aliasé avec AE et DF, mais A, D, E et F n'ont
aucun effet, donc l'interaction vient de BC.)

**Validation sur les 8 rounds** (`configs/candidates/sc65_apt01_pad25.conf`) :
80,52 contre 80,42 pour `soft_thr65`, soit **+0,11**. Par round : +0,11,
+0,14, +0,23, −0,06, +0,17, +0,14 (réglage) et +0,07, +0,07 (mis de côté).
Il gagne sur 7 rounds sur 8, et aussi sur les 2 rounds jamais utilisés pour le
choisir. C'est le premier gain **confirmé** de toute l'étude, mais il reste
petit : l'écart avec le top passe d'environ −1,82 à −1,72.

Il ne passe pas le seuil de +0,3 de la règle de décision, mais il est
régulier et sans risque (même core, un peu moins de FP). Recommandation :
l'adopter comme nouvelle prod, en sachant qu'il ne change pas le classement à
lui seul.

### 3.10 Diagnostic des faux positifs (12/10 soir)

**Comparaison avec le top sur 5 rounds récents** (métriques publiques de
`/scoring/detailed`, notre `sc65_apt01_pad25` lancé sur les mêmes rounds) :
l'écart vient surtout des **FP indels**, +20 à +71 % vs le top (743 vs 464,
1 104 vs 644…). On trouve même un peu **plus** d'indels synthétiques que le
top (FN 8 vs 13). FP SNP : +10 à +25 %. het/hom et ti/tv quasi identiques.

**Autopsie** (`scripts/fp_autopsy.py`) : nos FP indels sont surtout des
appels **homozygotes à faible profondeur** (DP médiane 17 vs 23, GQ 65 vs 99,
55 % en 1/1), et **61 % sont dans des répétitions ≥ 10 pb** (homopolymères,
STR), contre 29 % des vrais indels. Beaucoup sont sans doute de vrais indels
hors des régions de confiance GIAB, mais ils comptent comme FP sur la région.

**Angle mort de la recherche** : `pcr_indel_model`, le paramètre qui règle
la pénalité des indels dans les répétitions, n'a été essayé qu'en
`AGGRESSIVE` sur les 42 essais Optuna.

**Test `HOSTILE`** (5 rounds) :

| Round | FP indels top / AGGR / HOST | Indels synth. trouvés top / AGGR / HOST | Score AGGR → HOST |
|---|---|---|---|
| chr13:21–29M | 464 / 743 / 495 | 192 / 197 / 182 | 84,62 → 84,22 |
| chr13:65–71M | 305 / 385 / 296 | 52 / 53 / 49 | 88,56 → 85,50 |
| chr13:26–40M | 841 / 1 079 / 801 | 192 / 198 / 176 | 81,75 → 78,77 |
| chr12:84–89M | 233 / 279 / 208 | 103 / 101 / 97 | 88,66 → 86,09 |
| chr12:118–130M | 644 / 1 104 / 725 | 110 / 115 / 102 | 78,68 → 78,39 |

`HOSTILE` ramène nos FP indels **au niveau du top**, mais perd trop de vrais
indels (score −1,9 en moyenne). Le top a les FP d'un `HOSTILE` avec presque
le rappel d'un `AGGRESSIVE` : il combine vraisemblablement `HOSTILE` avec des
réglages qui récupèrent du rappel. → Phase 3d : recherche Optuna dans la zone
`HOSTILE`.

**Résultat de la phase 3d (40 essais, 13 paramètres, terminée le 13/10 à
4 h UTC) : négatif.** Meilleur essai `HOSTILE` : 77,53 sur les 6 rounds de
réglage, contre 79,70 pour `sc65_apt01_pad25` (`AGGRESSIVE`). Les 8 variations
ciblées (pruning 2/3, seuil 45, dangling 10, gap 10, indel_het 0,0008, phred 45,
pruning 2 + seuil 50) restent toutes entre 77,2 et 77,4 ; 30 essais sur 40
abandonnés en route. Aucun réglage accessible ne récupère, en `HOSTILE`, les
vrais indels qu'il fait perdre. On ne reproduit donc pas le profil du top
(FP d'un `HOSTILE`, rappel d'un `AGGRESSIVE`) avec ces paramètres.

### 3.11 Peut-on égaler le top sur UN round ? (13/10)

Test décisif proposé : optimiser directement le vrai score d'un seul round,
sans validation. Si même ainsi on n'atteint pas le top, aucune stratégie
« une config par round » ne marchera avec ces paramètres.

Round : chr12:84681001-89797000 (5 Mb). Top 90,37, rang 20 90,24.

**Recherche 1** (100 essais, 22 paramètres, base `sc65_apt01_pad25`) :
plafond à **89,33** (5 meilleurs essais entre 89,25 et 89,33).

| | Score | Core | Germline | FP SNP/indel | Synth. trouvés SNP/indel |
|---|---|---|---|---|---|
| Top | 90,37 | 0,904 | 0,904 | 187 / 233 | 403 / 103 |
| Notre meilleur (#12) | 89,33 | 0,895 | 0,889 | 222 / 264 | 407 / 103 |

Mais le leaderboard de ce round montre **51 miners GATK sur 66 au-dessus de
89,75**, tous avec une seule config soumise. Une config fixe meilleure que
notre plafond « oracle » existe donc dans l'espace des paramètres. Nos
recherches, toutes parties de la même famille de configs, ne l'ont pas
trouvée.

**Les valeurs hors bornes sont-elles appliquées ?** `soft_thr65` brut (soft
clips `true`, 201 haplotypes) donne 77,42 / 77,17 / 82,62 / 86,64 sur 4
rounds, sa version remise à défaut 77,37 / 77,28 / 82,75 / 86,97. Les scores
réseau les plus proches sont 77,39 / 77,29 / 82,77 / 86,97 : ils collent à
la version remise à défaut. L'hypothèse « la plateforme remet à défaut »
tient.

**Recherche 2** (en cours) : même round, base = défauts GATK, 8 points de
départ de familles différentes (tes configs + les nôtres), 40 essais au
hasard avant que l'optimiseur se concentre, 250 essais.

## 4. D'où peut venir l'écart — hypothèses, de la plus à la moins testable

1. **Une combinaison que les tests un-paramètre-à-la-fois ne voient pas.** Le
   journal a surtout fait des ablations sur un seul round. Une recherche
   bayésienne multi-rounds explore les interactions (ex. `min_pruning` ×
   seuil × MQ) et moyenne le bruit de région. *Testable : phase 2.*
2. **Le cluster de tête fait tourner une config « remise à défaut ».** Beaucoup
   de miners soumettent sans doute des valeurs hors bornes. Leur config
   effective a alors les défauts sur ces paramètres (apt=0,002, haplotypes=128,
   reads=50, soft clips=false…). Les 31 miners identiques peuvent venir de là.
   *Testable : partir de la meilleure config et remettre chaque paramètre à son
   défaut, un par un, sur les 6 rounds (phase 3).*
3. **Un régime différent sur un paramètre peu exploré avec les leviers actuels**
   (`phred_scaled_global_read_mismapping_rate`, `pair_hmm_gap_continuation_penalty`,
   `base_quality_score_threshold`, `min_base_quality_score`) : inclus dans la
   recherche.
4. **Rien dans les paramètres GATK** ne permet de rattraper le cluster. Si la
   phase 2 et la phase 3 plafonnent sous le top, c'est la conclusion honnête.
   Il restera alors à chercher ce que fait le cluster autrement (outil,
   version), pas plus de configs.

## 5. Plan

| Phase | Quoi | Statut | Critère pour passer à la suite |
|---|---|---|---|
| 0 | Outillage : `optimize_config.py` (bornes live, cache, reprise, held-out), `fetch-rounds` (vrais rounds + score top 1/top 20) | ✅ fait | — |
| 1 | Mesurer tes configs exactes, en version live, sur les 8 rounds | ✅ fait (§3.3) : `soft_thr65` est la base | Choisir la meilleure comme base |
| 2 | Recherche bayésienne sur les 6 rounds de tuning, partant de `soft_thr65`. D'abord 8 variations ciblées (§3.5), puis les propositions d'Optuna, 40 essais (~20 h sur 4 cœurs) | ✅ fait (§3.6) : aucun gain confirmé sur les rounds mis de côté, `soft_thr65` reste la prod | Gain ≥ +0,3 point sur les 6 rounds **et** confirmé sur les 2 rounds mis de côté |
| 3 | Ablation « vers les défauts » de `soft_thr65` (hypothèse 4.2) : 10 paramètres remis à défaut un par un + 1 groupe (les paramètres « sans effet » du journal), sur 4 rounds (2 chr12, 2 chr13), puis confirmation sur 8 rounds de ce qui ne perd rien. Configs dans `configs/candidates/ablation/` | ✅ fait (§3.7) : hypothèse rejetée | Idem |
| 3b | Sonde « faux positifs » : 8 paramètres jamais poussés dans le sens qui coupe les FP (max_reads 25, apt 0,01/0,05, padding 25, max_alt 2, contamination 0,03, indel_het 3e-5, recover=true), sur les 2 rounds où les FP du top sont connus. Configs dans `configs/candidates/fp_probe/` | ✅ fait (§3.9) : seule l'interaction padding 25 × apt 0,01 ressort (+0,18 sur 2 rounds), confirmée sur 8 rounds : +0,11, 7 rounds sur 8 | Couper ≥ 10 % des FP à core égal |
| 3d | Recherche Optuna en `HOSTILE` (base `sc65_hostile`) sur 13 paramètres de rappel/précision, 8 variations ciblées en tête de file, 8 rounds dont 2 mis de côté (`JOB=hostile`) | ✅ fait : négatif, meilleur `HOSTILE` 77,53 vs 79,70 (§3.10) | Gain confirmé sur les rounds mis de côté |
| 4 | Déploiement : une seule config globale (pas de configs par chromosome), passée au `check`, dans `configs/gatk.conf`. Supprimer les `configs/gatk_chr*.conf` | à venir | `check` OK |
| 5 | Suivi live : comparer ton score réel au score prédit sur 5 rounds (leaderboard public). Un écart > 0,2 point veut dire que le pipeline local diverge et doit être expliqué avant d'aller plus loin | à venir | Prédiction fiable |
| 6 | Entretien : une fois par semaine, `fetch-rounds` sur des rounds neufs, puis `baseline` sur la config en prod pour détecter une dérive | continu | — |

**Règle de décision.** Une config ne remplace la prod que si elle gagne en
moyenne sur les 6 rounds de tuning **et** sur les 2 rounds mis de côté.
Un gain sur un seul round ou un seul chromosome ne compte pas.

## 6. Ce qui n'est pas garanti

- Aucune recherche ne prouve qu'il n'existe pas mieux ailleurs : on teste une
  centaine de configs sur des milliards possibles.
- Si l'écart vient d'autre chose que des paramètres GATK (hypothèse 4.4), le
  tuning ne suffira pas et il faudra le dire plutôt que de multiplier les
  essais.
- Les validateurs n'ont pas tous le même matériel. Le rapport donne le temps
  d'exécution maximal de chaque config : une config proche du timeout chez toi
  peut échouer chez eux.

## 7. Commandes

```bash
# Récupérer les 8 derniers vrais rounds chr12/chr13 (+ score top 1 / top 20)
python scripts/optimize_config.py fetch-rounds --count 8

# Vérifier qu'une config est bien celle que les validateurs exécutent
python scripts/optimize_config.py check --tool gatk --conf configs/gatk.conf

# Mesurer des configs sur ces rounds
python scripts/optimize_config.py baseline --manifest datasets/rounds/manifest.json \
    --holdout 2 --run-dir optimize_runs/real --conf configs/candidates/user/soft_thr65.conf

# Recherche (reprend là où elle s'est arrêtée)
python scripts/optimize_config.py search --manifest datasets/rounds/manifest.json \
    --holdout 2 --run-dir optimize_runs/real --trials 40 \
    --conf configs/candidates/user/soft_thr65.conf \
    --params min_pruning,standard_min_confidence_threshold_for_calling,heterozygosity,min_dangling_branch_length,min_mapping_quality_score,min_base_quality_score,base_quality_score_threshold,pcr_indel_model,phred_scaled_global_read_mismapping_rate,pair_hmm_gap_continuation_penalty,pruning_lod_threshold

# Classement final + écriture de la meilleure config
python scripts/optimize_config.py report --manifest datasets/rounds/manifest.json \
    --holdout 2 --run-dir optimize_runs/real
```

## 8. Notes d'exécution

- La session cloud se met en veille dès qu'aucune tâche n'y tourne, ce qui tue
  les processus mais garde les fichiers. La recherche a été coupée ainsi une
  fois (perte d'environ 3 h). Depuis, une tâche de garde de 2 h lance
  `scripts/resume_search.sh` toutes les 5 minutes et elle est relancée à son
  expiration. Le script redémarre Docker et la recherche si besoin, sans rien
  recalculer de ce qui est fini.
- Le script a montré deux bugs pendant les premiers runs, corrigés depuis : une
  panne Docker était enregistrée comme un score de 0, et les identifiants de
  round contenant `:` cassaient les montages Docker. Les résultats affectés ont
  été effacés et recalculés.
