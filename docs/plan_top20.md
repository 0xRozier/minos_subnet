# Minos SN107 — Observations et plan pour le top 20

*Rédigé le 2026-10-11. Les chiffres viennent de runs réels : GATK 4.5.0.0 →
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

| Config | Moyenne | Écart moyen vs top 1 | Pire écart |
|---|---|---|---|
| Défauts GATK | 75,99 | −6,25 | −12,38 |
| `soft_thr65` reconstruite (`journal_best`) | **80,42** | **−1,82** | −2,25 |
| idem + `contamination_fraction_to_filter=0.05` | ≈80,35 | −1,89 | −2,38 |
| `soft_thr65`, `soft_base`, `ei_best`, `gatk_defaults_plus`, `ei_balanced` (fichiers exacts) | *en cours* | | |

Détail de `journal_best` :

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
| 1 | Mesurer tes configs exactes, en version live, sur les 8 rounds | 🔄 en cours | Choisir la meilleure comme base |
| 2 | Recherche bayésienne sur les 6 rounds de tuning, partant des meilleures configs, restreinte aux leviers (§3.4) + paramètres peu explorés (§4.3), ~40 essais (~15–20 h sur 4 cœurs) | à venir | Gain ≥ +0,3 point sur les 6 rounds **et** confirmé sur les 2 rounds mis de côté |
| 3 | Ablation « vers les défauts » de la meilleure config (hypothèse 4.2), et ajustement des FP `indel_8bp` / `snp_het` | à venir | Idem |
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
