varBscore: calculate varBscore using vcf input
==============================================

## Usage

---

### basic usage

```bash
snpScore-mp2 \
    -p '{"mutant":["mutan"],"wild":["wild"],"mutant_parent":[],"wild_parent":[],"background":[],"qtlseqr_window":1000000,"qtlseqr":"qtlseqr","ed":"ed","min_depth":5,"snp_number_window":20,"snp_number_step":5,"ref_freq":0.3,"p_ref_freq":0.3,"background_ref_freq":0.3,"qtlseqr_min_depth":5,"qtlseqr_ref_freq":0.3,"pop_stru":"RIL","snp_density_window":100000,"snp_density_step":100000}' \
    --vcf_dir vcf.dir1 \
    --vcf_dir vcf.dir2 \
    --chr_size chr.size \
    --thread 4 \
    --plant \
    -o out.dir \
    --snpeff_db snpeff.db.name \
    --snpeff_cfg snpeff.config.path
```
