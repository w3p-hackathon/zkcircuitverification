use std::string;

struct GeneticMarker {
    // rsid: Field,           // Hash of rsid string
    chromosome: u8,        // 1-22, 23=X, 24=Y, 25=MT
    position: u32,         // Genomic position
    allele1: u8,          // Encoded alleles: A=1, T=2, G=3, C=4, 0=0 (missing)
    allele2: u8,
}

struct QualityMetrics {
    call_rate: Field,          // As percentage * 100 (e.g., 9500 = 95.00%)
    heterozygosity_rate: Field, // As percentage * 100
    ti_tv_ratio: Field,        // Multiplied by 100 (e.g., 210 = 2.10)
    total_snps: u32,
}
// MAX_SNPS 0xab435
global MAX_SNPS: u32 = 1000;
fn calculate_rate(dna: [GeneticMarker; MAX_SNPS]) -> [u32; 3]{
    let mut total_valid: u32 = 0;
    let mut missing_count: u32 = 0;
    let mut heteroz: u32 = 0;
    let mut tis: u32 = 0;
    let mut tvs: u32 = 0;
    for i in 0 .. dna.len() {
        if (dna[i].allele1 == 0) | (dna[i].allele2 == 0) {
            missing_count += 1;
        }else if dna[i].allele1 != dna[i].allele2 {
            heteroz += 1;
        }
        if (dna[i].allele1 == 1) & (dna[i].allele2 == 3) | (dna[i].allele1 == 2) & (dna[i].allele2 == 4) {
            tis += 1;
        }else if (dna[i].allele1 != 0) | (dna[i].allele1 != 0) {
            tvs += 1;
        }
        total_valid += 1; 
    }
    let mut call_rate = 0;
    if total_valid > 0 {
        call_rate = ((total_valid-  missing_count) * 1000 ) / total_valid;
    }

    let mut heterozygosity_rate: u32 = 0;
    if (total_valid -missing_count) > 0 {
        heterozygosity_rate = (heteroz  * 100) / (total_valid - missing_count);
    }

    let mut ti_tv_rate: u32 = 0;

    if tvs > 0 {
        ti_tv_rate = (tis * 100) / tvs;
    }

    [call_rate, heterozygosity_rate, ti_tv_rate]

}


fn main(dna: [GeneticMarker; MAX_SNPS], challenge_hash: pub Field, min_heterozygosity_rate: pub u32, ti_tv_ratio: pub u32, min_call_rate: u32){
    let rates = calculate_rate(dna);
    assert( rates[0] == min_call_rate, "call rates don't match");
    assert( rates[1] == min_heterozygosity_rate, "heterozygosity doesn't match");
    assert( rates[2] <= ti_tv_ratio, "ti_tv doesn't match");


}

// #[test]
// fn test_main() {
//     main(1, 2);

    // Uncomment to make test fail
    // main(1, 1);
// }
