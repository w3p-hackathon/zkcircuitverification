// {
//     1 "European": {
//         "rs4988235": "T",
//         "rs12913832": "G",
//         "rs1426654": "A",
//         "rs16891982": "G"
//     },
//     2"Sub-Saharan African": {
//         "rs2814778": "C",
//         "rs334": "T",
//         "rs1426654": "G"
//     },
//     3 "South Asian": {
//         "rs1042602": "A",
//         "rs1426654": "A",
//         "rs16891982": "G",
//         "rs11190870": "T"
//     },
//     4 "East Asian": {
//         "rs3827760": "G",
//         "rs1229984": "A",
//         "rs17822931": "G"
//     },
//     5 "Amerindian": {
//         "rs174537": "G",
//         "rs3827760": "G",
//         "rs1426654": "A"
//     },
//     6 "Papuan": {
//         "rs3827760": "A",
//         "rs1426654": "G",
//         "rs2814778": "T",
//         "rs12913832": "A"
//     },
//     7 "West Asian/North African": {
//         "rs1426654": "A",
//         "rs16891982": "G",
//         "rs12913832": "G",
//         "rs2814778": "C"
//     }
// }

struct GeneticMarker {
    rsid: Field,
    chromosome: u8,
    position: u32,
    allele1: u8,
    allele2: u8,
}
global MAX_SNPS: u32 = 1000;

fn main(dna: [GeneticMarker; MAX_SNPS], ancestor: pub u32) {

    let ancestries = [
        [4988235,12913832,1426654], 
        [2814778,334,1426654], 
        [1042602,1426654,16891982], 
        [3827760,1229984,17822931], 
        [174537,3827760,1426654], 
        [3827760,1426654,2814778], 
        [1426654,16891982,12913832]
    ];
    let mut ratings = [0, 0, 0, 0,0 ,0, 0];
    let mut winner = 0;
    for i in 0 .. dna.len() {
        for j in 0 .. ancestries.len() {
            if (dna[i].rsid == ancestries[j][0]) | (dna[i].rsid == ancestries[j][1]) | (dna[i].rsid == ancestries[j][2]){
                ratings[j] += 1;
                if ratings[j] > ratings[winner] {
                    winner = j;
                }
            }
        }
    }

    assert(ancestor == winner);

}