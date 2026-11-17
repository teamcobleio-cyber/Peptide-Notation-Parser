/**
 * Residue vocabulary: the twenty proteinogenic amino acids, the two genetically
 * encoded extras, and the non-standard residues that turn up in research
 * peptide notation.
 *
 * Every entry is an identity, not a mass. This library deliberately carries no
 * mass table for non-standard residues: a mass needs a stated source and a
 * stated form (free acid, residue, salt), and a guessed one is worse than none.
 */
const P = (abbreviation, oneLetter, name, chiral = true) => ({
    abbreviation,
    oneLetter,
    name,
    proteinogenic: true,
    chiral,
});
export const PROTEINOGENIC = [
    P("Ala", "A", "Alanine"),
    P("Arg", "R", "Arginine"),
    P("Asn", "N", "Asparagine"),
    P("Asp", "D", "Aspartic acid"),
    P("Cys", "C", "Cysteine"),
    P("Gln", "Q", "Glutamine"),
    P("Glu", "E", "Glutamic acid"),
    P("Gly", "G", "Glycine", false),
    P("His", "H", "Histidine"),
    P("Ile", "I", "Isoleucine"),
    P("Leu", "L", "Leucine"),
    P("Lys", "K", "Lysine"),
    P("Met", "M", "Methionine"),
    P("Phe", "F", "Phenylalanine"),
    P("Pro", "P", "Proline"),
    P("Ser", "S", "Serine"),
    P("Thr", "T", "Threonine"),
    P("Trp", "W", "Tryptophan"),
    P("Tyr", "Y", "Tyrosine"),
    P("Val", "V", "Valine"),
];
const N = (abbreviation, name, options = {}) => ({
    abbreviation,
    oneLetter: options.oneLetter ?? null,
    name,
    proteinogenic: false,
    chiral: options.chiral ?? true,
    parent: options.parent,
    aliases: options.aliases,
});
/** Non-standard residues. Extend at runtime with `registerResidue`. */
export const NON_STANDARD = [
    N("Sec", "Selenocysteine", { oneLetter: "U", aliases: ["SEC"] }),
    N("Pyl", "Pyrrolysine", { oneLetter: "O", aliases: ["PYL"] }),
    N("Aib", "2-Aminoisobutyric acid (alpha-methylalanine)", { chiral: false, aliases: ["AIB", "2-methylalanine", "alpha-Me-Ala", "alphaMeAla"] }),
    N("Nle", "Norleucine (2-aminohexanoic acid)", { aliases: ["NLE"] }),
    N("Nva", "Norvaline", { aliases: ["NVA"] }),
    N("Orn", "Ornithine", { aliases: ["ORN"] }),
    N("Cit", "Citrulline", { aliases: ["CIT"] }),
    N("Abu", "2-Aminobutyric acid", { aliases: ["ABU"] }),
    N("Dab", "2,4-Diaminobutyric acid", { aliases: ["DAB"] }),
    N("Dap", "2,3-Diaminopropionic acid", { aliases: ["DAP", "Dpr"] }),
    N("Sar", "Sarcosine (N-methylglycine)", { chiral: false, aliases: ["SAR", "MeGly", "N-Me-Gly"] }),
    N("Hyp", "4-Hydroxyproline", { parent: "Pro", aliases: ["HYP", "4-Hyp"] }),
    N("Hcy", "Homocysteine", { aliases: ["HCY"] }),
    N("Ac4c", "1-Aminocyclobutane-1-carboxylic acid", { chiral: false, aliases: ["AC4C"] }),
    N("Ac5c", "1-Aminocyclopentane-1-carboxylic acid", { chiral: false, aliases: ["AC5C", "Cle"] }),
    N("Ac6c", "1-Aminocyclohexane-1-carboxylic acid", { chiral: false, aliases: ["AC6C"] }),
    N("2-Nal", "3-(2-Naphthyl)alanine", { aliases: ["2Nal", "Nal(2)", "NAL", "2-Nal"] }),
    N("1-Nal", "3-(1-Naphthyl)alanine", { aliases: ["1Nal", "Nal(1)"] }),
    N("Dmt", "2',6'-Dimethyltyrosine", { parent: "Tyr", aliases: ["DMT"] }),
    N("pGlu", "Pyroglutamic acid (5-oxoproline)", { aliases: ["Pyr", "PCA", "pE", "Glp", "<Glu"] }),
    N("gGlu", "Glutamic acid linked through its gamma-carboxyl", {
        parent: "Glu",
        aliases: ["gamma-Glu", "γ-Glu", "gammaGlu", "g-Glu", "isoGlu"],
    }),
    N("bAla", "beta-Alanine", { chiral: false, aliases: ["beta-Ala", "β-Ala", "betaAla", "BAL"] }),
    N("Ahx", "6-Aminohexanoic acid", { chiral: false, aliases: ["AHX", "eAhx", "epsilon-Ahx"] }),
    N("Tic", "1,2,3,4-Tetrahydroisoquinoline-3-carboxylic acid", { aliases: ["TIC"] }),
    N("Cha", "Cyclohexylalanine", { aliases: ["CHA"] }),
    N("Phg", "Phenylglycine", { aliases: ["PHG"] }),
    N("Bip", "4,4'-Biphenylalanine", { aliases: ["BIP"] }),
    N("Me-Leu", "alpha-Methylleucine (2-methylleucine)", {
        parent: "Leu",
        aliases: ["alpha-Me-Leu", "alphaMeLeu", "aMeLeu", "2ML", "alpha-methyl-Leu"],
    }),
    N("Me-Phe", "alpha-Methylphenylalanine", { parent: "Phe", aliases: ["alpha-Me-Phe", "alphaMePhe", "aMePhe"] }),
    N("Me-Lys", "alpha-Methyllysine", { parent: "Lys", aliases: ["alpha-Me-Lys", "alphaMeLys", "aMeLys"] }),
    N("2-Me-Trp", "2-Methyltryptophan (methyl on indole C2)", {
        parent: "Trp",
        aliases: ["2-methyl-Trp", "2MeTrp", "2-Me-Trp"],
    }),
];
const registry = new Map();
const oneLetterIndex = new Map();
function key(text) {
    return text.toLowerCase();
}
/** Add a residue (or replace one with the same abbreviation). */
export function registerResidue(definition) {
    registry.set(key(definition.abbreviation), definition);
    for (const alias of definition.aliases ?? [])
        registry.set(key(alias), definition);
    if (definition.oneLetter)
        oneLetterIndex.set(definition.oneLetter, definition);
}
for (const definition of [...PROTEINOGENIC, ...NON_STANDARD])
    registerResidue(definition);
/** Look up a residue by abbreviation or alias, case-insensitively. */
export function lookupResidue(text) {
    return registry.get(key(text));
}
/** Look up a one-letter code (upper case). */
export function lookupOneLetter(letter) {
    return oneLetterIndex.get(letter.toUpperCase());
}
