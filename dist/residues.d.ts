/**
 * Residue vocabulary: the twenty proteinogenic amino acids, the two genetically
 * encoded extras, and the non-standard residues that turn up in research
 * peptide notation.
 *
 * Every entry is an identity, not a mass. This library deliberately carries no
 * mass table for non-standard residues: a mass needs a stated source and a
 * stated form (free acid, residue, salt), and a guessed one is worse than none.
 */
export interface ResidueDefinition {
    /** Canonical abbreviation used in output, e.g. "Ala", "Aib", "2-Nal". */
    abbreviation: string;
    /** One-letter code when one exists and is unambiguous, otherwise null. */
    oneLetter: string | null;
    /** Full name. */
    name: string;
    /** True for the twenty proteinogenic amino acids only. */
    proteinogenic: boolean;
    /** False for residues with no alpha stereocentre (Gly, Aib, Ac4c, Sar). */
    chiral: boolean;
    /**
     * For a derivative of a proteinogenic residue, the parent's abbreviation
     * (Nle has none; alpha-methyl-Leu has "Leu").
     */
    parent?: string;
    /** Other spellings that resolve to this entry, matched case-insensitively. */
    aliases?: string[];
}
export declare const PROTEINOGENIC: ResidueDefinition[];
/** Non-standard residues. Extend at runtime with `registerResidue`. */
export declare const NON_STANDARD: ResidueDefinition[];
/** Add a residue (or replace one with the same abbreviation). */
export declare function registerResidue(definition: ResidueDefinition): void;
/** Look up a residue by abbreviation or alias, case-insensitively. */
export declare function lookupResidue(text: string): ResidueDefinition | undefined;
/** Look up a one-letter code (upper case). */
export declare function lookupOneLetter(letter: string): ResidueDefinition | undefined;
