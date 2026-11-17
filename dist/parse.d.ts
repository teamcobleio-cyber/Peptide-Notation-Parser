/**
 * A parser for peptide notation as it is actually written: one-letter runs,
 * three-letter hyphenated chains, stereo prefixes, terminal caps, non-standard
 * residues, side-chain modifications, cyclo[...] rings and disulfide notes.
 *
 * The rule that shapes everything here: nothing is dropped. A token the parser
 * cannot place is returned in `unrecognised`, and an ambiguity it resolved by
 * convention is returned in `warnings`. A parser that quietly skipped what it
 * did not understand would report part of a molecule as the whole of it.
 */
export type Stereo = "L" | "D" | null;
export interface ParsedResidue {
    /** 1-based position across the whole result. */
    position: number;
    /** 0-based index into `chains`. */
    chain: number;
    /** Canonical abbreviation ("Ala", "Aib", "2-Nal"), or "Xaa" for a placeholder. */
    abbreviation: string;
    /** One-letter code when the residue has one, otherwise null. */
    oneLetter: string | null;
    name: string;
    /** "L" or "D" for a chiral residue, null for an achiral one (Gly, Aib). */
    stereo: Stereo;
    /** True when the notation stated the configuration (a D- prefix, a lower-case letter). */
    stereoExplicit: boolean;
    /** True for the twenty proteinogenic amino acids, whatever their configuration. */
    standard: boolean;
    /** Backbone or side-chain modifications, e.g. "side chain: acyl", "N-methyl". */
    modifications: string[];
    /** The text this residue was read from. */
    source: string;
}
export type NTerminusType = "unspecified" | "free" | "acetyl" | "acyl";
export type CTerminusType = "unspecified" | "acid" | "amide" | "substituted-amide" | "group";
export interface Terminus<T extends string> {
    type: T;
    /** The text as written, e.g. "Ac", "(3E)-hex-3-enoyl", "NH2". */
    text: string | null;
}
export interface Chain {
    /** "Alpha chain" and similar, when the input labels its chains. */
    label: string | null;
    /** 1-based positions of the first and last residue in this chain. */
    start: number;
    end: number;
    nTerminus: Terminus<NTerminusType>;
    cTerminus: Terminus<CTerminusType>;
}
export interface Cyclisation {
    kind: "cyclo" | "disulfide";
    /** 1-based positions. */
    start: number;
    end: number;
    source: string;
}
export interface Unrecognised {
    token: string;
    /** 0-based chain index the token was found in. */
    chain: number;
}
export interface ParseResult {
    input: string;
    residues: ParsedResidue[];
    residueCount: number;
    chains: Chain[];
    /** The first chain's N-terminus, for the common single-chain case. */
    nTerminus: Terminus<NTerminusType>;
    /** The last chain's C-terminus, for the common single-chain case. */
    cTerminus: Terminus<CTerminusType>;
    cyclisations: Cyclisation[];
    unrecognised: Unrecognised[];
    warnings: string[];
}
export interface ParseOptions {
    /**
     * Read a lower-case one-letter code as the D-residue (the convention used for
     * retro-inverso and all-D peptides). Default true.
     */
    lowercaseIsD?: boolean;
    /**
     * What a one-letter "U" means. IUPAC assigns it to selenocysteine (the
     * default); some registries use it for pyroglutamate.
     */
    u?: "selenocysteine" | "pyroglutamate";
}
/** Parse peptide notation into residues, termini, rings and anything it could not place. */
export declare function parsePeptide(input: string, options?: ParseOptions): ParseResult;
/**
 * The plain one-letter sequence, or null.
 *
 * Returns a string only when every residue is a standard L (or achiral)
 * residue with no modification, nothing went unrecognised, and the termini are
 * free or unstated. Pass `{ ignoreTermini: true }` to accept capped termini.
 * A multi-chain result returns one sequence per chain, joined by "/".
 */
export declare function toOneLetter(result: ParseResult, options?: {
    ignoreTermini?: boolean;
}): string | null;
export interface CompositionSummary {
    total: number;
    proteinogenic: number;
    nonStandard: number;
    dResidues: number;
    modified: number;
    placeholders: number;
    unrecognised: number;
    /** Counts keyed by canonical abbreviation, configuration ignored. */
    byResidue: Record<string, number>;
}
export declare function composition(result: ParseResult): CompositionSummary;
/** Re-write a parse as normalised three-letter notation. */
export declare function toThreeLetter(result: ParseResult): string;
