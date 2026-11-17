import { lookupOneLetter, lookupResidue, type ResidueDefinition } from "./residues.js";

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

type Draft = Omit<ParsedResidue, "position" | "chain">;

interface BodyRead {
  residues: Draft[];
  /** Cyclisations with 0-based offsets relative to this body. */
  cycles: { kind: "cyclo"; start: number; end: number; source: string }[];
  unrecognised: string[];
  warnings: string[];
  nCap?: Terminus<NTerminusType>;
  cCap?: Terminus<CTerminusType>;
}

/** PDB component codes for D-residues, as deposited structures write them. */
const PDB_D_RESIDUES: Record<string, string> = {
  DAL: "Ala", DAR: "Arg", DSG: "Asn", DAS: "Asp", DCY: "Cys", DGN: "Gln",
  DGL: "Glu", DHI: "His", DIL: "Ile", DLE: "Leu", DLY: "Lys", MED: "Met",
  DPN: "Phe", DPR: "Pro", DSN: "Ser", DTH: "Thr", DTR: "Trp", DTY: "Tyr", DVA: "Val",
};

const N_CAP = /(?:^(?:Ac|Boc|Fmoc|For|Bz|Pal|Myr)$|o?yl$)/i;

const MODIFIER_NAMES: Record<string, string> = {
  alphame: "alpha-methyl",
  ame: "alpha-methyl",
  nme: "N-methyl",
  me: "N-methyl",
};

function splitTopLevel(text: string, separator: string): string[] {
  const parts: string[] = [];
  let depth = 0;
  let current = "";
  for (const char of text) {
    if (char === "(" || char === "[" || char === "{") depth += 1;
    if (char === ")" || char === "]" || char === "}") depth = Math.max(0, depth - 1);
    if (char === separator && depth === 0) {
      parts.push(current);
      current = "";
    } else {
      current += char;
    }
  }
  parts.push(current);
  return parts.filter((part) => part.length > 0);
}

function fromDefinition(definition: ResidueDefinition, source: string, stereo?: Stereo): Draft {
  return {
    abbreviation: definition.abbreviation,
    oneLetter: definition.oneLetter,
    name: definition.name,
    stereo: definition.chiral ? (stereo ?? "L") : null,
    stereoExplicit: Boolean(stereo) && definition.chiral,
    standard: definition.proteinogenic,
    modifications: [],
    source,
  };
}

function placeholder(source: string): Draft {
  return {
    abbreviation: "Xaa",
    oneLetter: "X",
    name: "Unspecified residue",
    stereo: null,
    stereoExplicit: false,
    standard: false,
    modifications: [],
    source,
  };
}

/** A single residue token: "Lys", "Aib", "(alpha-Me-Leu)", "Lys(acyl)", "(NMe)Asn", "DPN". */
function readResidueToken(token: string): Draft | null {
  const whole = /^\((.+)\)$/.exec(token);
  if (whole) return readResidueToken(whole[1]);

  const prefixed = /^\(([^()]+)\)(.+)$/.exec(token);
  if (prefixed) {
    const base = readResidueToken(prefixed[2]);
    if (!base) return null;
    const modifier = MODIFIER_NAMES[prefixed[1].toLowerCase().replace(/[^a-z]/g, "")] ?? prefixed[1];
    return { ...base, modifications: [...base.modifications, modifier], source: token };
  }

  const sideChain = /^(.+?)\((.+)\)$/.exec(token);
  if (sideChain) {
    const base = readResidueToken(sideChain[1]);
    if (!base) return null;
    return { ...base, modifications: [...base.modifications, `side chain: ${sideChain[2]}`], source: token };
  }

  const pdbD = PDB_D_RESIDUES[token.toUpperCase()];
  if (pdbD && token === token.toUpperCase()) {
    const definition = lookupResidue(pdbD);
    if (definition) return fromDefinition(definition, token, "D");
  }

  const definition = lookupResidue(token);
  return definition ? fromDefinition(definition, token) : null;
}

/** True when a bare token could be a run of one-letter codes. */
function looksLikeOneLetterRun(token: string): boolean {
  if (!/^[A-Za-z]+(?:\([A-Za-z0-9]+\)[A-Za-z]*)*$/.test(token)) return false;
  const letters = token.replace(/\([^)]*\)/g, "");
  if (letters.length === 0) return true;
  if (letters === letters.toLowerCase()) return true;
  // Mixed case is allowed only as isolated lower-case D-residues ("YaDAIF"),
  // never as a word ("Hexanoyl").
  return !/[a-z]{2}/.test(letters);
}

function readOneLetterRun(token: string, options: Required<ParseOptions>): BodyRead | null {
  if (!looksLikeOneLetterRun(token)) return null;
  const read: BodyRead = { residues: [], cycles: [], unrecognised: [], warnings: [] };
  let i = 0;
  while (i < token.length) {
    const char = token[i];
    if (char === "(") {
      const close = token.indexOf(")", i);
      const code = token.slice(i + 1, close);
      i = close + 1;
      if (/^NH2$/i.test(code) && i === token.length) {
        read.cCap = { type: "amide", text: code };
        continue;
      }
      if (/^ACE$/i.test(code) && read.residues.length === 0) {
        read.nCap = { type: "acetyl", text: code };
        continue;
      }
      const residue = readResidueToken(code);
      if (residue) read.residues.push(residue);
      else read.unrecognised.push(`(${code})`);
      continue;
    }
    i += 1;
    const upper = char.toUpperCase();
    if (upper === "X" || upper === "B" || upper === "Z") {
      read.residues.push(placeholder(char));
      read.warnings.push(`"${char}" is a placeholder or ambiguity code, not a defined residue`);
      continue;
    }
    const definition = lookupOneLetter(upper);
    if (!definition) return null;
    const isLower = char !== upper;
    if (upper === "U") {
      if (options.u === "pyroglutamate") {
        const pGlu = lookupResidue("pGlu");
        if (pGlu) {
          read.residues.push(fromDefinition(pGlu, char));
          continue;
        }
      }
      read.warnings.push(
        `"U" was read as selenocysteine; some sources use it for pyroglutamate (pass { u: "pyroglutamate" })`,
      );
    }
    read.residues.push(fromDefinition(definition, char, isLower && options.lowercaseIsD ? "D" : undefined));
  }
  return read;
}

/** Read one residue (or a cyclo[...] group or a run) starting at `parts[index]`. */
function readAt(
  parts: string[],
  index: number,
  options: Required<ParseOptions>,
): { read: BodyRead; consumed: number } | null {
  const part = parts[index];
  if (part === undefined) return null;

  const cyclo = /^cyclo\s*[[(](.*)[\])]$/i.exec(part);
  if (cyclo) {
    const inner = readBody(splitTopLevel(cyclo[1], "-"), options);
    if (inner.residues.length === 0) return null;
    inner.cycles.push({ kind: "cyclo", start: 0, end: inner.residues.length - 1, source: part });
    return { read: inner, consumed: 1 };
  }

  if ((part === "D" || part === "L") && index + 1 < parts.length) {
    const next = readSingle(parts, index + 1);
    if (next) {
      const residue = next.residue;
      if (residue.stereo !== null) {
        residue.stereo = part;
        residue.stereoExplicit = true;
      }
      residue.source = `${part}-${residue.source}`;
      return { read: { residues: [residue], cycles: [], unrecognised: [], warnings: [] }, consumed: next.consumed + 1 };
    }
  }

  const single = readSingle(parts, index);
  if (single) {
    return {
      read: { residues: [single.residue], cycles: [], unrecognised: [], warnings: [] },
      consumed: single.consumed,
    };
  }

  const run = readOneLetterRun(part, options);
  if (run && run.residues.length > 0) return { read: run, consumed: 1 };
  return null;
}

/** A named residue, trying multi-part names ("2-methyl-Trp", "gamma-Glu") first. */
function readSingle(parts: string[], index: number): { residue: Draft; consumed: number } | null {
  for (let span = Math.min(4, parts.length - index); span >= 1; span -= 1) {
    const text = parts.slice(index, index + span).join("-");
    // A bare upper-case word is a one-letter run ("GHK", "AED"), not an
    // abbreviation, unless it is written in title case.
    if (span === 1 && /^[A-Z]{3,}$/.test(text) && !PDB_D_RESIDUES[text]) continue;
    const residue = readResidueToken(text);
    if (residue) return { residue, consumed: span };
  }
  return null;
}

function readBody(parts: string[], options: Required<ParseOptions>): BodyRead {
  const body: BodyRead = { residues: [], cycles: [], unrecognised: [], warnings: [] };
  let index = 0;
  while (index < parts.length) {
    const next = readAt(parts, index, options);
    if (!next) {
      body.unrecognised.push(parts[index]);
      index += 1;
      continue;
    }
    const offset = body.residues.length;
    body.residues.push(...next.read.residues);
    for (const cycle of next.read.cycles) {
      body.cycles.push({ ...cycle, start: cycle.start + offset, end: cycle.end + offset });
    }
    body.unrecognised.push(...next.read.unrecognised);
    body.warnings.push(...next.read.warnings);
    if (next.read.nCap && offset === 0 && !body.nCap) body.nCap = next.read.nCap;
    if (next.read.cCap) body.cCap = next.read.cCap;
    index += next.consumed;
  }
  return body;
}

interface ChainRead extends BodyRead {
  nTerminus: Terminus<NTerminusType>;
  cTerminus: Terminus<CTerminusType>;
}

function readChain(text: string, options: Required<ParseOptions>): ChainRead {
  const parts = splitTopLevel(text.replace(/\s+/g, ""), "-");
  let nTerminus: Terminus<NTerminusType> = { type: "unspecified", text: null };
  let cTerminus: Terminus<CTerminusType> = { type: "unspecified", text: null };
  let start = 0;
  let end = parts.length;

  if (parts.length > 1 && parts[0] === "H") {
    nTerminus = { type: "free", text: "H" };
    start = 1;
  } else if (!readSingle(parts, 0) && parts[0] !== "D" && parts[0] !== "L") {
    // Only look for a cap when the chain does not open on a named residue:
    // "His-D-2-methyl-Trp" ends in "yl" after three parts and is not a cap.
    for (let k = 0; k < Math.min(5, parts.length - 1); k += 1) {
      const candidate = parts.slice(0, k + 1).join("-");
      if (N_CAP.test(candidate) && readAt(parts, k + 1, options)) {
        nTerminus = /^(?:Ac|acetyl)$/i.test(candidate)
          ? { type: "acetyl", text: candidate }
          : { type: "acyl", text: candidate };
        start = k + 1;
        break;
      }
    }
  }

  const last = parts[end - 1];
  if (end - start > 1 && /^OH$/i.test(last)) {
    cTerminus = { type: "acid", text: last };
    end -= 1;
  } else if (end - start > 1 && /^NH2$/i.test(last)) {
    cTerminus = { type: "amide", text: last };
    end -= 1;
  } else if (end - start > 1 && /^NH(?:Et|Me|Pr|iPr|Bu)$/i.test(last)) {
    cTerminus = { type: "substituted-amide", text: last };
    end -= 1;
  } else {
    const nh = parts.indexOf("NH", start + 1);
    if (nh > start) {
      cTerminus = { type: "substituted-amide", text: parts.slice(nh).join("-") };
      end = nh;
    } else if (end - start > 1 && /^\[.+\]$/.test(last)) {
      cTerminus = { type: "group", text: last };
      end -= 1;
    }
  }

  const body = readBody(parts.slice(start, end), options);
  if (body.nCap && nTerminus.type === "unspecified") nTerminus = body.nCap;
  if (body.cCap && cTerminus.type === "unspecified") cTerminus = body.cCap;
  return { ...body, nTerminus, cTerminus };
}

const CHAIN_LABEL = /([A-Za-z][A-Za-z0-9 ]*?\bchain)\s*(?:\(\s*\d+\s*residues?\s*\))?\s*:\s*/gi;

const DISULFIDE_NOTES: RegExp[] = [
  /\(\s*(?:cyclic\s+)?(\d+)\s*-\s*(\d+)\s+disulfide\s*\)/gi,
  /\(\s*disulfide\s*:?\s*(\d+)\s*-\s*(\d+)\s*\)/gi,
];

/** Parse peptide notation into residues, termini, rings and anything it could not place. */
export function parsePeptide(input: string, options: ParseOptions = {}): ParseResult {
  const settings: Required<ParseOptions> = {
    lowercaseIsD: options.lowercaseIsD ?? true,
    u: options.u ?? "selenocysteine",
  };
  let text = input.trim();

  const disulfides: { start: number; end: number; source: string }[] = [];
  for (const pattern of DISULFIDE_NOTES) {
    text = text.replace(pattern, (match, a: string, b: string) => {
      disulfides.push({ start: Number(a), end: Number(b), source: match.trim() });
      return "";
    }).trim();
  }

  const segments: { label: string | null; text: string }[] = [];
  const labels = [...text.matchAll(CHAIN_LABEL)];
  if (labels.length > 0) {
    labels.forEach((match, i) => {
      const from = (match.index ?? 0) + match[0].length;
      const to = i + 1 < labels.length ? labels[i + 1].index ?? text.length : text.length;
      segments.push({ label: match[1].trim(), text: text.slice(from, to).replace(/[\s|;,]+$/, "") });
    });
  } else {
    segments.push({ label: null, text });
  }

  const result: ParseResult = {
    input,
    residues: [],
    residueCount: 0,
    chains: [],
    nTerminus: { type: "unspecified", text: null },
    cTerminus: { type: "unspecified", text: null },
    cyclisations: [],
    unrecognised: [],
    warnings: [],
  };

  segments.forEach((segment, chainIndex) => {
    const read = readChain(segment.text, settings);
    const offset = result.residues.length;
    read.residues.forEach((residue, i) => {
      result.residues.push({ ...residue, position: offset + i + 1, chain: chainIndex });
    });
    for (const cycle of read.cycles) {
      result.cyclisations.push({ kind: "cyclo", start: offset + cycle.start + 1, end: offset + cycle.end + 1, source: cycle.source });
    }
    for (const token of read.unrecognised) result.unrecognised.push({ token, chain: chainIndex });
    result.warnings.push(...read.warnings);
    result.chains.push({
      label: segment.label,
      start: offset + 1,
      end: offset + read.residues.length,
      nTerminus: read.nTerminus,
      cTerminus: read.cTerminus,
    });
  });

  for (const bond of disulfides) {
    result.cyclisations.push({ kind: "disulfide", ...bond });
    for (const position of [bond.start, bond.end]) {
      const residue = result.residues[position - 1];
      if (!residue || residue.abbreviation !== "Cys") {
        result.warnings.push(`disulfide note names position ${position}, which is not a cysteine`);
      }
    }
  }

  result.residueCount = result.residues.length;
  result.nTerminus = result.chains[0]?.nTerminus ?? result.nTerminus;
  result.cTerminus = result.chains[result.chains.length - 1]?.cTerminus ?? result.cTerminus;
  return result;
}

/**
 * The plain one-letter sequence, or null.
 *
 * Returns a string only when every residue is a standard L (or achiral)
 * residue with no modification, nothing went unrecognised, and the termini are
 * free or unstated. Pass `{ ignoreTermini: true }` to accept capped termini.
 * A multi-chain result returns one sequence per chain, joined by "/".
 */
export function toOneLetter(result: ParseResult, options: { ignoreTermini?: boolean } = {}): string | null {
  if (result.unrecognised.length > 0 || result.residues.length === 0) return null;
  if (!options.ignoreTermini) {
    for (const chain of result.chains) {
      if (!["unspecified", "free"].includes(chain.nTerminus.type)) return null;
      if (!["unspecified", "acid"].includes(chain.cTerminus.type)) return null;
    }
  }
  const chains: string[] = [];
  for (const chain of result.chains) {
    let out = "";
    for (let p = chain.start; p <= chain.end; p += 1) {
      const residue = result.residues[p - 1];
      if (!residue.standard || !residue.oneLetter) return null;
      if (residue.stereo === "D" || residue.modifications.length > 0) return null;
      out += residue.oneLetter;
    }
    chains.push(out);
  }
  return chains.join("/");
}

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

export function composition(result: ParseResult): CompositionSummary {
  const summary: CompositionSummary = {
    total: result.residueCount,
    proteinogenic: 0,
    nonStandard: 0,
    dResidues: 0,
    modified: 0,
    placeholders: 0,
    unrecognised: result.unrecognised.length,
    byResidue: {},
  };
  for (const residue of result.residues) {
    if (residue.standard) summary.proteinogenic += 1;
    else summary.nonStandard += 1;
    if (residue.stereo === "D") summary.dResidues += 1;
    if (residue.modifications.length > 0) summary.modified += 1;
    if (residue.abbreviation === "Xaa") summary.placeholders += 1;
    summary.byResidue[residue.abbreviation] = (summary.byResidue[residue.abbreviation] ?? 0) + 1;
  }
  return summary;
}

/** Re-write a parse as normalised three-letter notation. */
export function toThreeLetter(result: ParseResult): string {
  return result.chains
    .map((chain) => {
      const pieces: string[] = [];
      for (let p = chain.start; p <= chain.end; p += 1) {
        const r = result.residues[p - 1];
        const prefix = r.stereo === "D" ? "D-" : "";
        const sideChain = r.modifications.filter((m) => m.startsWith("side chain: ")).map((m) => m.slice(12));
        const backbone = r.modifications.filter((m) => !m.startsWith("side chain: "));
        let text = `${prefix}${backbone.length ? `(${backbone.join(",")})` : ""}${r.abbreviation}`;
        if (sideChain.length) text += `(${sideChain.join(",")})`;
        const opens = result.cyclisations.filter((c) => c.kind === "cyclo" && c.start === p).length;
        const closes = result.cyclisations.filter((c) => c.kind === "cyclo" && c.end === p).length;
        pieces.push(`${"cyclo[".repeat(opens)}${text}${"]".repeat(closes)}`);
      }
      let body = pieces.join("-");
      if (chain.nTerminus.text) body = `${chain.nTerminus.text}-${body}`;
      if (chain.cTerminus.text) body = `${body}-${chain.cTerminus.text}`;
      return chain.label ? `${chain.label}: ${body}` : body;
    })
    .join(" | ");
}
