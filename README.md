# Peptide Notation Parser

> This library is maintained by **[Volta Peptides](https://voltapeptides.com/)**, a Canadian supplier of research peptides. Every worked example below is a compound in the [Volta Peptides catalogue](https://voltapeptides.com/catalog), linked to its product page.

Peptides are written in half a dozen notations, often in the same document: one-letter runs (`GEPPPGKPADDAGLV`), three-letter chains (`Lys-Pro-Val`), stereo prefixes (`D-Phe`), terminal caps (`Ac-`, `-NH2`), non-standard residues (`Aib`, `Nle`, `D-2-Nal`), side-chain conjugates (`Lys(acyl)`) and rings (`cyclo[...]`). This parser reads all of them into one structure: an ordered residue list, the two termini, ring and disulfide ranges, and every token it could not place.

**Nothing is dropped.** A token the parser does not recognise comes back in `unrecognised`; an ambiguity it resolved by convention comes back in `warnings`. A parser that skipped what it did not understand would report part of a molecule as the whole of it.

Dependency-free TypeScript, ESM, Node 18+. The test suite parses all 53 sequences in the Volta Peptides compound library, checks each against its published residue count, and requires zero unrecognised tokens.

## Install

```bash
npm install github:teamcobleio-cyber/Peptide-Notation-Parser
# or
npm install git+https://gitlab.com/sofiazugas/peptide-notation-parser.git
```

## Quick start

```js
import { parsePeptide, toOneLetter, toThreeLetter, composition } from "peptide-notation-parser";

const pt141 = parsePeptide("Ac-Nle-cyclo[Asp-His-D-Phe-Arg-Trp-Lys]-OH");

pt141.residueCount;        // 7
pt141.nTerminus;           // {"type":"acetyl","text":"Ac"}
pt141.cTerminus;           // {"type":"acid","text":"OH"}
pt141.cyclisations;        // [{"kind":"cyclo","start":2,"end":7}]
pt141.residues[0];         // {"abbreviation":"Nle","stereo":"L","standard":false} (and more fields)
pt141.residues[3].stereo;  // "D" (D-Phe)
toOneLetter(pt141);        // null, because Nle is not a standard residue
```

That is [PT-141](https://voltapeptides.com/product/pt-141-10mg); its [compound hub](https://voltapeptides.com/peptides/pt-141) and [library entry](https://voltapeptides.com/research/compounds/pt-141) explain the hormone-based numbering its literature uses.

### Worked examples from the catalogue

**[Ipamorelin](https://voltapeptides.com/product/ipamorelin-5mg)**: an achiral residue, a non-standard aromatic, two D-residues and an amidated C-terminus.

```js
const ipa = parsePeptide("Aib-His-D-2-Nal-D-Phe-Lys-NH2");
ipa.residues.map((r) => r.abbreviation); // ["Aib","His","2-Nal","Phe","Lys"]
ipa.residues.map((r) => r.stereo);       // [null,"L","D","D","L"]
ipa.cTerminus.type;                      // "amide"
```

**[Melanotan II](https://voltapeptides.com/product/melanotan-ii-10mg)** differs from PT-141 only at the C-terminus, and the parser says so:

```js
parsePeptide("Ac-Nle-cyclo[Asp-His-D-Phe-Arg-Trp-Lys]-NH2").cTerminus; // {"type":"amide","text":"NH2"}
parsePeptide("Ac-Nle-cyclo[Asp-His-D-Phe-Arg-Trp-Lys]-OH").cTerminus;  // {"type":"acid","text":"OH"}
```

**[SS-31](https://voltapeptides.com/product/ss-31-10mg)**: `D-Arg-Dmt-Lys-Phe-NH2`, where Dmt is 2',6'-dimethyltyrosine.

```js
composition(parsePeptide("D-Arg-Dmt-Lys-Phe-NH2"));
// {"total":4,"proteinogenic":3,"nonStandard":1,"dResidues":1,"modified":0,"placeholders":0,"unrecognised":0,"byResidue":{"Arg":1,"Dmt":1,"Lys":1,"Phe":1}}
```

**[TB-500](https://voltapeptides.com/product/tb-500-5mg)**: an acetylated heptapeptide, so it has a one-letter backbone but is not a plain one-letter peptide.

```js
const tb = parsePeptide("Ac-Leu-Lys-Lys-Thr-Glu-Thr-Gln-OH");
toOneLetter(tb);                          // null (capped N-terminus)
toOneLetter(tb, { ignoreTermini: true }); // "LKKTETQ"
```

**[Tesamorelin](https://voltapeptides.com/product/tesamorelin-10mg)**: a 44-residue chain behind an acyl cap that itself contains hyphens and digits.

```js
parsePeptide("(3E)-hex-3-enoyl-Tyr-Ala-Asp-Ala-Ile-Phe-...-Arg-Leu-NH2").nTerminus;
// {"type":"acyl","text":"(3E)-hex-3-enoyl"}
```

**[GLP-3 RT](https://voltapeptides.com/product/glp-3-rt-5mg)**: 39 residues with two Aib, an alpha-methyl leucine written as `(alpha-Me-Leu)`, and an acylated lysine.

```js
const rt = parsePeptide(sequence); // the full sequence is in test/fixtures
composition(rt).nonStandard;  // 3
composition(rt).modified;     // 1
rt.residues.filter((r) => !r.standard).map((r) => r.abbreviation + r.position);
// ["Aib2","Me-Leu13","Aib20"]
```

## What it reads

| Notation | Example | Read as |
|---|---|---|
| One-letter run | `GEPPPGKPADDAGLV` | 15 L-residues |
| Lower-case one-letter | `YaDAIF` | a is D-Ala (turn off with `{ lowercaseIsD: false }`) |
| Three-letter chain | `Lys-Pro-Val` | 3 residues |
| Stereo prefix | `D-Phe`, `L-Ala` | configuration set explicitly |
| N-terminal caps | `H-`, `Ac-`, `Hexanoyl-`, `(3E)-hex-3-enoyl-` | free, acetyl, acyl |
| C-terminal groups | `-OH`, `-NH2`, `-NHEt`, `-NH-(CH2)5-CONH2`, `-[...]` | acid, amide, substituted amide, group |
| Non-standard residues | `Aib`, `Nle`, `Orn`, `2-Nal`, `Dmt`, `Ac4c`, `gamma-Glu`, `pGlu`, `Sar`, `Hyp` and more | named, flagged non-standard |
| Backbone modifiers | `(alphaMe)Phe`, `(NMe)Asn`, `(alpha-Me-Leu)` | modification recorded |
| Side-chain modifications | `Lys(acyl)`, `Lys(Ac)`, `Lys(Nepsilon-3-maleimidopropionyl)` | modification recorded |
| Rings | `cyclo[Asp-...-Lys]` | cyclisation range |
| Disulfide notes | `(cyclic 1-6 disulfide)`, `(disulfide 1-6)` | disulfide range, checked against Cys positions |
| PDB component codes in a run | `SYS(NLE)EH(DPN)RWGKPV(NH2)` | Nle, D-Phe, amide |
| Labelled chains | `Alpha chain (92 residues): ...` then `Beta chain ...` | separate chains |

`U` is read as selenocysteine, as the IUPAC alphabet defines it, with a warning, because some sources use it for pyroglutamate; pass `{ u: "pyroglutamate" }` when yours does. `X`, `B` and `Z` are placeholders and come back as `Xaa` with a warning.

Add a residue the table does not know:

```js
import { registerResidue } from "peptide-notation-parser";
registerResidue({ abbreviation: "Hph", oneLetter: null, name: "Homophenylalanine", proteinogenic: false, chiral: true });
```

## No masses for non-standard residues

This library identifies residues; it does not weigh them. A mass for a non-standard residue needs a stated source and a stated form (residue, free acid or salt), and a guessed one is worse than none. For a sequence of the twenty standard residues, pass `toOneLetter(result)` to the companion **[Peptide Chemistry Toolkit](https://github.com/teamcobleio-cyber/Peptide-Chemistry-Toolkit)**, which computes average and monoisotopic mass, pI, net charge and fragment ions from sourced tables. For a modified peptide, start from its published formula: every Volta Peptides product page lists the molecular formula and weight.

## Catalogue peptides with modified notation

Each of these parses with zero unrecognised tokens. They are the sequences `toOneLetter` correctly refuses, because a one-letter string would hide a non-standard residue, a D-residue, a modification or a capped terminus.

| Compound | Residues | Non-standard | D-residues | Modified | Termini and rings | Volta Peptides |
|---|---|---|---|---|---|---|
| Adipotide | 25 | none | 14 | none | N: free or unstated; C: free or unstated | [Product](https://voltapeptides.com/product/adipotide-5mg) |
| ARA-290 (U read as pyroglutamate) | 11 | pGlu | 0 | none | N: free or unstated; C: free or unstated | [Product](https://voltapeptides.com/product/ara-290-10mg) |
| Cagrilintide | 37 | none | 0 | Lys1 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/cagrilintide-5mg) · [Hub](https://voltapeptides.com/peptides/cagrilintide) · [Library](https://voltapeptides.com/research/compounds/cagrilintide) |
| CJC-1295 with DAC | 30 | none | 1 | Lys30 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/cjc-1295-dac-5mg) · [Hub](https://voltapeptides.com/peptides/cjc-1295-dac) · [Library](https://voltapeptides.com/research/compounds/cjc-1295-dac) |
| CJC-1295 without DAC | 29 | none | 1 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/cjc-1295-no-dac-10mg) · [Hub](https://voltapeptides.com/peptides/cjc-1295-no-dac) · [Library](https://voltapeptides.com/research/compounds/cjc-1295-no-dac) |
| Dihexa | 2 | none | 0 | none | N: Hexanoyl; C: NH-(CH2)5-CONH2 | [Product](https://voltapeptides.com/product/dihexa-10mg) |
| Eloralintide | 37 | gGlu, Orn | 0 | Phe15, Asn22, Lys26 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/eloralintide-10mg) |
| FOXO4-DRI | 46 | none | 42 | none | N: free or unstated; C: free or unstated | [Product](https://voltapeptides.com/product/foxo4-dri-10mg) |
| GHRP-2 | 6 | 2-Nal | 3 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/ghrp-2-5mg) |
| GHRP-6 | 6 | none | 2 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/ghrp-6-5mg) |
| GLP-1 SM | 31 | Aib | 0 | Lys20 | N: free or unstated; C: free or unstated | [Product](https://voltapeptides.com/product/glp-1-sm-10mg) · [Hub](https://voltapeptides.com/peptides/glp-1-sm) |
| GLP-2 TZ | 39 | Aib | 0 | Lys20 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/glp-2-tz-10mg) · [Hub](https://voltapeptides.com/peptides/glp-2-tz) |
| GLP-3 RT | 39 | Aib, Me-Leu | 0 | Lys17 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/glp-3-rt-5mg) · [Hub](https://voltapeptides.com/peptides/glp-3-rt) |
| Glutathione | 3 | gGlu | 0 | none | N: free or unstated; C: free or unstated | [Product](https://voltapeptides.com/product/glutathione-1500mg) · [Hub](https://voltapeptides.com/peptides/glutathione) · [Library](https://voltapeptides.com/research/compounds/glutathione) |
| Gonadorelin | 10 | pGlu | 0 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/gonadorelin-acetate-5mg) |
| Hexarelin | 6 | 2-Me-Trp | 2 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/hexarelin-acetate-5mg) |
| Ipamorelin | 5 | Aib, 2-Nal | 2 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/ipamorelin-5mg) · [Hub](https://voltapeptides.com/peptides/ipamorelin) · [Library](https://voltapeptides.com/research/compounds/ipamorelin) |
| Mazdutide | 34 | Aib | 0 | Lys20 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/mazdutide-10mg) |
| Melanotan I | 13 | Nle | 1 | none | N: Ac; C: amide | [Product](https://voltapeptides.com/product/melanotan-1-10mg) |
| Melanotan II | 7 | Nle | 1 | none | N: Ac; C: amide; cyclo 2 to 7 | [Product](https://voltapeptides.com/product/melanotan-ii-10mg) |
| Oxytocin | 9 | none | 0 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/oxytocin-acetate-10mg) |
| P021 | 4 | none | 0 | none | N: Ac; C: [3-amino-adamantane-1-carboxamide] | [Product](https://voltapeptides.com/product/p21-5mg) |
| PT-141 | 7 | Nle | 1 | none | N: Ac; C: acid; cyclo 2 to 7 | [Product](https://voltapeptides.com/product/pt-141-10mg) · [Hub](https://voltapeptides.com/peptides/pt-141) · [Library](https://voltapeptides.com/research/compounds/pt-141) |
| Sermorelin | 29 | none | 0 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/sermorelin-5mg) · [Hub](https://voltapeptides.com/peptides/sermorelin) · [Library](https://voltapeptides.com/research/compounds/sermorelin) |
| SS-31 | 4 | Dmt | 1 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/ss-31-10mg) · [Hub](https://voltapeptides.com/peptides/ss-31) · [Library](https://voltapeptides.com/research/compounds/ss-31) |
| Survodutide | 29 | Ac4c | 0 | Lys24 | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/survodutide-10mg) |
| TB-500 | 7 | none | 0 | none | N: Ac; C: acid | [Product](https://voltapeptides.com/product/tb-500-5mg) · [Hub](https://voltapeptides.com/peptides/tb-500) · [Library](https://voltapeptides.com/research/compounds/tb-500) |
| Tesamorelin | 44 | none | 0 | none | N: (3E)-hex-3-enoyl; C: amide | [Product](https://voltapeptides.com/product/tesamorelin-10mg) · [Hub](https://voltapeptides.com/peptides/tesamorelin) · [Library](https://voltapeptides.com/research/compounds/tesamorelin) |
| Thymosin alpha-1 | 28 | none | 0 | none | N: Ac; C: free or unstated | [Product](https://voltapeptides.com/product/thymosin-alpha-1-10mg) · [Hub](https://voltapeptides.com/peptides/thymosin-alpha-1) · [Library](https://voltapeptides.com/research/compounds/thymosin-alpha-1) |
| Vasoactive intestinal peptide | 28 | none | 0 | none | N: free or unstated; C: amide | [Product](https://voltapeptides.com/product/vip-10mg) |

Browse by composition: [research peptides](https://voltapeptides.com/catalog?type=peptides), [Melanotan peptides](https://voltapeptides.com/catalog?type=melanotan-peptides), [peptide blends](https://voltapeptides.com/catalog?type=peptide-blends) and [bioregulators](https://voltapeptides.com/catalog?type=bioregulators).

## API

| Export | Returns |
|---|---|
| `parsePeptide(text, options?)` | `{ residues, residueCount, chains, nTerminus, cTerminus, cyclisations, unrecognised, warnings }` |
| `toOneLetter(result, { ignoreTermini? })` | One-letter string, or `null` unless every residue is a standard, unmodified L (or achiral) residue. Chains are joined by `/` |
| `toThreeLetter(result)` | Normalised three-letter notation that parses back to the same residues |
| `composition(result)` | Totals, non-standard, D-residues, modified, placeholders and counts by residue |
| `registerResidue(definition)` | Adds a residue or alias to the lookup table |
| `lookupResidue(text)`, `lookupOneLetter(letter)` | The residue definition, or `undefined` |
| `PROTEINOGENIC`, `NON_STANDARD` | The built-in residue tables |

Each residue carries `position`, `chain`, `abbreviation`, `oneLetter`, `name`, `stereo` (`"L"`, `"D"` or `null` for achiral), `stereoExplicit`, `standard`, `modifications` and the `source` text it was read from. An unprefixed chiral residue is reported as L, the usual convention; `stereoExplicit` tells you whether the notation actually said so.

## Development

```bash
npm install
npm run build
npm test
```

The fixture `test/fixtures/library-sequences.json` is a snapshot of the sequences in the [Volta Peptides compound library](https://voltapeptides.com/research/compounds).

## Licence

MIT. Built and maintained by [Volta Peptides](https://voltapeptides.com/). For research use only. Nothing here is medical advice.
