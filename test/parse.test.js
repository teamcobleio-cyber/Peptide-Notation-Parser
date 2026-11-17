import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { parsePeptide, toOneLetter, toThreeLetter, composition, registerResidue, lookupResidue } from "../dist/index.js";

const corpus = JSON.parse(readFileSync(new URL("./fixtures/library-sequences.json", import.meta.url), "utf8"));

test("every library sequence parses to its published residue count", () => {
  const wrong = [];
  for (const [slug, entry] of Object.entries(corpus)) {
    if (entry.residueCount == null) continue;
    const result = parsePeptide(entry.sequence);
    if (result.residueCount !== entry.residueCount) wrong.push(`${slug}: ${result.residueCount} vs ${entry.residueCount}`);
  }
  assert.deepEqual(wrong, []);
});

test("no library sequence leaves an unrecognised token", () => {
  const left = [];
  for (const [slug, entry] of Object.entries(corpus)) {
    const result = parsePeptide(entry.sequence);
    if (result.unrecognised.length > 0) left.push(`${slug}: ${result.unrecognised.map((u) => u.token).join(", ")}`);
  }
  // Printed so a regression names the entries rather than just a count.
  if (left.length) console.log("unrecognised tokens:\n" + left.join("\n"));
  assert.deepEqual(left, []);
});

test("normalised three-letter output parses back to the same residues", () => {
  for (const [slug, entry] of Object.entries(corpus)) {
    const first = parsePeptide(entry.sequence);
    const second = parsePeptide(toThreeLetter(first));
    assert.equal(second.residueCount, first.residueCount, slug);
    assert.deepEqual(
      second.residues.map((r) => [r.abbreviation, r.stereo, r.modifications.join(";")]),
      first.residues.map((r) => [r.abbreviation, r.stereo, r.modifications.join(";")]),
      slug,
    );
  }
});

test("PT-141: acetyl cap, cyclo ring over residues 2 to 7, free acid", () => {
  const r = parsePeptide("Ac-Nle-cyclo[Asp-His-D-Phe-Arg-Trp-Lys]-OH");
  assert.equal(r.residueCount, 7);
  assert.equal(r.nTerminus.type, "acetyl");
  assert.equal(r.cTerminus.type, "acid");
  assert.deepEqual(r.cyclisations.map((c) => [c.kind, c.start, c.end]), [["cyclo", 2, 7]]);
  assert.equal(r.residues[0].abbreviation, "Nle");
  assert.equal(r.residues[0].standard, false);
  assert.equal(r.residues[3].stereo, "D");
  assert.equal(toOneLetter(r), null);
});

test("Ipamorelin: Aib is achiral, 2-Nal and Phe are D, lysine amide", () => {
  const r = parsePeptide("Aib-His-D-2-Nal-D-Phe-Lys-NH2");
  assert.deepEqual(r.residues.map((x) => x.abbreviation), ["Aib", "His", "2-Nal", "Phe", "Lys"]);
  assert.deepEqual(r.residues.map((x) => x.stereo), [null, "L", "D", "D", "L"]);
  assert.equal(r.cTerminus.type, "amide");
});

test("side-chain and backbone modifications are kept, not dropped", () => {
  const acyl = parsePeptide("Tyr-Aib-Glu-Lys(acyl)-NH2");
  assert.deepEqual(acyl.residues[3].modifications, ["side chain: acyl"]);
  const nMe = parsePeptide("Ser-(NMe)Asn-Phe");
  assert.deepEqual(nMe.residues[1].modifications, ["N-methyl"]);
  const aMe = parsePeptide("Glu-(alphaMe)Phe-Leu");
  assert.deepEqual(aMe.residues[1].modifications, ["alpha-methyl"]);
  const wholeParen = parsePeptide("Ser-Ile-(alpha-Me-Leu)-Leu");
  assert.equal(wholeParen.residues[2].abbreviation, "Me-Leu");
  assert.equal(wholeParen.residues[2].name.startsWith("alpha-Methylleucine"), true);
});

test("acyl caps that contain hyphens and digits", () => {
  const r = parsePeptide("(3E)-hex-3-enoyl-Tyr-Ala-Asp-NH2");
  assert.equal(r.nTerminus.type, "acyl");
  assert.equal(r.nTerminus.text, "(3E)-hex-3-enoyl");
  assert.equal(r.residueCount, 3);
  // A residue name ending in "yl" is not a cap.
  const hex = parsePeptide("His-D-2-methyl-Trp-Ala-Trp-D-Phe-Lys-NH2");
  assert.equal(hex.nTerminus.type, "unspecified");
  assert.equal(hex.residues[1].abbreviation, "2-Me-Trp");
  assert.equal(hex.residues[1].stereo, "D");
});

test("lower-case one-letter codes are D-residues and glycine stays achiral", () => {
  const r = parsePeptide("YaDAg");
  assert.deepEqual(r.residues.map((x) => x.stereo), ["L", "D", "L", "L", null]);
  const off = parsePeptide("YaDA", { lowercaseIsD: false });
  assert.equal(off.residues[1].stereo, "L");
});

test("one-letter run with a C-terminal amide and a disulfide note", () => {
  const r = parsePeptide("CYIQNCPLG-NH2 (cyclic 1-6 disulfide)");
  assert.equal(r.residueCount, 9);
  assert.equal(r.cTerminus.type, "amide");
  assert.deepEqual(r.cyclisations.map((c) => [c.kind, c.start, c.end]), [["disulfide", 1, 6]]);
  assert.deepEqual(r.warnings, []);
  assert.equal(toOneLetter(r), null);
  assert.equal(toOneLetter(r, { ignoreTermini: true }), "CYIQNCPLG");
});

test("a disulfide note on a non-cysteine position is reported", () => {
  const r = parsePeptide("AYIQNCPLG (disulfide 1-6)");
  assert.ok(r.warnings.some((w) => w.includes("position 1")));
});

test("PDB-style codes inside a one-letter run", () => {
  const r = parsePeptide("SYS(NLE)EH(DPN)RWGKPV(NH2)");
  assert.equal(r.residueCount, 13);
  assert.equal(r.residues[3].abbreviation, "Nle");
  assert.equal(r.residues[6].abbreviation, "Phe");
  assert.equal(r.residues[6].stereo, "D");
  assert.equal(r.cTerminus.type, "amide");
});

test("labelled chains are kept apart", () => {
  const r = parsePeptide("Alpha chain (3 residues): GHK | Beta chain (2 residues): AE");
  assert.deepEqual(r.chains.map((c) => [c.label, c.start, c.end]), [["Alpha chain", 1, 3], ["Beta chain", 4, 5]]);
  assert.equal(toOneLetter(r), "GHK/AE");
});

test("anything unreadable is returned, never dropped", () => {
  const r = parsePeptide("Ala-Foo-Gly");
  assert.equal(r.residueCount, 2);
  assert.deepEqual(r.unrecognised, [{ token: "Foo", chain: 0 }]);
  assert.equal(toOneLetter(r), null);
  const x = parsePeptide("GXK");
  assert.equal(x.residues[1].abbreviation, "Xaa");
  assert.ok(x.warnings.length > 0);
});

test("U is read as selenocysteine, with a warning", () => {
  const r = parsePeptide("UEQLERALNSS");
  assert.equal(r.residues[0].abbreviation, "Sec");
  assert.ok(r.warnings.some((w) => w.includes("pyroglutamate")));
  const pyro = parsePeptide("UEQLERALNSS", { u: "pyroglutamate" });
  assert.equal(pyro.residues[0].abbreviation, "pGlu");
  assert.deepEqual(pyro.warnings, []);
});

test("plain sequences round-trip to one-letter code", () => {
  assert.equal(toOneLetter(parsePeptide("GEPPPGKPADDAGLV")), "GEPPPGKPADDAGLV");
  assert.equal(toOneLetter(parsePeptide("Lys-Pro-Val")), "KPV");
  assert.equal(toOneLetter(parsePeptide("H-Ala-His-Lys-OH")), "AHK");
  assert.equal(toOneLetter(parsePeptide("Ac-Leu-Lys-Lys-Thr-Glu-Thr-Gln-OH")), null);
});

test("composition summary", () => {
  const c = composition(parsePeptide("D-Arg-Dmt-Lys-Phe-NH2"));
  assert.equal(c.total, 4);
  assert.equal(c.proteinogenic, 3);
  assert.equal(c.nonStandard, 1);
  assert.equal(c.dResidues, 1);
  assert.equal(c.byResidue.Dmt, 1);
});

test("the residue table is extensible", () => {
  assert.equal(lookupResidue("Hph"), undefined);
  registerResidue({ abbreviation: "Hph", oneLetter: null, name: "Homophenylalanine", proteinogenic: false, chiral: true, aliases: ["HPH"] });
  const r = parsePeptide("Ala-Hph-Gly");
  assert.equal(r.unrecognised.length, 0);
  assert.equal(r.residues[1].name, "Homophenylalanine");
});
